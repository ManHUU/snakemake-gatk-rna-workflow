#!/bin/bash
# =============================================================================
# check_prerequisites.sh   (reviewer comment C-04)
#
# Inspects HOST readiness for the GATK RNA-seq workflow and reports PASS / WARN
# / FAIL for each check. It only inspects and reports — it installs nothing and
# runs no pipeline step. (The separate smoke test, workflow/scripts/run_smoke_test.sh,
# demonstrates real container/workflow EXECUTION on a miniature fixture.)
#
# Checks:
#   - Operating system is Linux (the scientific tools run inside Linux
#     Apptainer/Singularity containers).
#   - Snakemake on PATH, OR an environment manager (conda/mamba/micromamba)
#     that can create snakemake_env from workflow/envs/snakemake.yaml — either
#     one is sufficient (if Snakemake is already available you do not need a
#     manager, e.g. on clusters that provide `module load snakemake`).
#   - git, and an Apptainer OR Singularity runtime.
#   - Readable configuration (config/config.yaml) and sample sheet (units.tsv).
#   - Writable output_dir and a usable temporary directory.
#   - Available disk space for references + intermediates.
#   - SLURM: if a scheduler is present, that squeue/sinfo work and that
#     SLURM_ACCOUNT / SLURM_PARTITION are exported (required for HPC runs;
#     not needed for a local workstation run — use FORCE_LOCAL=1).
#
# Usage:
#   bash workflow/scripts/check_prerequisites.sh
#
# Exit code:
#   0 — no FAILs (WARNs may be present; read them)
#   1 — one or more mandatory requirements failed
# =============================================================================

set -uo pipefail

fail=0
warn=0

pass_msg() { printf "[PASS]  %s\n" "$1"; }
warn_msg() { printf "[WARN]  %s\n" "$1"; warn=$((warn + 1)); }
fail_msg() { printf "[FAIL]  %s\n" "$1"; fail=$((fail + 1)); }
info_msg() { printf "[INFO]  %s\n" "$1"; }

# Resolve a few config values (best-effort; defaults if parsing is unavailable).
CFG="config/config.yaml"
read_cfg() {  # read_cfg <dotted.key> <default>
    local key="$1" default="$2"
    if command -v python3 >/dev/null 2>&1 && [ -r "$CFG" ]; then
        python3 - "$key" "$default" <<'PY' 2>/dev/null || echo "$default"
import sys, yaml
key, default = sys.argv[1], sys.argv[2]
try:
    cfg = yaml.safe_load(open("config/config.yaml")) or {}
    cur = cfg
    for part in key.split("."):
        cur = cur[part]
    print(cur if cur not in (None, "") else default)
except Exception:
    print(default)
PY
    else
        echo "$default"
    fi
}

echo "GATK RNA-seq workflow — host prerequisite check"
echo "================================================"

# ── 1. Operating system ──────────────────────────────────────────────────────
os="$(uname -s 2>/dev/null || echo unknown)"
if [ "$os" = "Linux" ]; then
    pass_msg "Operating system: Linux ($(uname -r))"
else
    fail_msg "Operating system is '$os'; Linux is required (containers run on Linux)."
fi

# ── 2. Snakemake OR an environment manager that can install it ─────────────────
#    Accept either: if snakemake is already on PATH (e.g. via `module load
#    snakemake`) no manager is needed; otherwise conda/mamba/micromamba can
#    create snakemake_env from workflow/envs/snakemake.yaml.
if command -v snakemake >/dev/null 2>&1; then
    pass_msg "Snakemake: $(snakemake --version 2>&1 | head -1)"
else
    mgr=""
    for m in conda mamba micromamba; do
        if command -v "$m" >/dev/null 2>&1; then mgr="$m"; break; fi
    done
    if [ -n "$mgr" ]; then
        pass_msg "Environment manager: $mgr ($("$mgr" --version 2>&1 | head -1)) — use it to create snakemake_env"
    else
        fail_msg "Snakemake not on PATH, and no conda/mamba/micromamba to create snakemake_env."
    fi
fi

# ── 3. git ────────────────────────────────────────────────────────────────────
if command -v git >/dev/null 2>&1; then
    pass_msg "git: $(git --version 2>&1 | head -1)"
else
    fail_msg "git not found on PATH (needed to clone the repository)."
fi

# ── 4. Apptainer / Singularity ────────────────────────────────────────────────
if command -v apptainer >/dev/null 2>&1; then
    pass_msg "Container runtime: $(apptainer --version 2>&1 | head -1)"
elif command -v singularity >/dev/null 2>&1; then
    pass_msg "Container runtime: $(singularity --version 2>&1 | head -1)"
else
    fail_msg "Neither apptainer nor singularity found (required to run all tools)."
fi

# ── 5. Readable configuration and sample sheet ────────────────────────────────
if [ -r "$CFG" ]; then
    pass_msg "Configuration readable: $CFG"
else
    fail_msg "Cannot read $CFG (run this from the repository root)."
fi
# The workflow reads samples from the sample sheet (units.tsv), NOT by globbing
# fastq_dir; fastq_dir is used only by make_units.sh when you build the sheet.
UNITS="$(read_cfg units config/units.tsv)"
if [ -r "$UNITS" ]; then
    n=$(( $(wc -l < "$UNITS" 2>/dev/null || echo 1) - 1 ))
    pass_msg "Sample sheet present: $UNITS (${n} sample row(s))"
else
    warn_msg "Sample sheet $UNITS not found yet — build it with 'bash workflow/scripts/make_units.sh <fastq_dir> <paired|single>'."
fi

# ── 6. Writable output and temporary paths ────────────────────────────────────
OUTPUT_DIR="$(read_cfg output_dir results)"
out_parent="$(dirname "$OUTPUT_DIR")"
[ -d "$OUTPUT_DIR" ] && probe="$OUTPUT_DIR" || probe="$out_parent"
if [ -w "$probe" ] 2>/dev/null; then
    pass_msg "Output location writable: $probe"
else
    fail_msg "Output location not writable: $probe (output_dir=$OUTPUT_DIR)."
fi
tmp="${TMPDIR:-/tmp}"
if [ -d "$tmp" ] && [ -w "$tmp" ]; then
    pass_msg "Temporary directory writable: $tmp"
else
    warn_msg "Temporary directory '$tmp' not writable; set TMPDIR or HPC_SCRATCH_DIR."
fi

# Node-local scratch for heavy temp (STAR BAM sort). Routing heavy temp to a
# node-local SSD avoids the mid-sort truncation seen when STAR's intermediates
# land on a quota-limited project filesystem (see README "Disk space and runtime
# expectations"). On HPC this scratch is usually allocated PER JOB, so it may
# only be visible from inside a compute allocation — absence on a login node is
# normal, hence INFO rather than FAIL.
if [ -n "${HPC_SCRATCH_DIR:-}" ]; then
    if [ -d "$HPC_SCRATCH_DIR" ] && [ -w "$HPC_SCRATCH_DIR" ]; then
        pass_msg "Node-local scratch (HPC_SCRATCH_DIR) writable: $HPC_SCRATCH_DIR — heavy temp will route here."
    else
        warn_msg "HPC_SCRATCH_DIR='$HPC_SCRATCH_DIR' is set but not a writable directory."
    fi
else
    scratch_found=""
    for s in /scratch-local /scratch-shared; do
        if [ -d "$s" ] && [ -w "$s" ] 2>/dev/null; then scratch_found="$s"; break; fi
    done
    if [ -n "$scratch_found" ]; then
        info_msg "Node-local scratch present: $scratch_found — run.sh binds it so STAR's heavy temp routes there."
    else
        info_msg "No node-local scratch (/scratch-local, /scratch-shared) or HPC_SCRATCH_DIR detected here."
        info_msg "        On HPC this is usually allocated per job and appears only inside a compute allocation."
        info_msg "        If STAR's heavy temp lands on a quota-limited project filesystem it can truncate the"
        info_msg "        BAM sort; set HPC_SCRATCH_DIR to a fast node-local path, or ensure /scratch-local exists."
    fi
fi

# ── 7. Disk space (references ~40-50 GB + intermediates) ──────────────────────
avail_gb() { df -PBG "$1" 2>/dev/null | awk 'NR==2{gsub("G","",$4); print $4}'; }
space="$(avail_gb "${OUTPUT_DIR:-.}")"
[ -z "$space" ] && space="$(avail_gb .)"
if [ -n "$space" ]; then
    if [ "$space" -ge 100 ]; then
        pass_msg "Free disk at run location: ${space} GB"
    elif [ "$space" -ge 50 ]; then
        warn_msg "Free disk at run location: ${space} GB (>=100 GB recommended; refs ~40-50 GB + intermediates)."
    else
        fail_msg "Free disk at run location: ${space} GB (too low; need room for ~40-50 GB refs + intermediates)."
    fi
else
    warn_msg "Could not determine free disk space at the run location."
fi

# ── 8. Outbound HTTPS (reference downloads + container pulls) ──────────────────
if command -v curl >/dev/null 2>&1; then
    if curl --silent --fail --max-time 10 --head https://github.com >/dev/null 2>&1; then
        pass_msg "Outbound HTTPS reachable (https://github.com)."
    else
        warn_msg "Could not reach https://github.com in 10s; run downloads from a login/data-transfer node."
    fi
else
    warn_msg "curl not available; skipping outbound HTTPS probe."
fi

# ── 9. SLURM (only relevant for HPC runs) ─────────────────────────────────────
if command -v sbatch >/dev/null 2>&1; then
    if command -v squeue >/dev/null 2>&1 && command -v sinfo >/dev/null 2>&1; then
        pass_msg "SLURM client present (sbatch, squeue, sinfo)."
    else
        warn_msg "sbatch found but squeue/sinfo missing; SLURM client looks incomplete."
    fi
    if [ -n "${SLURM_ACCOUNT:-}" ] && [ -n "${SLURM_PARTITION:-}" ]; then
        pass_msg "SLURM_ACCOUNT and SLURM_PARTITION are set (account=$SLURM_ACCOUNT, partition=$SLURM_PARTITION)."
        if command -v sinfo >/dev/null 2>&1; then
            if sinfo -h -p "$SLURM_PARTITION" >/dev/null 2>&1 \
               && [ -n "$(sinfo -h -p "$SLURM_PARTITION" 2>/dev/null)" ]; then
                pass_msg "Partition '$SLURM_PARTITION' is known to the scheduler."
            else
                warn_msg "Partition '$SLURM_PARTITION' not listed by sinfo; check the name."
            fi
        fi
    else
        warn_msg "SLURM detected but SLURM_ACCOUNT/SLURM_PARTITION are not exported (required for HPC runs; export them before run.sh, or use FORCE_LOCAL=1 for a local run)."
    fi
else
    info_msg "No SLURM client (sbatch) found — this host will run the workflow locally."
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo "================================================"
if [ "$fail" -eq 0 ]; then
    echo "Result: PASS (${warn} warning(s)). Host prerequisites satisfied."
    if ! command -v snakemake >/dev/null 2>&1; then
        echo "Next: conda env create -f workflow/envs/snakemake.yaml && conda activate snakemake_env"
    fi
    exit 0
else
    echo "Result: FAIL — ${fail} mandatory requirement(s) missing, ${warn} warning(s)."
    echo "Resolve the [FAIL] items above. System-level tools (conda, git, apptainer)"
    echo "may be installed but hidden behind your cluster's module system — try e.g."
    echo "'module load apptainer' / 'module load Miniconda3' and re-run this check."
    exit 1
fi
