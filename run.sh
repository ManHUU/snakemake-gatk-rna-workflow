#!/bin/bash
#SBATCH --job-name=gatk_pipeline
#SBATCH --output=logs/gatk_pipeline_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=72:00:00

# =============================================================================
# Unified entrypoint for the GATK RNA-seq Snakemake pipeline.
#
#   bash   run.sh                # local; or self-submits via sbatch on HPC
#   bash   run.sh --dry-run      # preview the DAG (works on head nodes)
#   sbatch run.sh                # explicit SLURM submission
#
# Requirements on PATH:
#   snakemake               (install however you like: conda/mamba/micromamba/
#                            pip in a venv, or `module load snakemake` on HPC)
#   apptainer or singularity (every scientific tool runs from a container)
#
# ─ HPC settings (everything that varies from one cluster to the next) ────────
# These are ENVIRONMENT VARIABLES, consumed by this script / sbatch. Set them
# once per shell (or in your ~/.bashrc / job script). Local users set nothing.
#
#   REQUIRED on HPC:
#     export SLURM_ACCOUNT=<your_slurm_account>
#     export SLURM_PARTITION=<your_slurm_partition>
#
#   RECOMMENDED on HPC — your cluster's scratch path, where apptainer builds and
#   caches its (multi-GB) container images on fast disk. The path varies by site;
#   check your HPC's docs, or `env | grep -iE 'scratch|tmp'`. Common conventions:
#     export HPC_SCRATCH_DIR=/scratch-shared/$USER   # Snellius / SURF
#     export HPC_SCRATCH_DIR=$SCRATCH                # TACC and sites that set $SCRATCH
#     export HPC_SCRATCH_DIR=/scratch/$USER          # many university clusters
#   If set, it ALWAYS wins. If unset, the script auto-detects scratch from
#   $SCRATCH, /scratch-shared/$USER, /scratch/$USER (first writable, disk-backed
#   path wins) and only falls back to resources/containers/{tmp,cache} in the
#   repo if none exist. An inherited APPTAINER_TMPDIR that is tmpfs (e.g. /tmp on
#   Snellius compute nodes) or anywhere under /tmp is rejected — it would
#   re-trigger the image-build OOM bug. The resolved path is printed at startup.
#
#   OPTIONAL — make extra host dirs visible inside containers (comma-separated);
#   only needed if fastq_dir/output_dir point outside the repo (repo is auto-bound):
#     export EXTRA_BIND_PATHS=/scratch-shared,/tmp
#
# The one per-site value that is NOT an env var is `partition_max_runtime` in
# config/config.yaml (your partition's MaxTime, in minutes). It lives in config
# because Snakemake — not this script — consumes it. No other code edits needed.
# =============================================================================

set -euo pipefail

# ── Locate and enter the repository root ────────────────────────────────────
# Note: when SLURM submits a script it copies the file to a spool directory,
# so $0 inside the job no longer points to the repo.  We use SLURM_SUBMIT_DIR
# inside SLURM jobs and fall back to $0 for direct `bash` invocation.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/Snakefile" ]]; then
        cd "$SLURM_SUBMIT_DIR"
    else
        echo "ERROR: cannot locate Snakefile from"
        echo "       SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}."
        echo "       Please submit run.sh from the repository root."
        exit 1
    fi
else
    cd "$(cd "$(dirname "$0")" && pwd)"
fi

mkdir -p logs results resources/containers

# ── 1. Verify required tools on PATH ────────────────────────────────────────
# We do NOT require `conda` itself. The user may have activated their
# snakemake env with mamba/micromamba (or used pip in a venv, or
# `module load snakemake` on HPC) — `conda` may not exist in that shell.
# The only thing that matters is that `snakemake` is callable here.
if ! command -v snakemake >/dev/null 2>&1; then
    echo "ERROR: 'snakemake' not found on PATH."
    echo "Activate a snakemake environment first, for example:"
    echo "  conda activate snakemake_env       # if you used conda"
    echo "  mamba activate snakemake_env       # if you used mamba"
    echo "  micromamba activate snakemake_env  # if you used micromamba"
    echo "  module load snakemake              # on HPC with environment modules"
    echo "  source <your-venv>/bin/activate    # if you pip-installed snakemake"
    exit 1
fi
echo "Snakemake : $(snakemake --version)"

# ── Housekeeping fast-paths (no sbatch, no apptainer needed) ────────────────
# Operations like --unlock and --cleanup-metadata only touch .snakemake/
# bookkeeping files. Routing them through sbatch would waste a job slot
# (and queue wait) on a sub-second operation, so we short-circuit here.
# Typical use: a previous run was killed un-gracefully (SLURM timeout, OOM,
# scancel, node failure) and left a stale lock. Running `bash run.sh --unlock`
# on a head node clears it in seconds; the next `bash run.sh` then resumes
# from the first missing output.
for arg in "$@"; do
    case "$arg" in
        --unlock|--cleanup-metadata|--report|--report=*)
            # --report builds report.html (DAG + runtimes + provenance) from the
            # .snakemake metadata of a completed run; like --unlock it only reads
            # bookkeeping, so it runs locally with no sbatch/apptainer needed.
            # Example: bash run.sh --report report.html
            echo "Housekeeping: running snakemake $arg locally (no sbatch)."
            exec snakemake --snakefile Snakefile "$@"
            ;;
    esac
done

if command -v apptainer >/dev/null 2>&1; then
    echo "Apptainer : $(apptainer --version | head -1)"
elif command -v singularity >/dev/null 2>&1; then
    echo "Singularity : $(singularity --version | head -1)"
else
    echo "ERROR: neither apptainer nor singularity found on PATH."
    echo "       Every scientific tool in this pipeline runs from a container,"
    echo "       so one of them is required. Install apptainer (preferred) or"
    echo "       singularity, or 'module load apptainer' on HPC."
    exit 1
fi

# ── Apptainer scratch space (must be disk-backed) ───────────────────────────
# On some HPCs (e.g. Snellius) /tmp is tmpfs (RAM-backed), so apptainer's
# default build location eats the job's memory budget and large images
# (GATK, STAR) OOM-kill the driver during `mksquashfs`.
#
# Resolution order (first match wins). Explicit user choice beats inherited
# environment, which is the most-suspect input:
#   1. HPC_SCRATCH_DIR set by the user → ${HPC_SCRATCH_DIR}/apptainer-{tmp,cache}.
#      This is the explicit override and ALWAYS wins, even over an inherited
#      APPTAINER_TMPDIR.
#   2. A *safe* inherited APPTAINER_TMPDIR — used as-is, UNLESS it is tmpfs
#      (e.g. /tmp on Snellius compute nodes) or anywhere under /tmp, in which case
#      we discard it with a warning. /tmp is small, shared, often auto-purged, and
#      is tmpfs on compute nodes, so it would re-trigger the OOM bug during
#      mksquashfs.
#   3. Auto-detect from common HPC conventions: $SCRATCH, /scratch-shared/$USER,
#      /scratch/$USER. First path that is writable AND not tmpfs wins.
#   4. In-repo fallback: $(pwd)/resources/containers/{tmp,cache}.

# Returns 0 if $1 is a directory on a tmpfs filesystem.
is_tmpfs() {
    [[ -d "$1" ]] || return 1
    local t
    t=$(stat -f -c %T "$1" 2>/dev/null || echo "")
    [[ "$t" == "tmpfs" ]]
}

# Returns 0 if $1 is /tmp or a path under it. /tmp is never a safe place to build
# multi-GB images: small, shared, often auto-purged, and tmpfs on many compute nodes.
is_under_tmp() {
    case "$1" in
        /tmp|/tmp/*) return 0 ;;
        *)           return 1 ;;
    esac
}

# Step 1: an explicit HPC_SCRATCH_DIR is the user's override and always wins.
# Drop any inherited APPTAINER_TMPDIR so Step 4 derives it from scratch instead.
if [[ -n "${HPC_SCRATCH_DIR:-}" && -n "${APPTAINER_TMPDIR:-}" ]]; then
    echo "Note: HPC_SCRATCH_DIR is set; ignoring inherited APPTAINER_TMPDIR=$APPTAINER_TMPDIR."
    unset APPTAINER_TMPDIR APPTAINER_CACHEDIR
fi

# Step 2: reject an inherited APPTAINER_TMPDIR that is unsafe for image builds —
# tmpfs (RAM-backed, OOMs the driver) or anything under /tmp (small/shared/purged,
# and tmpfs on compute nodes). Either way, fall through to auto-detect below.
if [[ -n "${APPTAINER_TMPDIR:-}" ]] && { is_tmpfs "$APPTAINER_TMPDIR" || is_under_tmp "$APPTAINER_TMPDIR"; }; then
    echo "WARNING: ignoring APPTAINER_TMPDIR=$APPTAINER_TMPDIR"
    echo "         (tmpfs or under /tmp — unsafe for multi-GB image builds; would"
    echo "         risk OOM or a full /tmp). Falling back to scratch auto-detect."
    unset APPTAINER_TMPDIR APPTAINER_CACHEDIR
fi

# Step 3: if HPC_SCRATCH_DIR and APPTAINER_TMPDIR are both unset, probe for
# scratch using common HPC conventions.
if [[ -z "${APPTAINER_TMPDIR:-}" && -z "${HPC_SCRATCH_DIR:-}" ]]; then
    for _cand in "${SCRATCH:-}" "/scratch-shared/${USER:-}" "/scratch/${USER:-}"; do
        [[ -z "$_cand" || "$_cand" == "/scratch-shared/" || "$_cand" == "/scratch/" ]] && continue
        if mkdir -p "$_cand" 2>/dev/null && [[ -w "$_cand" ]] && ! is_tmpfs "$_cand"; then
            HPC_SCRATCH_DIR="$_cand"
            echo "Auto-detected scratch: $HPC_SCRATCH_DIR"
            echo "         (Override with: export HPC_SCRATCH_DIR=<path>)"
            break
        fi
    done
    unset _cand
fi

# Step 4: derive final APPTAINER_TMPDIR / APPTAINER_CACHEDIR.
if [[ -n "${HPC_SCRATCH_DIR:-}" ]]; then
    : "${APPTAINER_TMPDIR:=${HPC_SCRATCH_DIR}/apptainer-tmp}"
    : "${APPTAINER_CACHEDIR:=${HPC_SCRATCH_DIR}/apptainer-cache}"
else
    : "${APPTAINER_TMPDIR:=$(pwd)/resources/containers/tmp}"
    : "${APPTAINER_CACHEDIR:=$(pwd)/resources/containers/cache}"
fi
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
export HPC_SCRATCH_DIR APPTAINER_TMPDIR APPTAINER_CACHEDIR
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR" SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
echo "Apptainer scratch : $APPTAINER_TMPDIR"

# ── 2. Build the singularity --bind argument ────────────────────────────────
SING_BIND="$(pwd)"
if [[ -n "${EXTRA_BIND_PATHS:-}" ]]; then
    SING_BIND="${SING_BIND},${EXTRA_BIND_PATHS}"
fi
# Bind node-local and shared scratch into containers so tools that write large
# intermediate files (notably STAR's _STARtmp/ during BAM sort, ~20-40 GB) can
# use Snakemake's `tmpdir` resource. On SLURM Snakemake auto-resolves tmpdir to
# /scratch-local/<user>.<jobid>, which keeps scratch off the project filesystem.
for _scratch in /scratch-local /scratch-shared; do
    [[ -d "$_scratch" ]] && SING_BIND="${SING_BIND},${_scratch}"
done
unset _scratch

# ── 3. Decide execution mode ────────────────────────────────────────────────
DRY=0
for arg in "$@"; do
    [[ "$arg" == "--dry-run" || "$arg" == "-n" ]] && DRY=1
done

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    MODE="slurm-job"
elif command -v sbatch >/dev/null 2>&1; then
    # On a head node with sbatch available.
    # Real runs self-submit via sbatch; --dry-run plans locally.
    if [[ "$DRY" == "1" ]]; then
        MODE="slurm-plan"
    else
        MODE="slurm-submit"
    fi
else
    MODE="local"
fi

# ── 4. Run ──────────────────────────────────────────────────────────────────
case "$MODE" in
slurm-submit)
    if [[ -z "${SLURM_ACCOUNT:-}" || -z "${SLURM_PARTITION:-}" ]]; then
        echo "ERROR: please export SLURM_ACCOUNT and SLURM_PARTITION before running:"
        echo "  export SLURM_ACCOUNT=<your_account>"
        echo "  export SLURM_PARTITION=<your_partition>"
        exit 1
    fi
    if [[ -z "${HPC_SCRATCH_DIR:-}" ]]; then
        echo "WARNING: HPC_SCRATCH_DIR is not set; apptainer will build images"
        echo "         under $(pwd)/resources/containers. This works, but on HPC"
        echo "         the recommended setup is to point it at fast scratch, e.g."
        echo "           export HPC_SCRATCH_DIR=/scratch-shared/\$USER"
        echo "         See the header of run.sh for site-specific examples."
    fi
    echo "Self-submitting via sbatch (account=$SLURM_ACCOUNT, partition=$SLURM_PARTITION)"
    exec sbatch \
        --chdir="$(pwd)" \
        --account="$SLURM_ACCOUNT" \
        --partition="$SLURM_PARTITION" \
        --export=ALL,SLURM_ACCOUNT,SLURM_PARTITION,EXTRA_BIND_PATHS,HPC_SCRATCH_DIR,APPTAINER_TMPDIR,APPTAINER_CACHEDIR \
        "$0" "$@"
    ;;
slurm-job|slurm-plan)
    if [[ -z "${SLURM_ACCOUNT:-}" || -z "${SLURM_PARTITION:-}" ]]; then
        echo "ERROR: SLURM_ACCOUNT and SLURM_PARTITION must be set in the environment."
        exit 1
    fi
    snakemake \
        --snakefile Snakefile \
        --profile workflow/profiles/slurm \
        --singularity-args "--bind ${SING_BIND}" \
        --default-resources \
            "slurm_account=${SLURM_ACCOUNT}" \
            "slurm_partition=${SLURM_PARTITION}" \
            "runtime=600" "mem_mb=8000" \
        "$@"
    ;;
local)
    snakemake \
        --snakefile Snakefile \
        --profile workflow/profiles/local \
        --singularity-args "--bind ${SING_BIND}" \
        "$@"
    ;;
esac
