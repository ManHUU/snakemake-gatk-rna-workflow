#!/bin/bash
#SBATCH --job-name=gatk_pipeline
#SBATCH --output=logs/gatk_pipeline_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00

# =============================================================================
# Unified entrypoint for the GATK RNA-seq Snakemake pipeline.
#
#   bash   run.sh                # local; or self-submits via sbatch on HPC
#   bash   run.sh --dry-run      # preview the DAG (works on head nodes)
#   sbatch run.sh                # explicit SLURM submission
#
# HPC users only need to set:
#   export SLURM_ACCOUNT=<your_slurm_account>
#   export SLURM_PARTITION=<your_slurm_partition>
#
# Optional: bind extra paths into Singularity containers (comma-separated):
#   export EXTRA_BIND_PATHS=/scratch-shared,/tmp
#
# Local users need set nothing — the script binds $(pwd) automatically.
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

# ── 1. Activate the snakemake conda env ─────────────────────────────────────
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: 'conda' not found. Install Miniconda or Mambaforge first."
    exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda activate snakemake_env 2>/dev/null; then
    echo "ERROR: conda env 'snakemake_env' not found."
    echo "Create it once with:"
    echo "  conda env create -f workflow/envs/snakemake.yaml"
    exit 1
fi
echo "Snakemake : $(snakemake --version)"
if command -v apptainer >/dev/null 2>&1; then
    echo "Apptainer : $(apptainer --version | head -1)"
elif command -v singularity >/dev/null 2>&1; then
    echo "Singularity : $(singularity --version | head -1)"
else
    echo "WARNING: neither apptainer nor singularity found on PATH."
fi

# ── 2. Build the singularity --bind argument ────────────────────────────────
SING_BIND="$(pwd)"
if [[ -n "${EXTRA_BIND_PATHS:-}" ]]; then
    SING_BIND="${SING_BIND},${EXTRA_BIND_PATHS}"
fi

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
    echo "Self-submitting via sbatch (account=$SLURM_ACCOUNT, partition=$SLURM_PARTITION)"
    exec sbatch \
        --chdir="$(pwd)" \
        --account="$SLURM_ACCOUNT" \
        --partition="$SLURM_PARTITION" \
        --export=ALL,SLURM_ACCOUNT,SLURM_PARTITION,EXTRA_BIND_PATHS \
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
