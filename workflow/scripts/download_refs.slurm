#!/bin/bash
#SBATCH --job-name=download_refs
#SBATCH --output=logs/download_refs_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00

# =============================================================================
# Reference download — auto-detecting wrapper around download_refs.sh.
#
#   bash   workflow/scripts/download_refs.slurm    # local OR self-submits on HPC
#   sbatch workflow/scripts/download_refs.slurm    # explicit SLURM submission
#
# HPC users only need to set (same env vars used by run.sh):
#   export SLURM_ACCOUNT=<your_slurm_account>
#   export SLURM_PARTITION=<your_slurm_partition>
#
# Local users need set nothing.
# =============================================================================

set -euo pipefail

# ── Locate and enter the repository root ────────────────────────────────────
# Note: when SLURM submits a script it copies the file to a spool directory
# (e.g. /var/spool/slurm/slurmd/job<id>/), so $0 inside the job no longer
# points to workflow/scripts/.  We use SLURM_SUBMIT_DIR inside SLURM jobs and
# fall back to $0 for direct `bash` invocation.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/workflow/scripts/download_refs.sh" ]]; then
        cd "$SLURM_SUBMIT_DIR"
    else
        echo "ERROR: cannot locate workflow/scripts/download_refs.sh from"
        echo "       SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-<unset>}."
        echo "       Please submit this script from the repository root."
        exit 1
    fi
else
    cd "$(cd "$(dirname "$0")/../.." && pwd)"
fi

mkdir -p logs

# ── Decide execution mode ───────────────────────────────────────────────────
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    MODE="slurm-job"
elif command -v sbatch >/dev/null 2>&1; then
    MODE="slurm-submit"
else
    MODE="local"
fi

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
        --export=ALL,SLURM_ACCOUNT,SLURM_PARTITION \
        "$0" "$@"
    ;;
slurm-job|local)
    if command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate snakemake_env 2>/dev/null \
            || conda activate vcf_annotation 2>/dev/null \
            || true
    fi
    bash workflow/scripts/download_refs.sh
    ;;
esac
