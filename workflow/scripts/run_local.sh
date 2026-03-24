#!/bin/bash
# ==============================================================================
# run_local.sh
# Run the GATK RNA-seq variant calling pipeline locally (no SLURM).
# Uses Apptainer (Singularity) containers and the local Snakemake profile.
#
# Usage:
#   bash workflow/scripts/run_local.sh            # full run
#   bash workflow/scripts/run_local.sh --dry-run  # preview steps only
# ==============================================================================

set -euo pipefail

# ── 1. Activate conda environment ────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake_env

echo "Snakemake : $(snakemake --version)"
echo "Apptainer : $(apptainer --version)"

# ── 2. Create required directories ───────────────────────────────────────────
mkdir -p logs results resources/containers

# ── 3. Run Snakemake ─────────────────────────────────────────────────────────
snakemake \
    --snakefile Snakefile \
    --profile workflow/profiles/local \
    "$@"
