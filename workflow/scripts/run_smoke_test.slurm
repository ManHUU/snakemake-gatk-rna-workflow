#!/bin/bash
# =============================================================================
# run_smoke_test.slurm  (reviewer comments C-07 / C-08)
#
# HPC batch entry point for the miniature smoke test. Runs entirely on ONE
# compute node (do NOT run on a login node - alignment is too slow there):
#   1. builds the single-contig (chr22) fixture if it is missing;
#   2. runs the real workflow + containers on that fixture;
#   3. checks the expected outputs (non-empty final VCF, etc.).
#
# Submit with YOUR account/partition (and, if your cluster needs it, load the
# conda/apptainer modules first - see the MODULE LOADS section below):
#
#   sbatch --account=<your_account> --partition=<your_partition> \
#          workflow/scripts/run_smoke_test.slurm
#
# Optional overrides (export before sbatch, or edit the #SBATCH lines):
#   HPC_SCRATCH_DIR=/scratch/$USER/gatk   # container cache + temp on fast scratch
#   SMOKE_CONTIG=chr22                    # fixture contig
#   SMOKE_MAX_INPUT_PAIRS=5000000         # read pairs sliced per sample before align
# =============================================================================
#SBATCH --job-name=gatk_smoke
#SBATCH --output=logs/gatk_smoke_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00

set -euo pipefail

# Locate the repository root robustly and cd into it. This must work both when
# launched via `sbatch` (SLURM copies the script to a spool dir, so the real
# repo root is SLURM_SUBMIT_DIR = where you submitted from) AND when run directly
# as `bash workflow/scripts/run_smoke_test.slurm` inside an interactive srun
# session (where SLURM_SUBMIT_DIR points at wherever the srun was started, not
# the repo). We take the first candidate that actually contains a Snakefile:
#   1. SLURM_SUBMIT_DIR  (correct for sbatch)
#   2. $PWD              (correct when you cd'd into the repo first)
#   3. this script's own ../.. (correct when run in place)
_repo_root=""
for _cand in "${SLURM_SUBMIT_DIR:-}" "$PWD" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)"; do
    if [[ -n "$_cand" && -f "$_cand/Snakefile" ]]; then _repo_root="$_cand"; break; fi
done
if [[ -z "$_repo_root" ]]; then
    echo "ERROR: could not locate the repository root (no Snakefile in SLURM_SUBMIT_DIR, \$PWD, or the script's parent directory)." >&2
    exit 1
fi
cd "$_repo_root"
mkdir -p logs

# ── MODULE LOADS (uncomment/adapt for your cluster) ──────────────────────────
# Some clusters expose conda and apptainer only via environment modules:
#   module load 2024
#   module load Miniconda3
#   module load Apptainer      # or Singularity
# -----------------------------------------------------------------------------

# ── Conda environment ────────────────────────────────────────────────────────
# The batch job may not inherit your shell's PATH, so try common conda locations.
if ! command -v conda >/dev/null 2>&1; then
    for base in "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/mambaforge" /opt/conda; do
        if [[ -f "$base/etc/profile.d/conda.sh" ]]; then
            # shellcheck disable=SC1091
            source "$base/etc/profile.d/conda.sh"; break
        fi
    done
fi
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Load it via your cluster module system (edit the MODULE LOADS section above)." >&2
    exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake_env

# ── Container cache + temp (prefer fast scratch) ─────────────────────────────
SCRATCH="${HPC_SCRATCH_DIR:-$PWD/resources/containers}"
export APPTAINER_CACHEDIR="$SCRATCH/apptainer_cache"
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export APPTAINER_TMPDIR="$SCRATCH/apptainer_tmp"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

CONTIG="${SMOKE_CONTIG:-chr22}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

echo "============================================================"
echo "GATK RNA-seq smoke test on compute node $(hostname)"
echo "  contig=$CONTIG  threads=$THREADS  cache=$APPTAINER_CACHEDIR"
echo "============================================================"

# ── 1. Build the miniature fixture if missing ────────────────────────────────
if [[ ! -s resources/smoke/smoke.fa || -z "$(ls -A data/smoke 2>/dev/null || true)" ]]; then
    echo "[build] creating smoke fixture ..."
    SMOKE_THREADS="$THREADS" bash workflow/scripts/make_smoke_fixture.sh "$CONTIG"
else
    echo "[build] fixture already present - skipping."
fi

# ── 2. Run the smoke pipeline on THIS allocated node ─────────────────────────
# We run snakemake directly with the local profile so all rules execute on the
# allocated node (no nested SLURM submission). The full tutorial instead uses
# run.sh, which drives the per-rule SLURM executor.
BIND="$PWD"
[[ -n "${HPC_SCRATCH_DIR:-}" && -d "${HPC_SCRATCH_DIR}" ]] && BIND="${BIND},${HPC_SCRATCH_DIR}"
for p in /scratch-local /scratch-shared; do [[ -d "$p" ]] && BIND="${BIND},${p}"; done

echo "[run] starting workflow (wall-clock timing starts now) ..."
SECONDS=0
snakemake \
    --snakefile Snakefile \
    --profile workflow/profiles/local \
    --configfile config/smoke/config.yaml \
    --singularity-args "--bind ${BIND}" \
    --cores "$THREADS"
echo "[run] workflow finished in ${SECONDS}s of execution (excludes any first-time container pulls)."

# ── 3. Acceptance checks ─────────────────────────────────────────────────────
echo "[check] verifying expected outputs ..."
bash workflow/scripts/check_smoke_outputs.sh results_smoke config/smoke/expected_checkpoints.tsv
rc=$?

echo "============================================================"
if [[ $rc -eq 0 ]]; then
    echo "SMOKE TEST PASSED. Record the measured counts/timing above into"
    echo "config/smoke/expected_checkpoints.tsv and the manuscript [MEASURE] cells."
else
    echo "SMOKE TEST FAILED - see [FAIL] lines above and logs/ for details."
fi
exit $rc
