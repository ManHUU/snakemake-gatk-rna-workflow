#!/bin/bash
#SBATCH --job-name=download_GSE256519
#SBATCH --output=logs/download_GSE256519_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --account=gusr114115
#SBATCH --partition=genoa

set -euo pipefail

mkdir -p logs

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate sra_tools

export THREADS="${SLURM_CPUS_PER_TASK:-8}"

bash workflow/scripts/download_GSE256519.sh
