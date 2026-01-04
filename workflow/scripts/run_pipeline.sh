#!/usr/bin/bash

#!/bin/bash
#SBATCH --job-name=gatk_snakmake
#SBATCH --output=logs/gatk_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=genoa

# Load modules if required by your HPC
# module load singularity 

# --- 1. PREPARE ENVIRONMENT ---
# Load modules if required by your HPC
# module load singularity 


# A. Source the conda.sh script (Found using your $CONDA_EXE path)
# This enables the 'conda' command inside the batch job
source /home/mhu/miniconda3/etc/profile.d/conda.sh

# B. Activate your specific environment
conda activate snake310    #Note: Replace 'snake310' with your own environment for snakemake

# C. Check versions (Debugging)
# This writes to the log so you know it worked
echo "Active Conda Env: $CONDA_DEFAULT_ENV"
echo "Snakemake Version: $(snakemake --version)"
echo "Apptainer Version: $(apptainer --version)"


# Ensure logs directory exists
mkdir -p logs

# Execute Snakemake
# --use-conda: uses your yaml files in workflow/envs
# --use-singularity: pulls docker:// images and converts them to .sif automatically
snakemake \
    --snakefile Snakefile \
    --profile workflow/profiles/slurm \
    --use-conda \
    --use-singularity \
    --singularity-args "--bind /projects/prjs1616/gatk_paper_demo/root" \
    --latency-wait 60 \
    --rerun-incomplete



