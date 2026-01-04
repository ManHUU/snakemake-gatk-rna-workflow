#!/bin/bash
#SBATCH -A ausei14399 # Update it to you actual HPC account
#SBATCH -p genoa # Update it to you actual HPC partition name
#SBATCH --job-name=download_hg38
#SBATCH --output=logs/download_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00

# Move to your database directory
cd /resources  # Update it to you actual path
echo "Starting download at $(date)"

# 1. Download GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz

# 2. Download FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

echo "Downloads finished. Starting decompression..."

# 3. Unzip them (gunzip is faster than downloading but still takes a few minutes)
gunzip gencode.v44.primary_assembly.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

echo "Process complete at $(date)"
