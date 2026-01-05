#!/bin/bash
#SBATCH -A ausei14399 # Update it to your actual HPC account
#SBATCH -p genoa # Update it to your actual HPC partition name
#SBATCH --job-name=download_hg38
#SBATCH --output=logs/download_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00  # Increased time limit as dbSNP is large (~9GB)

# Move to your database directory
# Ensure the directory exists to avoid errors
mkdir -p /resources
cd /resources  # Update it to your actual path
echo "Starting download at $(date)"

# --- 1. Download GTF ---
echo "Downloading Gencode GTF..."
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz

# --- 2. Download FASTA ---
echo "Downloading Gencode FASTA..."
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz

# --- 3. Download GATK Resource Bundle (Mills & dbSNP) ---
# Source: Broad Institute Google Cloud Public Data (hg38/v0)
echo "Downloading GATK Resources..."

# Mills & 1000G Gold Standard Indels
# Note: Keeping this compressed (.vcf.gz) because the .tbi index requires it.
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# dbSNP 138
# Note: This file is uncompressed in the Broad bundle (~9GB).
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -c https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

echo "Downloads finished. Starting decompression of Gencode files..."

# --- 4. Unzip Gencode files ---
# (dbSNP is already unzipped, and Mills must stay zipped for the .tbi index to work)
gunzip -f gencode.v44.primary_assembly.annotation.gtf.gz
gunzip -f GRCh38.primary_assembly.genome.fa.gz

echo "Process complete at $(date)"
