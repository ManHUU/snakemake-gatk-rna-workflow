#!/bin/bash
#SBATCH -A ausei14399 # Update it to your actual HPC account
#SBATCH -p genoa # Update it to your actual HPC partition name
#SBATCH --job-name=download_hg38
#SBATCH --output=logs/download_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00

# =============================================================================
# download_refs.slurm
#
# Downloads ALL reference databases required by the pipeline:
#   1. Gencode v44 GTF annotation
#   2. GRCh38 primary assembly FASTA
#   3. GATK resource bundle (Mills indels + dbSNP138)
#   4. REDIportal v2.0 RNA editing database (hg38)
#   5. SnpEff hg38 annotation database
#
# Usage (HPC):  sbatch workflow/scripts/download_refs.slurm
# Usage (local): bash workflow/scripts/download_refs.slurm
#
# Update the USER CONFIGURATION section below before running.
# =============================================================================

# ── USER CONFIGURATION ────────────────────────────────────────────────────────
# Directory where GATK reference files will be stored (must match config.yaml)
RESOURCES_DIR="resources"

# Directory for REDIportal editing database (must match config.yaml rediportal_bed)
EDITING_DIR="resources/editing_db"

# Directory for SnpEff annotation database (must match config.yaml snpeff.data_dir)
SNPEFF_DATA_DIR="resources/snpeff_data"

# Path to SnpEff JAR (set after installing the vcf_annotation conda env)
# conda create -n vcf_annotation -c bioconda snpeff snpsift
SNPEFF_JAR=$(find "$HOME/miniforge3/envs/vcf_annotation" -name "snpEff.jar" 2>/dev/null | head -1)
JAVA_BIN=$(find "$HOME/miniforge3/envs/vcf_annotation" -name "java" 2>/dev/null | head -1)
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "$RESOURCES_DIR" "$EDITING_DIR" "$SNPEFF_DATA_DIR" logs

echo "Starting reference downloads at $(date)"
echo "Resources dir : $RESOURCES_DIR"
echo "Editing DB dir: $EDITING_DIR"
echo "SnpEff data   : $SNPEFF_DATA_DIR"
echo ""

# =============================================================================
# 1. GENCODE v44 GTF
# =============================================================================
echo "[1/5] Downloading Gencode v44 GTF..."
wget -c -P "$RESOURCES_DIR" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
gunzip -f "$RESOURCES_DIR/gencode.v44.primary_assembly.annotation.gtf.gz"
echo "  Done: $RESOURCES_DIR/gencode.v44.primary_assembly.annotation.gtf"


# =============================================================================
# 2. GRCh38 PRIMARY ASSEMBLY FASTA
# =============================================================================
echo "[2/5] Downloading GRCh38 primary assembly FASTA..."
wget -c -P "$RESOURCES_DIR" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
gunzip -f "$RESOURCES_DIR/GRCh38.primary_assembly.genome.fa.gz"
echo "  Done: $RESOURCES_DIR/GRCh38.primary_assembly.genome.fa"


# =============================================================================
# 3. GATK RESOURCE BUNDLE (Broad Institute Google Cloud, hg38/v0)
#    Mills & 1000G gold standard indels  — used by BQSR (Rule 6)
#    dbSNP138 hg38                       — used by BQSR (Rule 6) and rsID
#                                          annotation (Rule 15g)
# =============================================================================
echo "[3/5] Downloading GATK resource bundle (Mills indels + dbSNP138)..."

BROAD="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

# Mills & 1000G indels (must stay bgzipped — tabix index requires it)
wget -c -P "$RESOURCES_DIR" "$BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
wget -c -P "$RESOURCES_DIR" "$BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

# dbSNP138 (Broad bundle provides uncompressed VCF + GATK .idx)
# The pipeline also builds a bgzip+tabix version at runtime (prepare_dbsnp_gz rule)
wget -c -P "$RESOURCES_DIR" "$BROAD/Homo_sapiens_assembly38.dbsnp138.vcf"
wget -c -P "$RESOURCES_DIR" "$BROAD/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

echo "  Done: Mills indels and dbSNP138 in $RESOURCES_DIR/"


# =============================================================================
# 4. REDIPORTAL v2.0 — RNA EDITING DATABASE (hg38)
#    Source : http://srv00.recas.ba.infn.it/atlas/
#    Used by: remove_rna_editing (Rule 15b) to exclude known A-to-I editing
#             sites from the variant call set.
#
#    The raw table uses 1-based coordinates; we convert to 0-based BED format
#    required by bcftools view -T.
# =============================================================================
echo "[4/5] Downloading REDIportal v2.0 hg38 RNA editing database..."

wget -c \
    "http://srv00.recas.ba.infn.it/atlas/TABLE1_hg38.txt.gz" \
    -O "$EDITING_DIR/TABLE1_hg38.txt.gz"

echo "  Converting to 0-based BED format..."
zcat "$EDITING_DIR/TABLE1_hg38.txt.gz" \
    | awk 'NR > 1 { print $1 "\t" ($2 - 1) "\t" $2 }' \
    | sort -k1,1 -k2,2n \
    > "$EDITING_DIR/known_editing_hg38.bed"

rm -f "$EDITING_DIR/TABLE1_hg38.txt.gz"
echo "  Done: $EDITING_DIR/known_editing_hg38.bed"
echo "  Update config.yaml → rna_editing.rediportal_bed to point to this file."


# =============================================================================
# 5. SNPEFF hg38 ANNOTATION DATABASE
#    Source : built-in SnpEff download command (downloads from SourceForge CDN)
#    Used by: annotate_variants (Rule 15d) and health_check_cds (Rule 15f)
#
#    Requires the vcf_annotation conda environment:
#      conda create -n vcf_annotation -c bioconda snpeff snpsift
# =============================================================================
echo "[5/5] Downloading SnpEff hg38 annotation database (~3 GB)..."

if [ -z "$SNPEFF_JAR" ] || [ -z "$JAVA_BIN" ]; then
    echo "  ERROR: SnpEff JAR or Java not found."
    echo "  Install first:  conda create -n vcf_annotation -c bioconda snpeff snpsift"
    echo "  Then re-run this section manually:"
    echo "    java -jar /path/to/snpEff.jar download hg38 -dataDir $SNPEFF_DATA_DIR -v"
else
    "$JAVA_BIN" -Xmx8g \
        -jar "$SNPEFF_JAR" \
        download hg38 \
        -dataDir "$SNPEFF_DATA_DIR" \
        -v
    echo "  Done: SnpEff hg38 database in $SNPEFF_DATA_DIR/"
    echo "  Update config.yaml → snpeff.data_dir to point to this directory."
fi


# =============================================================================
echo ""
echo "========================================================="
echo "All downloads complete at $(date)"
echo ""
echo "Next: update config/config.yaml paths if needed:"
echo "  reference.fasta        → $RESOURCES_DIR/GRCh38.primary_assembly.genome.fa"
echo "  reference.gtf          → $RESOURCES_DIR/gencode.v44.primary_assembly.annotation.gtf"
echo "  reference.dbsnp        → $RESOURCES_DIR/Homo_sapiens_assembly38.dbsnp138.vcf"
echo "  reference.mills_indels → $RESOURCES_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
echo "  rna_editing.rediportal_bed → $EDITING_DIR/known_editing_hg38.bed"
echo "  snpeff.data_dir        → $SNPEFF_DATA_DIR"
echo "========================================================="
