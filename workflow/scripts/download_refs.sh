#!/bin/bash
# =============================================================================
# download_refs.sh
#
# Downloads ALL reference databases required by the pipeline:
#   1. Gencode v44 GTF annotation
#   2. GRCh38 primary assembly FASTA
#   3. GATK resource bundle (Mills indels + dbSNP138)
#   4. REDIportal v2.0 RNA editing database (hg38)
#   5. SnpEff hg38 annotation database
#
# Usage:
#   bash workflow/scripts/download_refs.sh        # run directly
#   sbatch workflow/scripts/download_refs.slurm   # submit via the SLURM wrapper
#
# Output paths match the defaults in config/config.yaml — no edits required.
# Requires the snakemake_env conda env (provides snpEff for step 5):
#   conda env create -f workflow/envs/snakemake.yaml
# Falls back to the vcf_annotation env if snpEff is not in snakemake_env.
# =============================================================================

set -euo pipefail

# ── Paths (must match config/config.yaml) ────────────────────────────────────
RESOURCES_DIR="resources"
EDITING_DIR="resources/editing_db"
SNPEFF_DATA_DIR="resources/snpeff_data"

mkdir -p "$RESOURCES_DIR" "$EDITING_DIR" "$SNPEFF_DATA_DIR" logs

echo "Starting reference downloads at $(date)"
echo "Resources    : $RESOURCES_DIR"
echo "Editing DB   : $EDITING_DIR"
echo "SnpEff data  : $SNPEFF_DATA_DIR"
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
#    Mills & 1000G gold standard indels  — used by BQSR
#    dbSNP138 hg38                       — used by BQSR and rsID annotation
# =============================================================================
echo "[3/5] Downloading GATK resource bundle (Mills indels + dbSNP138)..."

BROAD="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

wget -c -P "$RESOURCES_DIR" "$BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
wget -c -P "$RESOURCES_DIR" "$BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
wget -c -P "$RESOURCES_DIR" "$BROAD/Homo_sapiens_assembly38.dbsnp138.vcf"
wget -c -P "$RESOURCES_DIR" "$BROAD/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

echo "  Done: Mills indels and dbSNP138 in $RESOURCES_DIR/"


# =============================================================================
# 4. REDIPORTAL v2.0 — RNA EDITING DATABASE (hg38)
#    Source : https://rediportal.cloud.ba.infn.it/atlas/download.html
#    File   : TABLE1_hg38_v2.txt.gz  (~90 MB gzipped)
#    Output : 0-based BED for use with bcftools view -T
#    Note   : Older URLs at srv00.recas.ba.infn.it now 404; the host migrated
#             to rediportal.cloud.ba.infn.it/download/ in early 2025.
# =============================================================================
echo "[4/5] Downloading REDIportal v2.0 hg38 RNA editing database..."

# Remove any stale 0-byte / 404-HTML partial from previous URLs
[[ -s "$EDITING_DIR/TABLE1_hg38.txt.gz" ]] || rm -f "$EDITING_DIR/TABLE1_hg38.txt.gz"

wget -c \
    "http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v2.txt.gz" \
    -O "$EDITING_DIR/TABLE1_hg38.txt.gz"

echo "  Converting to 0-based BED format..."
zcat "$EDITING_DIR/TABLE1_hg38.txt.gz" \
    | awk 'NR > 1 { print $1 "\t" ($2 - 1) "\t" $2 }' \
    | sort -k1,1 -k2,2n \
    > "$EDITING_DIR/known_editing_hg38.bed"

rm -f "$EDITING_DIR/TABLE1_hg38.txt.gz"
echo "  Done: $EDITING_DIR/known_editing_hg38.bed"


# =============================================================================
# 5. SNPEFF hg38 ANNOTATION DATABASE
#    Strategy:
#      a) if snpEff is already on PATH, use it
#      b) else try `conda activate` for vcf_annotation, then snakemake_env
#      c) else, if conda + the yaml are available, create vcf_annotation
#         (one-time, ~5 min) and use it
#      d) else, print clear instructions and exit non-zero
#    NOTE: activation must happen in the *parent* shell — using a function
#    invoked via $(...) puts conda activate in a subshell, so PATH does not
#    survive past the function call.
# =============================================================================
echo "[5/5] Downloading SnpEff hg38 annotation database (~3 GB)..."

VCF_YAML="$(dirname "$0")/../envs/vcf_annotation.yaml"
SNPEFF_NEEDS_DEACTIVATE=0

if command -v snpEff >/dev/null 2>&1; then
    : # already available; nothing to do
elif command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    if [[ -n "$CONDA_BASE" ]]; then
        # shellcheck disable=SC1091
        source "$CONDA_BASE/etc/profile.d/conda.sh"

        # (b) try existing envs that are likely to contain snpEff
        for env in vcf_annotation snakemake_env; do
            if conda activate "$env" 2>/dev/null; then
                if command -v snpEff >/dev/null 2>&1; then
                    SNPEFF_NEEDS_DEACTIVATE=1
                    break
                fi
                conda deactivate
            fi
        done

        # (c) still no snpEff — create vcf_annotation from the pinned yaml
        if ! command -v snpEff >/dev/null 2>&1 && [[ -f "$VCF_YAML" ]]; then
            echo "  vcf_annotation env not found — creating it from $VCF_YAML"
            echo "  (one-time, ~5 min)..."
            conda env create -f "$VCF_YAML"
            if conda activate vcf_annotation 2>/dev/null \
                && command -v snpEff >/dev/null 2>&1; then
                SNPEFF_NEEDS_DEACTIVATE=1
            fi
        fi
    fi
fi

if command -v snpEff >/dev/null 2>&1; then
    echo "  Using snpEff from: $(command -v snpEff)"
    snpEff download hg38 -dataDir "$(realpath "$SNPEFF_DATA_DIR")" -v
    echo "  Done: SnpEff hg38 database in $SNPEFF_DATA_DIR/"
    [[ "$SNPEFF_NEEDS_DEACTIVATE" == "1" ]] && conda deactivate
else
    echo "  ERROR: snpEff is not available and could not be installed automatically."
    echo "  Resolve manually with:"
    echo "    conda env create -f workflow/envs/vcf_annotation.yaml"
    echo "    conda activate vcf_annotation"
    echo "    snpEff download hg38 -dataDir $(realpath "$SNPEFF_DATA_DIR") -v"
    exit 1
fi


# =============================================================================
echo ""
echo "========================================================="
echo "All downloads complete at $(date)"
echo "Default config/config.yaml paths already match — no edits needed."
echo "========================================================="
