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
# Requires apptainer or singularity on PATH (step 5/5 runs snpEff from the
# biocontainer image — no conda env needed).
# =============================================================================

set -euo pipefail

# ── Paths (must match config/config.yaml) ────────────────────────────────────
RESOURCES_DIR="resources"
EDITING_DIR="resources/editing_db"
SNPEFF_DATA_DIR="resources/snpeff_data"

mkdir -p "$RESOURCES_DIR" "$EDITING_DIR" "$SNPEFF_DATA_DIR" logs

# ── Expected sha256 checksums (generated once from verified reference copies) ─
# GTF / FASTA hashes are of the DECOMPRESSED files (upstream serves .gz; we
# verify after gunzip). The REDIportal BED is the converted/sorted output.
SHA_GTF="544c057e3c124d26c9c7f89d134318f4ca14ba3c4d6af71707b3dce93ad36812"
SHA_FASTA="e49b92b3e4f321bf254c042f25b726d9931c4d74c7523e8b6bb530e63b0cfd4b"
SHA_MILLS="e11195ef119927518f01965bfe63591ec4c761b727ed26e72f17fed23a6739e5"
SHA_MILLS_TBI="6856f5ad5690e1f095960ecee450e26a67f9486a8336aa08be4f599d958c1bac"
SHA_DBSNP="0ff368a3a4e16fc539be19ddbc58996d7e9a959a50e8597a7de07e217cce8527"
SHA_DBSNP_IDX="d0c185e0b5b43f3c2587e01c93f0c0dc81754649c0e0de69af64defa74cb33f8"
SHA_REDI_BED="2716906413d302982f191cbf2f8601d8661be00b756eac6b4278a9f03eadfbf7"
SHA_SNPEFF_TGZ="420dff0211786de1bd875c4f92b1490f99327c32952cde98c8fc9a8d4d953a59"

# Zenodo fallback mirror (record 21030946) for the two link-rot-prone files.
ZENODO="https://zenodo.org/records/21030946/files"

# ── Helper functions ─────────────────────────────────────────────────────────
# verify_sha256 <file> <expected>
#   Return 0 if <file>'s sha256 matches <expected>. Pass "SKIP" to disable.
verify_sha256() {
    local file="$1" expected="$2" actual
    [[ "$expected" == "SKIP" ]] && return 0
    [[ -f "$file" ]] || { echo "  [checksum] missing file: $file"; return 1; }
    actual=$(sha256sum "$file" | awk '{print $1}')
    if [[ "$actual" == "$expected" ]]; then
        echo "  [checksum] OK: $file"
        return 0
    fi
    echo "  [checksum] MISMATCH: $file"
    echo "             expected: $expected"
    echo "             actual:   $actual"
    return 1
}

# download_one <url> <out>
#   wget with a few retries; resumes a partial within a URL. No checksum.
download_one() {
    local url="$1" out="$2" attempt
    mkdir -p "$(dirname "$out")"
    for attempt in 1 2 3; do
        echo "  [fetch] $url (attempt $attempt/3)"
        if wget --tries=2 --timeout=60 -c -O "$out" "$url"; then
            return 0
        fi
        sleep 5
    done
    return 1
}

# fetch_with_fallback <primary_url> <backup_url> <out_path> <expected_sha256>
#   Download primary → out and verify. On any failure, try the backup URL.
#   Both sources must serve a file matching <expected_sha256>. Pass backup=""
#   to skip the fallback. Skips work if <out_path> already verifies.
fetch_with_fallback() {
    local primary="$1" backup="$2" out="$3" sha="$4" url
    if [[ -f "$out" ]] && verify_sha256 "$out" "$sha" >/dev/null 2>&1; then
        echo "  [fetch] already present and verified: $out"
        return 0
    fi
    for url in "$primary" "$backup"; do
        [[ -z "$url" ]] && continue
        rm -f "$out"
        if download_one "$url" "$out" && verify_sha256 "$out" "$sha"; then
            return 0
        fi
        echo "  [fetch] source failed (download or checksum): $url"
    done
    echo "  ERROR: all sources failed for $out"
    return 1
}

echo "Starting reference downloads at $(date)"
echo "Resources    : $RESOURCES_DIR"
echo "Editing DB   : $EDITING_DIR"
echo "SnpEff data  : $SNPEFF_DATA_DIR"
echo ""

# =============================================================================
# 1. GENCODE v44 GTF
# =============================================================================
echo "[1/5] Downloading Gencode v44 GTF..."
GTF_OUT="$RESOURCES_DIR/gencode.v44.primary_assembly.annotation.gtf"
if verify_sha256 "$GTF_OUT" "$SHA_GTF" >/dev/null 2>&1; then
    echo "  [fetch] already present and verified: $GTF_OUT"
else
    download_one \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" \
        "$GTF_OUT.gz"
    gunzip -f "$GTF_OUT.gz"
    verify_sha256 "$GTF_OUT" "$SHA_GTF"
fi
echo "  Done: $GTF_OUT"


# =============================================================================
# 2. GRCh38 PRIMARY ASSEMBLY FASTA
# =============================================================================
echo "[2/5] Downloading GRCh38 primary assembly FASTA..."
FASTA_OUT="$RESOURCES_DIR/GRCh38.primary_assembly.genome.fa"
if verify_sha256 "$FASTA_OUT" "$SHA_FASTA" >/dev/null 2>&1; then
    echo "  [fetch] already present and verified: $FASTA_OUT"
else
    download_one \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz" \
        "$FASTA_OUT.gz"
    gunzip -f "$FASTA_OUT.gz"
    verify_sha256 "$FASTA_OUT" "$SHA_FASTA"
fi
echo "  Done: $FASTA_OUT"


# =============================================================================
# 3. GATK RESOURCE BUNDLE (Broad Institute Google Cloud, hg38/v0)
#    Mills & 1000G gold standard indels  — used by BQSR
#    dbSNP138 hg38                       — used by BQSR and rsID annotation
# =============================================================================
echo "[3/5] Downloading GATK resource bundle (Mills indels + dbSNP138)..."

BROAD="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

# No Zenodo mirror needed — the Broad public-data bucket is a stable host.
fetch_with_fallback "$BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" "" \
    "$RESOURCES_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" "$SHA_MILLS"
fetch_with_fallback "$BROAD/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" "" \
    "$RESOURCES_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" "$SHA_MILLS_TBI"
fetch_with_fallback "$BROAD/Homo_sapiens_assembly38.dbsnp138.vcf" "" \
    "$RESOURCES_DIR/Homo_sapiens_assembly38.dbsnp138.vcf" "$SHA_DBSNP"
fetch_with_fallback "$BROAD/Homo_sapiens_assembly38.dbsnp138.vcf.idx" "" \
    "$RESOURCES_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx" "$SHA_DBSNP_IDX"

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

REDI_BED="$EDITING_DIR/known_editing_hg38.bed"
if verify_sha256 "$REDI_BED" "$SHA_REDI_BED" >/dev/null 2>&1; then
    echo "  [fetch] already present and verified: $REDI_BED"
else
    # Primary: upstream table → convert to 0-based BED (LC_ALL=C keeps the sort
    # byte-reproducible, so the result matches the pinned checksum / Zenodo copy).
    # Fallback: pull the ready-made BED from the Zenodo archive.
    rm -f "$EDITING_DIR/TABLE1_hg38.txt.gz"
    if download_one "http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v2.txt.gz" \
                    "$EDITING_DIR/TABLE1_hg38.txt.gz"; then
        echo "  Converting to 0-based BED format..."
        zcat "$EDITING_DIR/TABLE1_hg38.txt.gz" \
            | awk 'NR > 1 { print $1 "\t" ($2 - 1) "\t" $2 }' \
            | LC_ALL=C sort -k1,1 -k2,2n \
            > "$REDI_BED"
        rm -f "$EDITING_DIR/TABLE1_hg38.txt.gz"
    fi
    if ! verify_sha256 "$REDI_BED" "$SHA_REDI_BED"; then
        echo "  Upstream unavailable or checksum failed — falling back to Zenodo BED..."
        fetch_with_fallback "$ZENODO/known_editing_hg38.bed?download=1" "" \
            "$REDI_BED" "$SHA_REDI_BED"
    fi
fi
echo "  Done: $REDI_BED"


# =============================================================================
# 5. SNPEFF hg38 ANNOTATION DATABASE
#    Strategy: run snpEff from the pinned biocontainer via apptainer (or
#    singularity).  Same image used by the pipeline itself — no conda
#    env needed on the host.
# =============================================================================
echo "[5/5] Downloading SnpEff hg38 annotation database (~3 GB)..."

SNPEFF_IMAGE="docker://quay.io/biocontainers/snpeff:5.3.0a--hdfd78af_1"

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_RUN=(apptainer exec --bind "$(pwd)" "$SNPEFF_IMAGE")
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_RUN=(singularity exec --bind "$(pwd)" "$SNPEFF_IMAGE")
else
    echo "  ERROR: neither apptainer nor singularity found on PATH."
    echo "  One of them is required to run snpEff for the database download."
    exit 1
fi

# Skip if a usable database is already extracted.
if [[ -f "$SNPEFF_DATA_DIR/hg38/snpEffectPredictor.bin" ]]; then
    echo "  [fetch] SnpEff hg38 database already present: $SNPEFF_DATA_DIR/hg38/"
else
    # Primary: let snpEff pull the versioned db from its official host.
    echo "  Pulling/using snpEff container: $SNPEFF_IMAGE"
    if "${CONTAINER_RUN[@]}" snpEff download hg38 \
            -dataDir "$(realpath "$SNPEFF_DATA_DIR")" -v; then
        echo "  SnpEff download succeeded."
    else
        # Fallback: checksum-verified archive on Zenodo (extracts to hg38/).
        echo "  snpEff download failed — falling back to Zenodo archive..."
        SNPEFF_TGZ="$SNPEFF_DATA_DIR/snpeff_hg38.tar.gz"
        fetch_with_fallback "$ZENODO/snpeff_hg38.tar.gz?download=1" "" \
            "$SNPEFF_TGZ" "$SHA_SNPEFF_TGZ"
        echo "  Extracting $SNPEFF_TGZ ..."
        tar xzf "$SNPEFF_TGZ" -C "$SNPEFF_DATA_DIR"
        rm -f "$SNPEFF_TGZ"
    fi
fi
echo "  Done: SnpEff hg38 database in $SNPEFF_DATA_DIR/"


# =============================================================================
echo ""
echo "========================================================="
echo "All downloads complete at $(date)"
echo "Default config/config.yaml paths already match — no edits needed."
echo "========================================================="
