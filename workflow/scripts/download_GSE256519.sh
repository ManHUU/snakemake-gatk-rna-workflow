#!/bin/bash
# Download GSE256519 test data (paired-end, 2x150bp, human heart RNA-seq)
# SRA accessions: SRR28074879, SRR28074880
#
# Requires: SRA Toolkit (prefetch + fasterq-dump)
# Install:
#   conda env create -f workflow/envs/sra-tools.yaml
#   conda activate sra_tools
#
# Usage:
#   bash workflow/scripts/download_GSE256519.sh
#   THREADS=4 bash workflow/scripts/download_GSE256519.sh

set -euo pipefail

OUTDIR="data/GSE256519"
TMPDIR="data/GSE256519/sra_cache"
THREADS="${THREADS:-8}"
ACCESSIONS=("SRR28074879" "SRR28074880")

for cmd in prefetch fasterq-dump; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found in PATH."
        echo ""
        echo "Install the required SRA Toolkit environment with:"
        echo "  conda env create -f workflow/envs/sra-tools.yaml"
        echo "  conda activate sra_tools"
        exit 1
    fi
done

mkdir -p "$OUTDIR" "$TMPDIR"

for ACC in "${ACCESSIONS[@]}"; do
    R1="$OUTDIR/${ACC}_1.fastq"
    R2="$OUTDIR/${ACC}_2.fastq"
    SRA="$TMPDIR/$ACC/$ACC.sra"

    echo "========================================="
    echo "Downloading $ACC ..."
    echo "========================================="

    if [[ -s "$R1" && -s "$R2" ]]; then
        echo "Skipping $ACC: FASTQ files already exist."
        echo "  $R1"
        echo "  $R2"
        echo ""
        continue
    fi

    # Step 1: prefetch — downloads compressed .sra file (resumable)
    prefetch \
        --output-directory "$TMPDIR" \
        --max-size 200G \
        "$ACC"

    if [[ ! -s "$SRA" ]]; then
        echo "ERROR: Expected SRA file was not created: $SRA"
        exit 1
    fi

    # Step 2: fasterq-dump — converts .sra to paired FASTQ
    # --split-files produces {ACC}_1.fastq and {ACC}_2.fastq
    fasterq-dump \
        --outdir "$OUTDIR" \
        --temp "$TMPDIR" \
        --split-files \
        --threads "$THREADS" \
        --progress \
        "$SRA"

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "ERROR: FASTQ conversion failed for $ACC."
        echo "Expected:"
        echo "  $R1"
        echo "  $R2"
        exit 1
    fi

    echo "Done: $R1  $R2"
    echo ""
done

# Clean up .sra cache to save space
echo "Cleaning up SRA cache..."
rm -rf "$TMPDIR"

echo ""
echo "All downloads complete. Files in $OUTDIR:"
ls -lh "$OUTDIR"
