#!/bin/bash
# Download GSE256519 test data (paired-end, 2x150bp, human heart RNA-seq)
# SRA accessions: SRR28074880 (fetal heart), SRR28074881 (adult heart)
#
# Requires: SRA Toolkit (prefetch + fasterq-dump)
# Install: conda install -c bioconda sra-tools
#
# Usage: bash workflow/scripts/download_GSE256519.sh

set -euo pipefail

OUTDIR="data/GSE256519"
TMPDIR="data/GSE256519/sra_cache"
ACCESSIONS=("SRR28074879")

mkdir -p "$OUTDIR" "$TMPDIR"

for ACC in "${ACCESSIONS[@]}"; do
    echo "========================================="
    echo "Downloading $ACC ..."
    echo "========================================="

    # Step 1: prefetch — downloads compressed .sra file (resumable)
    prefetch \
        --output-directory "$TMPDIR" \
        --max-size 200G \
        "$ACC"

    # Step 2: fasterq-dump — converts .sra to paired FASTQ
    # --split-files produces {ACC}_1.fastq and {ACC}_2.fastq
    fasterq-dump \
        --outdir "$OUTDIR" \
        --temp "$TMPDIR" \
        --split-files \
        --threads 8 \
        --progress \
        "$TMPDIR/$ACC/$ACC.sra"

    echo "Done: $OUTDIR/${ACC}_1.fastq  $OUTDIR/${ACC}_2.fastq"
    echo ""
done

# Clean up .sra cache to save space
echo "Cleaning up SRA cache..."
rm -rf "$TMPDIR"

echo ""
echo "All downloads complete. Files in $OUTDIR:"
ls -lh "$OUTDIR"
