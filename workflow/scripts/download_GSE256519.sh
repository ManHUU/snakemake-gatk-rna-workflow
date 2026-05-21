#!/bin/bash
# Download GSE256519 test data (paired-end, 2x150bp, human heart RNA-seq)
# SRA accessions: SRR28074879, SRR28074880
#
# Requires: apptainer or singularity on PATH (the SRA Toolkit runs from the
# pinned biocontainer image — no conda env needed on the host).
#
# Usage:
#   bash workflow/scripts/download_GSE256519.sh
#   THREADS=4 bash workflow/scripts/download_GSE256519.sh

set -euo pipefail

OUTDIR="data/GSE256519"
TMPDIR="data/GSE256519/sra_cache"
THREADS="${THREADS:-8}"
ACCESSIONS=("SRR28074879" "SRR28074880")

SRA_IMAGE="docker://quay.io/biocontainers/sra-tools:3.1.1--h4304569_2"

# Ensure prefetch + fasterq-dump are available either directly or via container.
have_sra_tools() {
    command -v prefetch >/dev/null 2>&1 && command -v fasterq-dump >/dev/null 2>&1
}

if have_sra_tools; then
    # Already on PATH (e.g. user has sra-tools installed system-wide).
    RUN=()
elif command -v apptainer >/dev/null 2>&1; then
    RUN=(apptainer exec --bind "$(pwd)" "$SRA_IMAGE")
elif command -v singularity >/dev/null 2>&1; then
    RUN=(singularity exec --bind "$(pwd)" "$SRA_IMAGE")
else
    echo "ERROR: SRA Toolkit not on PATH, and neither apptainer nor singularity"
    echo "       was found to run the biocontainer image."
    echo "Install apptainer (preferred) or singularity, or 'module load apptainer'"
    echo "on HPC, and re-run this script."
    exit 1
fi

if [[ "${#RUN[@]}" -gt 0 ]]; then
    echo "Using sra-tools from container: $SRA_IMAGE"
else
    echo "Using sra-tools from PATH: $(command -v prefetch)"
fi

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
    "${RUN[@]+"${RUN[@]}"}" prefetch \
        --output-directory "$TMPDIR" \
        --max-size 200G \
        "$ACC"

    if [[ ! -s "$SRA" ]]; then
        echo "ERROR: Expected SRA file was not created: $SRA"
        exit 1
    fi

    # Step 2: fasterq-dump — converts .sra to paired FASTQ
    # --split-files produces {ACC}_1.fastq and {ACC}_2.fastq
    "${RUN[@]+"${RUN[@]}"}" fasterq-dump \
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
