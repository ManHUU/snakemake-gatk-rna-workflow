#!/bin/bash
# =============================================================================
# make_units.sh — generate config/units.tsv from a FASTQ directory.
#
# Run this ONCE to (re)build the sample sheet; it is NOT part of the pipeline
# runtime. The pipeline itself reads config/units.tsv and never globs, so the
# exact set of analysed samples is explicit and version-controlled.
#
# Usage:
#   bash workflow/scripts/make_units.sh <fastq_dir> [paired|single] [out.tsv]
#
# Examples:
#   bash workflow/scripts/make_units.sh data/GSE256519 paired
#   bash workflow/scripts/make_units.sh data/my_se single config/units.tsv
#
# Naming convention assumed (same as before):
#   paired  → {sample}_1.fastq + {sample}_2.fastq
#   single  → {sample}.fastq
#
# Review the generated file before running the pipeline — this is exactly the
# point where a half-downloaded or mis-named FASTQ becomes visible.
# =============================================================================

set -euo pipefail

FASTQ_DIR="${1:?Usage: make_units.sh <fastq_dir> [paired|single] [out.tsv]}"
SEQ_TYPE="${2:-paired}"
OUT="${3:-config/units.tsv}"

mkdir -p "$(dirname "$OUT")"
printf 'sample\tunit\tfq1\tfq2\n' > "$OUT"

count=0
if [[ "$SEQ_TYPE" == "paired" ]]; then
    for r1 in "$FASTQ_DIR"/*_1.fastq; do
        [[ -e "$r1" ]] || { echo "ERROR: no *_1.fastq in $FASTQ_DIR"; exit 1; }
        sample=$(basename "$r1" _1.fastq)
        r2="$FASTQ_DIR/${sample}_2.fastq"
        [[ -e "$r2" ]] || { echo "ERROR: missing mate $r2"; exit 1; }
        printf '%s\t1\t%s\t%s\n' "$sample" "$r1" "$r2" >> "$OUT"
        count=$((count + 1))
    done
elif [[ "$SEQ_TYPE" == "single" ]]; then
    for fq in "$FASTQ_DIR"/*.fastq; do
        [[ -e "$fq" ]] || { echo "ERROR: no *.fastq in $FASTQ_DIR"; exit 1; }
        sample=$(basename "$fq" .fastq)
        printf '%s\t1\t%s\t\n' "$sample" "$fq" >> "$OUT"
        count=$((count + 1))
    done
else
    echo "ERROR: sequencing type must be 'paired' or 'single' (got '$SEQ_TYPE')."
    exit 1
fi

echo "Wrote $count sample(s) to $OUT"
echo "Review it, then make sure config.yaml sequencing_type matches '$SEQ_TYPE'."
