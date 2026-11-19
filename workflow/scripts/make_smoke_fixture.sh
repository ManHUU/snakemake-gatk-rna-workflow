#!/bin/bash
# =============================================================================
# make_smoke_fixture.sh  (reviewer comments C-07 / C-08 / C-09)
#
# Builds a MINIATURE, single-contig fixture for the smoke test from the full
# references already downloaded under resources/ and the demo FASTQs under
# data/GSE256519/. The contig is subset WHOLE so assembly coordinates are
# preserved and all derived resources (sequence dictionary, indexes, GTF,
# known-sites, REDIportal) stay mutually compatible.
#
# Outputs (consumed by config/smoke/config.yaml):
#   resources/smoke/smoke.fa (+ .fai, .dict)     miniature reference
#   resources/smoke/smoke.gtf                     annotation for the contig
#   resources/smoke/smoke.dbsnp138.vcf            dbSNP subset
#   resources/smoke/smoke.mills.indels.vcf.gz     Mills indels subset
#   resources/smoke/smoke_known_editing.bed       REDIportal subset
#   resources/smoke/star_index/                   small-genome STAR index
#   data/smoke/<sample>_1.fastq / _2.fastq        reads covering the contig
#   resources/smoke/fixture_manifest.sha256       checksums of the above
#
# Usage:
#   bash workflow/scripts/make_smoke_fixture.sh [CONTIG]     # default chr22
#
# Requires: apptainer/singularity on PATH and the full references present
# (run workflow/scripts/download_refs.sh and download_GSE256519.sh first).
# This builds inputs only; it runs no pipeline step.
# =============================================================================

set -euo pipefail

CONTIG="${1:-chr22}"
SA_INDEX_NBASES="${SMOKE_SA_INDEX_NBASES:-11}"   # keep in sync with config/smoke/config.yaml
SJDB_OVERHANG="${SMOKE_SJDB_OVERHANG:-149}"
# Keep the fixture LIGHT: only the first MAX_INPUT_PAIRS read pairs of each demo
# sample are aligned to the one-contig mini index, and only the pairs that map
# to the contig are kept. The demo FASTQs are ~200M pairs each, so this avoids
# aligning the full dataset while still yielding reads that cover the contig.
MAX_INPUT_PAIRS="${SMOKE_MAX_INPUT_PAIRS:-5000000}"
THREADS="${SMOKE_THREADS:-4}"

# Full references (must match config/config.yaml defaults)
REF_FA="resources/GRCh38.primary_assembly.genome.fa"
REF_GTF="resources/gencode.v44.primary_assembly.annotation.gtf"
DBSNP="resources/Homo_sapiens_assembly38.dbsnp138.vcf"
MILLS="resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
EDIT_BED="resources/editing_db/known_editing_hg38.bed"
DEMO_DIR="data/GSE256519"

OUT="resources/smoke"
FQ_OUT="data/smoke"
mkdir -p "$OUT" "$FQ_OUT"

# Container runtime
if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else echo "ERROR: need apptainer or singularity on PATH." >&2; exit 1; fi

STAR_C="docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_3"
SAMTOOLS_C="docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"
GATK_C="docker://broadinstitute/gatk:4.6.1.0"
BCFTOOLS_C="docker://quay.io/biocontainers/bcftools:1.21--h3a4d415_1"

BIND="$(pwd)"
[[ -n "${HPC_SCRATCH_DIR:-}" && -d "${HPC_SCRATCH_DIR}" ]] && BIND="${BIND},${HPC_SCRATCH_DIR}"
run() { local img="$1"; shift; "$RT" exec --bind "$BIND" "$img" "$@"; }

for f in "$REF_FA" "$REF_GTF" "$DBSNP" "$MILLS" "$EDIT_BED"; do
    [[ -e "$f" ]] || { echo "ERROR: missing reference '$f' (run download_refs.sh first)." >&2; exit 1; }
done

echo "[1/7] Subsetting reference FASTA to $CONTIG ..."
run "$SAMTOOLS_C" samtools faidx "$REF_FA" "$CONTIG" > "$OUT/smoke.fa"
run "$SAMTOOLS_C" samtools faidx "$OUT/smoke.fa"
rm -f "$OUT/smoke.dict"
# samtools dict (not gatk CreateSequenceDictionary): equivalent .dict output and
# avoids the gatk python-wrapper dependency during fixture prep. The pipeline's
# prepare_dict rule still uses gatk at run time (covered by the Snakefile shim).
run "$SAMTOOLS_C" samtools dict "$OUT/smoke.fa" -o "$OUT/smoke.dict"

echo "[2/7] Subsetting GTF to $CONTIG ..."
awk -v c="$CONTIG" '/^#/ {print; next} $1==c {print}' "$REF_GTF" > "$OUT/smoke.gtf"

echo "[3/7] Subsetting dbSNP138 to $CONTIG ..."
# Index a bgzipped copy if needed, then region-select.
if [[ ! -e "${DBSNP}.gz.tbi" ]]; then
    run "$BCFTOOLS_C" bash -c "bgzip -c '$DBSNP' > '${DBSNP}.gz' && tabix -p vcf '${DBSNP}.gz'"
fi
run "$BCFTOOLS_C" bcftools view -r "$CONTIG" "${DBSNP}.gz" -O v -o "$OUT/smoke.dbsnp138.vcf"

echo "[4/7] Subsetting Mills indels to $CONTIG ..."
[[ -e "${MILLS}.tbi" ]] || run "$BCFTOOLS_C" bcftools index -t "$MILLS"
run "$BCFTOOLS_C" bcftools view -r "$CONTIG" "$MILLS" -O z -o "$OUT/smoke.mills.indels.vcf.gz"
run "$BCFTOOLS_C" bcftools index -t "$OUT/smoke.mills.indels.vcf.gz"

echo "[5/7] Subsetting REDIportal editing BED to $CONTIG ..."
awk -v c="$CONTIG" '$1==c' "$EDIT_BED" > "$OUT/smoke_known_editing.bed"

echo "[6/7] Building small-genome STAR index (SAindexNbases=$SA_INDEX_NBASES) ..."
rm -rf "$OUT/star_index"; mkdir -p "$OUT/star_index"
run "$STAR_C" STAR --runThreadN 4 --runMode genomeGenerate \
    --genomeDir "$OUT/star_index" \
    --genomeFastaFiles "$OUT/smoke.fa" \
    --sjdbGTFfile "$OUT/smoke.gtf" \
    --sjdbOverhang "$SJDB_OVERHANG" \
    --genomeSAindexNbases "$SA_INDEX_NBASES"

echo "[7/7] Extracting a LIGHT set of reads covering $CONTIG from the demo FASTQs ..."
echo "      (slicing the first $MAX_INPUT_PAIRS pairs/sample, then keeping only $CONTIG-mapped pairs)"
tmpd="$(mktemp -d)"
trap 'rm -rf "$tmpd"' EXIT
shopt -s nullglob
slice_lines=$(( MAX_INPUT_PAIRS * 4 ))
for r1 in "$DEMO_DIR"/*_1.fastq; do
    sample="$(basename "$r1" _1.fastq)"
    r2="$DEMO_DIR/${sample}_2.fastq"
    [[ -e "$r2" ]] || { echo "  skip $sample (no R2)"; continue; }

    echo "  [$sample] slicing first $MAX_INPUT_PAIRS pairs ..."
    head -n "$slice_lines" "$r1" > "$tmpd/${sample}_1.slice.fastq"
    head -n "$slice_lines" "$r2" > "$tmpd/${sample}_2.slice.fastq"

    echo "  [$sample] aligning slice to the $CONTIG mini index ..."
    run "$STAR_C" STAR --runThreadN "$THREADS" --genomeDir "$OUT/star_index" \
        --readFilesIn "$tmpd/${sample}_1.slice.fastq" "$tmpd/${sample}_2.slice.fastq" \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix "$tmpd/${sample}."

    # Keep read pairs where both mates mapped to the contig (-f 0x2 = proper pair),
    # then stream back to paired FASTQ. This is the LIGHT fixture input.
    run "$SAMTOOLS_C" bash -c \
        "samtools collate -@${THREADS} -O '$tmpd/${sample}.Aligned.out.bam' '$tmpd/${sample}.collate' \
         | samtools fastq -f 0x2 -1 '$FQ_OUT/${sample}_1.fastq' -2 '$FQ_OUT/${sample}_2.fastq' -0 /dev/null -s /dev/null -n"

    n1=$(( $(wc -l < "$FQ_OUT/${sample}_1.fastq") / 4 ))
    sz=$(du -h "$FQ_OUT/${sample}_1.fastq" | cut -f1)
    echo "  [$sample] kept $n1 read pairs covering $CONTIG (R1 $sz)"
    rm -f "$tmpd/${sample}"*.slice.fastq "$tmpd/${sample}.Aligned.out.bam"
done

echo "Writing checksum manifest ..."
( cd "$OUT" && find . -type f ! -name fixture_manifest.sha256 -exec sha256sum {} \; | sort -k2 ) \
    > "$OUT/fixture_manifest.sha256"
( cd "$FQ_OUT" && find . -type f -exec sha256sum {} \; | sort -k2 ) \
    >> "$OUT/fixture_manifest.sha256"

echo
echo "Smoke fixture ready."
echo "  reference : $OUT/"
echo "  reads     : $FQ_OUT/"
echo "Next: bash workflow/scripts/run_smoke_test.sh"
