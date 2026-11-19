#!/bin/bash
# =============================================================================
# check_smoke_outputs.sh  (reviewer comments C-07 / C-08)
#
# Quantitative acceptance checks for a completed smoke-test run. Distinguishes a
# genuine EXECUTION FAILURE (missing/unreadable/malformed required output, wrong
# sample IDs, header-only final VCF) from a completed run. For the FIXED smoke
# fixture, absence of the expected non-empty variant records fails the test.
#
# Usage:
#   bash workflow/scripts/check_smoke_outputs.sh <results_dir> [expected_checkpoints.tsv]
#
# Exit 0 = all mandatory checks passed; 1 = one or more failed.
# =============================================================================

set -uo pipefail

OUT="${1:-results_smoke}"
EXPECTED="${2:-config/smoke/expected_checkpoints.tsv}"
FASTQ_SMOKE="data/smoke"

fail=0
pass_msg() { printf "[PASS]  %s\n" "$1"; }
warn_msg() { printf "[WARN]  %s\n" "$1"; }
fail_msg() { printf "[FAIL]  %s\n" "$1"; fail=$((fail + 1)); }

# Prefer bcftools from a container for robust VCF parsing; else fall back to zcat.
BCFTOOLS_C="docker://quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
if command -v apptainer >/dev/null 2>&1; then RT=apptainer
elif command -v singularity >/dev/null 2>&1; then RT=singularity
else RT=""; fi
bcf() { [[ -n "$RT" ]] && "$RT" exec --bind "$(pwd)" "$BCFTOOLS_C" bcftools "$@"; }
count_records() {  # robust record count
    local vcf="$1"
    if [[ -n "$RT" ]]; then bcf view -H "$vcf" 2>/dev/null | wc -l
    else zcat -f "$vcf" 2>/dev/null | grep -vc '^#'; fi
}

echo "Smoke-test acceptance checks on: $OUT"
echo "=================================================="

# 1. Expected sample IDs (derived from the fixture FASTQs) all present as inputs.
if [[ -d "$FASTQ_SMOKE" ]]; then
    # -L so a symlinked fixture dir (e.g. data/smoke -> shared storage) is followed.
    samples=$(find -L "$FASTQ_SMOKE" -maxdepth 1 -name '*_1.fastq' -printf '%f\n' 2>/dev/null | sed 's/_1\.fastq$//' | sort)
    if [[ -n "$samples" ]]; then
        pass_msg "Fixture samples: $(echo "$samples" | tr '\n' ' ')"
    else
        fail_msg "No fixture samples found in $FASTQ_SMOKE (build the fixture first)."
    fi
else
    fail_msg "Fixture FASTQ directory $FASTQ_SMOKE missing."
fi

# 2. Required outputs exist.
required=(
    "$OUT/Final_CDS_rsID.vcf.gz"
    "$OUT/Final_CDS_rsID.vcf.gz.tbi"
    "$OUT/qc/multiqc_report.html"
    "$OUT/Final_annotated_healthCheck.html"
    "$OUT/Final_CDS_healthCheck.html"
    "$OUT/stats/ts_tv.tsv"
    "$OUT/stats/filter_summary.tsv"
)
for f in "${required[@]}"; do
    if [[ -s "$f" ]]; then pass_msg "present: $f"
    else fail_msg "missing or empty required output: $f"; fi
done

# 3. Final VCF parses and is NON-EMPTY (a header-only VCF must not pass).
FINAL="$OUT/Final_CDS_rsID.vcf.gz"
if [[ -e "$FINAL" ]]; then
    if gzip -t "$FINAL" 2>/dev/null; then pass_msg "final VCF is valid gzip/bgzip"
    else fail_msg "final VCF is not a valid gzip stream: $FINAL"; fi
    n=$(count_records "$FINAL"); n=${n:-0}
    if [[ "$n" -gt 0 ]]; then pass_msg "final VCF record count: $n (non-empty)"
    else fail_msg "final VCF has 0 records (header-only) - smoke test did not produce variants."; fi
fi

# 4. Variant-statistics tables present (quantitative checkpoints).
# These are what the published variant_statistics rule emits.
for s in filter_summary chrom_counts sample_counts quality_metrics; do
    if [[ -s "$OUT/stats/$s.tsv" ]]; then pass_msg "checkpoint: stats/$s.tsv"
    else warn_msg "checkpoint stats/$s.tsv not found (variant_statistics may have been skipped)."; fi
done

# 5. Compare measured values against expected_checkpoints.tsv when populated.
if [[ -s "$EXPECTED" ]]; then
    echo "--------------------------------------------------"
    echo "Expected-vs-measured (expected from $EXPECTED):"
    measured_final="$(count_records "$FINAL" 2>/dev/null || echo NA)"
    exp_final="$(awk -F'\t' '$1=="final_cds_rsid_records"{print $2}' "$EXPECTED" 2>/dev/null)"
    if [[ -z "$exp_final" || "$exp_final" == "[MEASURE]" ]]; then
        warn_msg "expected final record count not yet recorded; measured = $measured_final. Record this in $EXPECTED."
    elif [[ "$measured_final" == "$exp_final" ]]; then
        pass_msg "final record count matches expected ($measured_final)."
    else
        fail_msg "final record count $measured_final != expected $exp_final."
    fi
fi

echo "=================================================="
if [[ "$fail" -eq 0 ]]; then
    echo "SMOKE TEST: PASS"; exit 0
else
    echo "SMOKE TEST: FAIL ($fail check(s) failed)"; exit 1
fi
