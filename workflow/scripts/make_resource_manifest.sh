#!/bin/bash
# =============================================================================
# make_resource_manifest.sh  (reviewer comment G-04)
#
# Records the exact reference resources used by the workflow into a versioned
# manifest: role, assembly/release identifier, source URL, chromosome naming
# convention, and a sha256 checksum of the file actually present. The manuscript
# resource-compatibility table links to this manifest so the example run is
# reproducible and users can tell when a swapped resource is incompatible.
#
# Usage:
#   bash workflow/scripts/make_resource_manifest.sh
# Output:
#   resources/resource_manifest.tsv
# =============================================================================

set -uo pipefail

CFG="config/config.yaml"
OUT="resources/resource_manifest.tsv"

sum() {  # sha256 of a file if present, else "absent"
    [ -e "$1" ] && sha256sum "$1" 2>/dev/null | awk '{print $1}' || echo "absent"
}
naming() {  # report chromosome naming of the first record of a FASTA/GTF/BED/VCF
    local f="$1"
    [ -e "$f" ] || { echo "n/a"; return; }
    case "$f" in
        *.fa|*.fasta) awk '/^>/{print ($0 ~ /^>chr/)?"chr-prefixed":"no-chr-prefix"; exit}' "$f" ;;
        *.gtf)        awk '!/^#/{print ($1 ~ /^chr/)?"chr-prefixed":"no-chr-prefix"; exit}' "$f" ;;
        *.bed)        awk '{print ($1 ~ /^chr/)?"chr-prefixed":"no-chr-prefix"; exit}' "$f" ;;
        *.vcf)        awk '!/^#/{print ($1 ~ /^chr/)?"chr-prefixed":"no-chr-prefix"; exit}' "$f" ;;
        *)            echo "unknown" ;;
    esac
}

printf 'resource\trole\tassembly_or_release\tsource_url\tchrom_naming\tfile\tsha256\n' > "$OUT"

row() {  # row <resource> <role> <assembly/release> <url> <file>
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$1" "$2" "$3" "$4" "$(naming "$5")" "$5" "$(sum "$5")" >> "$OUT"
}

row "Reference FASTA" "STAR index, GATK reference" "GRCh38 primary assembly, GENCODE release 44" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz" \
    "resources/GRCh38.primary_assembly.genome.fa"
row "Gene annotation GTF" "STAR splice junctions" "GENCODE release 44 (GRCh38)" \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" \
    "resources/gencode.v44.primary_assembly.annotation.gtf"
row "dbSNP138" "BQSR known-sites, rsID annotation" "Homo_sapiens_assembly38 (hg38, Broad bundle v0)" \
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf" \
    "resources/Homo_sapiens_assembly38.dbsnp138.vcf"
row "Mills & 1000G indels" "BQSR known-sites" "hg38 (Broad bundle v0)" \
    "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
    "resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
row "REDIportal editing BED" "RNA-editing exclusion" "REDIportal v2.0 (hg38), TABLE1_hg38_v2" \
    "http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v2.txt.gz" \
    "resources/editing_db/known_editing_hg38.bed"

# SnpEff database: record the installed database directory name/version if present.
SNPEFF_DIR="$(awk '/data_dir:/{gsub(/"/,"",$2); print $2; exit}' "$CFG" 2>/dev/null || echo resources/snpeff_data)"
snpeff_db="$(ls -d "${SNPEFF_DIR}"/hg38* 2>/dev/null | head -1)"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "SnpEff database" "Functional annotation" "hg38 (SnpEff; record bundle version from snpEff databases)" \
    "downloaded by snpEff (snpeff.blob.core.windows.net)" "chr-prefixed" \
    "${snpeff_db:-${SNPEFF_DIR}/hg38 (not yet downloaded)}" \
    "$([ -n "$snpeff_db" ] && echo "directory; see snpEffectPredictor.bin" || echo absent)" >> "$OUT"

echo "Wrote $OUT:"
column -t -s$'\t' "$OUT" 2>/dev/null || cat "$OUT"
echo
echo "NOTE: run this after download_refs.sh so checksums reflect the files in use."
