import os
import csv
from snakemake.utils import validate


# 1. Load the configurations
# We point to the specific paths in your workflow directory
configfile: "config/config.yaml"
# Fail fast at DAG-build time on a mistyped/missing config field.
validate(config, "workflow/schemas/config.schema.yaml")


# --- CONTAINER DEFINITIONS ---
# STAR Container (Updated to your specific version)
# Was: 2.7.10b -> Now: 2.7.11b
STAR_CONTAINER = "docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_3"

# Samtools Container
SAMTOOLS_CONTAINER = "docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0"

# GATK Container
GATK_CONTAINER = "docker://broadinstitute/gatk:4.6.1.0"


# -----------------------------

# 2. Path & Variable Definitions
OUTPUT_DIR = config["output_dir"]
LOG_DIR = config["log_dir"]
FASTQ_DIR = config["fastq_dir"]

# Reference files are in root/resources (defined in config.yaml)
REF_FASTA = config["reference"]["fasta"]
REF_FAI = REF_FASTA + ".fai"
REF_DICT = os.path.splitext(REF_FASTA)[0] + ".dict"
CHROMOSOMES = config["chromosomes"]


# 3. Sequencing mode (set in config: sequencing_type: "paired" | "single")
IS_PAIRED = config.get("sequencing_type", "paired") == "paired"

# 4. Sample sheet – samples and FASTQ paths come from config/units.tsv
#   (columns: sample, unit, fq1, fq2). This replaces runtime globbing so the
#   exact set of analysed samples is explicit, version-controlled and
#   schema-validated. One row per sample (the `unit` column is kept for
#   snakemake-workflows template compatibility; multi-lane merging is not wired
#   in). Regenerate the sheet with workflow/scripts/make_units.sh.
UNITS_TSV = config.get("units", "config/units.tsv")

UNITS = {}
with open(UNITS_TSV, newline="") as _fh:
    for _row in csv.DictReader(_fh, delimiter="\t"):
        _row = {k: (v.strip() if isinstance(v, str) else v) for k, v in _row.items()}
        validate(_row, "workflow/schemas/units.schema.yaml")
        _sample = _row["sample"]
        if _sample in UNITS:
            raise ValueError(
                f"Duplicate sample '{_sample}' in {UNITS_TSV}: this pipeline "
                "expects one row per sample (multi-lane merging is not wired in)."
            )
        if IS_PAIRED and not _row.get("fq2"):
            raise ValueError(
                f"sequencing_type is 'paired' but sample '{_sample}' has no fq2 "
                f"in {UNITS_TSV}."
            )
        UNITS[_sample] = _row

SAMPLES = sorted(UNITS)


# FASTQ path lookups (sample sheet is the single source of truth).
def fq1_for(sample):
    return UNITS[sample]["fq1"]


def fq2_for(sample):
    return UNITS[sample]["fq2"]


# 5. Target Rule
rule all:
    input:
        # Reference indices
        REF_FASTA + ".fai",
        # QC
        os.path.join(OUTPUT_DIR, "qc", "multiqc_report.html"),
        # Variant calling – final CDS-filtered VCF with rsIDs
        os.path.join(OUTPUT_DIR, "Final_CDS_rsID.vcf.gz"),
        # SnpEff health-check reports
        os.path.join(OUTPUT_DIR, "Final_annotated_healthCheck.html"),
        os.path.join(OUTPUT_DIR, "Final_CDS_healthCheck.html"),
        # Statistics
        os.path.join(OUTPUT_DIR, "stats", "quality_metrics.tsv"),
        os.path.join(OUTPUT_DIR, "stats", "ts_tv.tsv"),
        # Visualization
        os.path.join(OUTPUT_DIR, "plots")


# 6. Include QC rules
include: "workflow/rules/qc.smk"

# 7. Include the variant-calling rules
include: "workflow/rules/variant_calling_jointCall_per_Chr.smk"

# 8. Include filtering, statistics and visualization rules
include: "workflow/rules/filter_and_visualize.smk"


# Optional: On-start/On-error messaging
onsuccess:
    print("Workflow finished successfully!")

onerror:
    print("Workflow failed. Check logs in the logs/ directory.")
