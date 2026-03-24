import os
import glob


# 1. Load the configurations
# We point to the specific paths in your workflow directory
configfile: "config/config.yaml"


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

# 4. Sample Discovery Logic – file naming depends on sequencing type:
#   Paired-end  → {sample}_1.fastq  +  {sample}_2.fastq
#   Single-end  → {sample}.fastq
if IS_PAIRED:
    FOUND_R1 = glob.glob(os.path.join(FASTQ_DIR, "*_1.fastq"))
    SAMPLES = sorted([os.path.basename(f).replace("_1.fastq", "") for f in FOUND_R1])
else:
    FOUND_SE = glob.glob(os.path.join(FASTQ_DIR, "*.fastq"))
    SAMPLES = sorted([os.path.basename(f).replace(".fastq", "") for f in FOUND_SE])


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
