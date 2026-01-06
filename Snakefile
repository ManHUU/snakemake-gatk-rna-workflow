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


# 3. Sample Discovery Logic
# This runs once when snakemake starts
FOUND_FASTQ = glob.glob(os.path.join(FASTQ_DIR, "*.fastq"))
SAMPLES = sorted([os.path.basename(f).split('.')[0] for f in FOUND_FASTQ])


# 4. Target Rule
rule all:
    input:
        REF_FASTA + ".fai",
        os.path.join(OUTPUT_DIR, "Filtered.vcf.gz")


# 5. Include the "Worker" rules
include: "workflow/rules/variant_calling_jointCall_per_Chr.smk"


# Optional: On-start/On-error messaging
onsuccess:
    print("Workflow finished successfully!")

onerror:
    print("Workflow failed. Check logs in the logs/ directory.")
