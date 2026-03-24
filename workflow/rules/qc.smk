# ==============================================================================
# SECTION 0 – QUALITY CONTROL (QC)
#
# Rules QC1, QC2
# Runs FastQC on raw FASTQ inputs and aggregates results with MultiQC.
# These rules run in parallel with the alignment/variant-calling pipeline.
#
# Outputs:
#   results/qc/fastqc/{sample}_fastqc.html   – per-sample FastQC report
#   results/qc/fastqc/{sample}_fastqc.zip    – per-sample FastQC data
#   results/qc/multiqc_report.html           – aggregated MultiQC report
#
# MultiQC aggregates:
#   - FastQC reports (raw read quality)
#   - STAR alignment logs (mapping rate, mismatch rate, splice junctions)
#   - GATK MarkDuplicates metrics (duplication rate)
# ==============================================================================

FASTQC_CONTAINER  = "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
MULTIQC_CONTAINER = "docker://quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0"


# ── Rule QC1: FastQC on raw FASTQ files ───────────────────────────────────────
# Two versions are defined: one for paired-end, one for single-end.
# Snakemake loads only the branch matching IS_PAIRED (set from config).
#
# Paired-end  (sequencing_type: "paired"):
#   Input  → {sample}_1.fastq  +  {sample}_2.fastq
#   Output → {sample}_1_fastqc.{html,zip}  +  {sample}_2_fastqc.{html,zip}
#
# Single-end  (sequencing_type: "single"):
#   Input  → {sample}.fastq
#   Output → {sample}_fastqc.{html,zip}
# ─────────────────────────────────────────────────────────────────────────────

if IS_PAIRED:
    rule fastqc:
        """Runs FastQC on R1 and R2 for each paired-end sample."""
        input:
            r1 = lambda wildcards: os.path.join(FASTQ_DIR, f"{wildcards.sample}_1.fastq"),
            r2 = lambda wildcards: os.path.join(FASTQ_DIR, f"{wildcards.sample}_2.fastq")
        output:
            html_r1 = f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_1_fastqc.html",
            zip_r1  = f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_1_fastqc.zip",
            html_r2 = f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_2_fastqc.html",
            zip_r2  = f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_2_fastqc.zip"
        singularity: FASTQC_CONTAINER
        log: f"{LOG_DIR}/{{sample}}.fastqc.log"
        threads: 4
        resources:
            mem_mb  = 4096,
            runtime = 60
        shell:
            """
            exec &> {log}
            mkdir -p {OUTPUT_DIR}/qc/fastqc
            fastqc --threads {threads} --outdir {OUTPUT_DIR}/qc/fastqc {input.r1} {input.r2}
            """

    # FastQC ZIP paths expected by MultiQC (one per read per sample)
    FASTQC_ZIPS = expand(
        [f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_1_fastqc.zip",
         f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_2_fastqc.zip"],
        sample=SAMPLES
    )
    # fastp JSON paths for MultiQC (one per sample)
    FASTP_JSONS = expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_fastp.json", sample=SAMPLES)

else:
    rule fastqc:
        """Runs FastQC on the single FASTQ file for each single-end sample."""
        input:
            r1 = lambda wildcards: os.path.join(FASTQ_DIR, f"{wildcards.sample}.fastq")
        output:
            html = f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_fastqc.html",
            zip  = f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_fastqc.zip"
        singularity: FASTQC_CONTAINER
        log: f"{LOG_DIR}/{{sample}}.fastqc.log"
        threads: 4
        resources:
            mem_mb  = 4096,
            runtime = 60
        shell:
            """
            exec &> {log}
            mkdir -p {OUTPUT_DIR}/qc/fastqc
            fastqc --threads {threads} --outdir {OUTPUT_DIR}/qc/fastqc {input.r1}
            """

    # FastQC ZIP paths expected by MultiQC (one per sample)
    FASTQC_ZIPS = expand(
        f"{OUTPUT_DIR}/qc/fastqc/{{sample}}_fastqc.zip",
        sample=SAMPLES
    )
    # fastp JSON paths for MultiQC (one per sample)
    FASTP_JSONS = expand(f"{OUTPUT_DIR}/trimmed/{{sample}}_fastp.json", sample=SAMPLES)


# ── Rule QC2: MultiQC – aggregate FastQC + STAR logs + MarkDuplicates ────────
rule multiqc:
    """
    Aggregates quality reports from all samples into a single interactive HTML
    report using MultiQC.

    Inputs aggregated:
        - FastQC per-sample ZIP archives          (raw read quality)
        - STAR *Log.final.out files               (alignment rate, mismatches)
        - GATK MarkDuplicates metrics TXT files   (duplication rate)
    """
    input:
        fastqc_zips = FASTQC_ZIPS,
        fastp_jsons = FASTP_JSONS,
        star_logs = expand(
            f"{OUTPUT_DIR}/{{sample}}.Log.final.out",
            sample=SAMPLES
        ),
        markdup_metrics = expand(
            f"{OUTPUT_DIR}/{{sample}}.markdup_metrics.txt",
            sample=SAMPLES
        )
    output:
        html = f"{OUTPUT_DIR}/qc/multiqc_report.html",
        data = directory(f"{OUTPUT_DIR}/qc/multiqc_report_data")
    singularity: MULTIQC_CONTAINER
    log: f"{LOG_DIR}/multiqc.log"
    resources:
        mem_mb  = 4096,
        runtime = 30
    shell:
        """
        exec &> {log}
        multiqc \
            {OUTPUT_DIR}/qc/fastqc \
            {OUTPUT_DIR}/trimmed \
            {OUTPUT_DIR} \
            --filename multiqc_report.html \
            --outdir {OUTPUT_DIR}/qc \
            --force
        """
