# ==============================================================================
# SECTION 2 – TYPE-SPECIFIC FILTERING, CLEAN VCF, STATISTICS & VISUALIZATION
#
# Rules 14a, 14b, 15, 16, 17
# All rules consume Joint_all.vcf.gz produced by merge_joint_vcfs (Rule 12b).
#
# Workflow:
#   Joint_all.vcf.gz ──┬── [14a] filter_snps   ──┐
#                      └── [14b] filter_indels  ──┴── [15] select_pass_variants
#                                                           │
#                                                   Final_PASS.vcf.gz
#                                                           │
#                                               ┌──────────┴──────────┐
#                                         [16] variant_statistics   (same input)
#                                               │
#                                         stats/*.tsv
#                                               │
#                                         [17] visualize_variants
#                                               │
#                                          plots/*.png
# ==============================================================================

BCFTOOLS_CONTAINER   = "docker://quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
PYTHON_VIZ_CONTAINER = "docker://quay.io/biocontainers/matplotlib:3.5.1"

# dbSNP paths – plain VCF from config; bgzipped + tabix version built by prepare_dbsnp_gz
DBSNP_VCF = config["reference"]["dbsnp"]
DBSNP_GZ  = DBSNP_VCF + ".gz"

# SnpEff / SnpSift tool paths (vcf_annotation conda env, Java 25)
_VCF_ANN_ENV  = "/home/iman/miniforge3/envs/vcf_annotation"
JAVA_BIN      = f"{_VCF_ANN_ENV}/lib/jvm/bin/java"
SNPEFF_JAR    = f"{_VCF_ANN_ENV}/share/snpeff-5.4.0a-0/snpEff.jar"
SNPSIFT_JAR   = f"{_VCF_ANN_ENV}/share/snpsift-5.4.0a-0/SnpSift.jar"
BCFTOOLS_BIN  = "/home/iman/miniforge3/envs/gatk_env/bin/bcftools"


# ── Rule 14pre: bgzip + tabix the dbSNP VCF (one-time, reused by annotate_rsid) ─
rule prepare_dbsnp_gz:
    """
    Compresses the plain dbSNP VCF with bgzip and creates a tabix index.
    Required once so that bcftools annotate can use it for rsID lookup.
    Output lives alongside the source VCF in resources/.
    """
    input:  DBSNP_VCF
    output:
        gz  = DBSNP_GZ,
        tbi = DBSNP_GZ + ".tbi"
    singularity: BCFTOOLS_CONTAINER
    log: f"{LOG_DIR}/prepare_dbsnp_gz.log"
    resources:
        mem_mb  = 4096,
        runtime = 120
    shell:
        """
        exec &> {log}
        bgzip -c {input} > {output.gz}
        tabix -p vcf {output.gz}
        """


# ── Rule 14a: Select SNPs and apply RNA-seq SNP hard filters ──────────────────
rule filter_snps:
    """
    Extracts SNPs from the joint-called VCF and applies GATK best-practice
    hard filters tuned for RNA-seq data.
    Thresholds are configured in config: variant_filters.snp.*
    Output is temporary; it is consumed by select_pass_variants (Rule 15).
    """
    input:
        vcf = f"{OUTPUT_DIR}/Joint_all.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Joint_all.vcf.gz.tbi"
    output:
        vcf = temp(f"{OUTPUT_DIR}/SNPs_hardfilt.vcf.gz"),
        tbi = temp(f"{OUTPUT_DIR}/SNPs_hardfilt.vcf.gz.tbi")
    singularity: GATK_CONTAINER
    log: f"{LOG_DIR}/filter_snps.log"
    params:
        reference = config["reference"]["fasta"],
        qd        = config["variant_filters"]["snp"]["QD"],
        fs        = config["variant_filters"]["snp"]["FS"],
        mq        = config["variant_filters"]["snp"]["MQ"],
        mqrs      = config["variant_filters"]["snp"]["MQRankSum"],
        rprs      = config["variant_filters"]["snp"]["ReadPosRankSum"],
        dp        = config["variant_filters"]["snp"]["DP"],
        raw       = f"{OUTPUT_DIR}/SNPs_raw.vcf.gz"
    resources:
        mem_mb  = 16384,
        runtime = 180
    shell:
        """
        exec &> {log}

        # Step 1: Extract SNPs only
        gatk --java-options "-Xmx{resources.mem_mb}M" SelectVariants \
            -R {params.reference} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O {params.raw}

        # Step 2: Apply SNP-specific hard filters (GATK best practices for RNA-seq)
        gatk --java-options "-Xmx{resources.mem_mb}M" VariantFiltration \
            -R {params.reference} \
            -V {params.raw} \
            --filter-name "SNP_LowQD"        --filter-expression "QD < {params.qd}" \
            --filter-name "SNP_HighFS"        --filter-expression "FS > {params.fs}" \
            --filter-name "SNP_LowMQ"         --filter-expression "MQ < {params.mq}" \
            --filter-name "SNP_LowMQRankSum"  --filter-expression "MQRankSum < {params.mqrs}" \
            --filter-name "SNP_LowReadPos"    --filter-expression "ReadPosRankSum < {params.rprs}" \
            --filter-name "SNP_LowDP"         --filter-expression "DP < {params.dp}" \
            -O {output.vcf}

        rm -f {params.raw} {params.raw}.tbi
        """


# ── Rule 14b: Select INDELs and apply RNA-seq INDEL hard filters ──────────────
rule filter_indels:
    """
    Extracts INDELs from the joint-called VCF and applies GATK best-practice
    hard filters tuned for RNA-seq data.
    Thresholds are configured in config: variant_filters.indel.*
    Output is temporary; it is consumed by select_pass_variants (Rule 15).
    """
    input:
        vcf = f"{OUTPUT_DIR}/Joint_all.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Joint_all.vcf.gz.tbi"
    output:
        vcf = temp(f"{OUTPUT_DIR}/INDELs_hardfilt.vcf.gz"),
        tbi = temp(f"{OUTPUT_DIR}/INDELs_hardfilt.vcf.gz.tbi")
    singularity: GATK_CONTAINER
    log: f"{LOG_DIR}/filter_indels.log"
    params:
        reference = config["reference"]["fasta"],
        qd        = config["variant_filters"]["indel"]["QD"],
        fs        = config["variant_filters"]["indel"]["FS"],
        rprs      = config["variant_filters"]["indel"]["ReadPosRankSum"],
        dp        = config["variant_filters"]["indel"]["DP"],
        raw       = f"{OUTPUT_DIR}/INDELs_raw.vcf.gz"
    resources:
        mem_mb  = 16384,
        runtime = 180
    shell:
        """
        exec &> {log}

        # Step 1: Extract INDELs only
        gatk --java-options "-Xmx{resources.mem_mb}M" SelectVariants \
            -R {params.reference} \
            -V {input.vcf} \
            --select-type-to-include INDEL \
            -O {params.raw}

        # Step 2: Apply INDEL-specific hard filters (GATK best practices for RNA-seq)
        gatk --java-options "-Xmx{resources.mem_mb}M" VariantFiltration \
            -R {params.reference} \
            -V {params.raw} \
            --filter-name "INDEL_LowQD"      --filter-expression "QD < {params.qd}" \
            --filter-name "INDEL_HighFS"      --filter-expression "FS > {params.fs}" \
            --filter-name "INDEL_LowReadPos"  --filter-expression "ReadPosRankSum < {params.rprs}" \
            --filter-name "INDEL_LowDP"       --filter-expression "DP < {params.dp}" \
            -O {output.vcf}

        rm -f {params.raw} {params.raw}.tbi
        """


# ── Rule 15: Merge type-filtered VCFs and select PASS-only variants ───────────
rule select_pass_variants:
    """
    Merges the SNP-filtered and INDEL-filtered VCFs, then discards all variants
    carrying any filter flag, producing the final clean PASS-only VCF.
    This is the terminal variant-calling output of the pipeline.
    """
    input:
        snps     = f"{OUTPUT_DIR}/SNPs_hardfilt.vcf.gz",
        snps_tbi = f"{OUTPUT_DIR}/SNPs_hardfilt.vcf.gz.tbi",
        indels     = f"{OUTPUT_DIR}/INDELs_hardfilt.vcf.gz",
        indels_tbi = f"{OUTPUT_DIR}/INDELs_hardfilt.vcf.gz.tbi"
    output:
        vcf = f"{OUTPUT_DIR}/Final_PASS.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Final_PASS.vcf.gz.tbi"
    singularity: GATK_CONTAINER
    log: f"{LOG_DIR}/select_pass_variants.log"
    params:
        reference = config["reference"]["fasta"],
        merged    = f"{OUTPUT_DIR}/Joint_typed_filtered.vcf.gz"
    resources:
        mem_mb  = 16384,
        runtime = 120
    shell:
        """
        exec &> {log}

        # Step 1: Merge SNP and INDEL filtered VCFs
        gatk --java-options "-Xmx{resources.mem_mb}M" MergeVcfs \
            -I {input.snps} \
            -I {input.indels} \
            -O {params.merged}

        # Step 2: Remove all variants with a filter flag (keep PASS only)
        gatk --java-options "-Xmx{resources.mem_mb}M" SelectVariants \
            -R {params.reference} \
            -V {params.merged} \
            --exclude-filtered \
            -O {output.vcf}

        gatk IndexFeatureFile -I {output.vcf}

        rm -f {params.merged} {params.merged}.tbi
        """


# ── Rule 15b: Remove known RNA editing sites (REDIportal) ────────────────────
rule remove_rna_editing:
    """
    Removes known RNA A-to-I editing sites from the PASS VCF using the
    REDIportal hg38 database.  bcftools -T with the '^' prefix streams
    through the VCF and excludes any variant whose position overlaps an
    entry in the BED file (0-based coordinates handled automatically).
    """
    input:
        vcf = f"{OUTPUT_DIR}/Final_PASS.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Final_PASS.vcf.gz.tbi",
        bed = config["rna_editing"]["rediportal_bed"]
    output:
        vcf = f"{OUTPUT_DIR}/Final_PASS_noEdit.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Final_PASS_noEdit.vcf.gz.tbi"
    singularity: BCFTOOLS_CONTAINER
    log: f"{LOG_DIR}/remove_rna_editing.log"
    resources:
        mem_mb  = 8192,
        runtime = 90
    shell:
        """
        exec &> {log}

        # -T ^bed streams through VCF excluding positions in the BED file.
        # bcftools auto-detects 0-based BED coords from the .bed extension.
        bcftools view \
            -T ^{input.bed} \
            {input.vcf} \
            -O z -o {output.vcf}

        bcftools index -t {output.vcf}
        """


# ── Rule 15c: Remove dense variant clusters (RNA-seq artefacts) ──────────────
rule filter_clusters:
    """
    Tags and removes dense variant clusters (>= cluster_size variants
    within window_size bp) using GATK VariantFiltration.  Clusters are
    common artefacts at mis-spliced junctions and repetitive regions in
    RNA-seq data.
    """
    input:
        vcf = f"{OUTPUT_DIR}/Final_PASS_noEdit.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Final_PASS_noEdit.vcf.gz.tbi"
    output:
        vcf = f"{OUTPUT_DIR}/Final_clean.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Final_clean.vcf.gz.tbi"
    singularity: GATK_CONTAINER
    log: f"{LOG_DIR}/filter_clusters.log"
    params:
        reference = config["reference"]["fasta"],
        window    = config["cluster_filter"]["window_size"],
        size      = config["cluster_filter"]["cluster_size"],
        tagged    = f"{OUTPUT_DIR}/Final_PASS_noEdit_tagged.vcf.gz"
    resources:
        mem_mb  = 16384,
        runtime = 60
    shell:
        """
        exec &> {log}

        # Step 1: tag clustered variants
        gatk --java-options "-Xmx{resources.mem_mb}M" VariantFiltration \
            -R {params.reference} \
            -V {input.vcf} \
            --cluster-window-size {params.window} \
            --cluster-size {params.size} \
            -O {params.tagged}

        # Step 2: keep only PASS (untagged) variants
        gatk --java-options "-Xmx{resources.mem_mb}M" SelectVariants \
            -R {params.reference} \
            -V {params.tagged} \
            --exclude-filtered \
            -O {output.vcf}

        rm -f {params.tagged} {params.tagged}.tbi
        """


# ── Rule 15d: Annotate variants with SnpEff ───────────────────────────────────
rule annotate_variants:
    """
    Annotates all variants in Final_clean.vcf.gz with SnpEff (hg38), adding
    ANN fields describing the functional effect of each variant on every
    overlapping transcript.  The annotated VCF is the input for the CDS filter
    and the SnpEff HTML health-check report.
    """
    input:
        vcf = f"{OUTPUT_DIR}/Final_clean.vcf.gz"
    output:
        vcf    = f"{OUTPUT_DIR}/Final_annotated.vcf.gz",
        tbi    = f"{OUTPUT_DIR}/Final_annotated.vcf.gz.tbi",
        report = f"{OUTPUT_DIR}/Final_annotated_healthCheck.html"
    log: f"{LOG_DIR}/annotate_variants.log"
    params:
        genome   = config["snpeff"]["genome"],
        data_dir = config["snpeff"]["data_dir"]
    resources:
        mem_mb  = 16384,
        runtime = 120
    shell:
        """
        exec &> {log}

        # Annotate and compress in one pipe; write HTML report alongside
        {JAVA_BIN} -Xmx{resources.mem_mb}m \
            -jar {SNPEFF_JAR} \
            -v {params.genome} \
            -s {output.report} \
            -dataDir {params.data_dir} \
            {input.vcf} \
        | {BCFTOOLS_BIN} view -O z -o {output.vcf}

        {BCFTOOLS_BIN} index -t {output.vcf}
        """


# ── Rule 15e: CDS filter – keep only coding-region variants ──────────────────
rule filter_cds:
    """
    Retains only variants with a coding consequence (missense, synonymous,
    frameshift, stop/start changes, inframe indels) using SnpSift.
    Removes intronic, UTR, intergenic, and other non-coding MODIFIER variants
    that dominate Ribo-Zero total-RNA datasets.
    This is the terminal output of the variant-calling pipeline.
    """
    input:
        vcf = f"{OUTPUT_DIR}/Final_annotated.vcf.gz"
    output:
        vcf = f"{OUTPUT_DIR}/Final_CDS.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Final_CDS.vcf.gz.tbi"
    log: f"{LOG_DIR}/filter_cds.log"
    resources:
        mem_mb  = 8192,
        runtime = 30
    shell:
        """
        exec &> {log}

        {JAVA_BIN} -Xmx{resources.mem_mb}m \
            -jar {SNPSIFT_JAR} filter \
            "( ANN[*].EFFECT has 'missense_variant'       ) | \
             ( ANN[*].EFFECT has 'synonymous_variant'     ) | \
             ( ANN[*].EFFECT has 'stop_gained'            ) | \
             ( ANN[*].EFFECT has 'stop_lost'              ) | \
             ( ANN[*].EFFECT has 'start_lost'             ) | \
             ( ANN[*].EFFECT has 'frameshift_variant'     ) | \
             ( ANN[*].EFFECT has 'inframe_insertion'      ) | \
             ( ANN[*].EFFECT has 'inframe_deletion'       )" \
            {input.vcf} \
        | {BCFTOOLS_BIN} view -O z -o {output.vcf}

        {BCFTOOLS_BIN} index -t {output.vcf}
        """


# ── Rule 15g: Annotate rsIDs from dbSNP138 ───────────────────────────────────
rule annotate_rsid:
    """
    Adds dbSNP rsIDs to the ID column of Final_CDS.vcf.gz using bcftools annotate.
    Matches variants by chromosome + position + ref + alt against dbSNP138 hg38.
    Produces Final_CDS_rsID.vcf.gz — the final pipeline output used for reporting.
    """
    input:
        vcf       = f"{OUTPUT_DIR}/Final_CDS.vcf.gz",
        tbi       = f"{OUTPUT_DIR}/Final_CDS.vcf.gz.tbi",
        dbsnp     = DBSNP_GZ,
        dbsnp_tbi = DBSNP_GZ + ".tbi"
    output:
        vcf = f"{OUTPUT_DIR}/Final_CDS_rsID.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Final_CDS_rsID.vcf.gz.tbi"
    singularity: BCFTOOLS_CONTAINER
    log: f"{LOG_DIR}/annotate_rsid.log"
    resources:
        mem_mb  = 8192,
        runtime = 60
    shell:
        """
        exec &> {log}

        bcftools annotate \
            -a {input.dbsnp} \
            -c ID \
            -O z -o {output.vcf} \
            {input.vcf}

        bcftools index -t {output.vcf}
        """


# ── Rule 15f: Health check – SnpEff summary report on Final_CDS ──────────────
rule health_check_cds:
    """
    Runs SnpEff on the final CDS-filtered VCF to generate an HTML summary
    report (variant counts, effect categories, Ts/Tv, codon usage, etc.).
    The annotated VCF output is discarded — only the HTML report is kept.
    This mirrors the health-check pattern used in the annotation step.
    """
    input:
        vcf = f"{OUTPUT_DIR}/Final_CDS_rsID.vcf.gz"
    output:
        report = f"{OUTPUT_DIR}/Final_CDS_healthCheck.html"
    log: f"{LOG_DIR}/health_check_cds.log"
    params:
        genome   = config["snpeff"]["genome"],
        data_dir = config["snpeff"]["data_dir"]
    resources:
        mem_mb  = 16384,
        runtime = 60
    shell:
        """
        exec &> {log}

        {JAVA_BIN} -Xmx{resources.mem_mb}m \
            -jar {SNPEFF_JAR} \
            -v {params.genome} \
            -s {output.report} \
            -dataDir {params.data_dir} \
            {input.vcf} > /dev/null
        """


# ── Rule 16: Generate per-variant statistics TSVs via bcftools ────────────────
rule variant_statistics:
    """
    Extracts quality metrics and counts from the pre-filter (Joint_all) and
    post-filter (Final_PASS) VCFs, writing flat TSV files consumed by
    visualize_variants (Rule 17).

    Outputs:
        quality_metrics.tsv  – TYPE, FILTER, QD, FS, MQ, DP per variant (pre-filter)
        filter_summary.tsv   – count of variants per filter flag
        chrom_counts.tsv     – PASS variant counts per chromosome and type
        sample_counts.tsv    – per-sample SNP / INDEL counts (PASS only)
        ts_tv.tsv            – Ts, Tv counts and ratio before and after filtering
    """
    input:
        raw_vcf  = f"{OUTPUT_DIR}/Joint_all.vcf.gz",
        pass_vcf = f"{OUTPUT_DIR}/Final_CDS.vcf.gz"
    output:
        metrics_tsv = f"{OUTPUT_DIR}/stats/quality_metrics.tsv",
        filter_tsv  = f"{OUTPUT_DIR}/stats/filter_summary.tsv",
        chrom_tsv   = f"{OUTPUT_DIR}/stats/chrom_counts.tsv",
        sample_tsv  = f"{OUTPUT_DIR}/stats/sample_counts.tsv",
        tstv_tsv    = f"{OUTPUT_DIR}/stats/ts_tv.tsv"
    singularity: BCFTOOLS_CONTAINER
    log: f"{LOG_DIR}/variant_statistics.log"
    resources:
        mem_mb  = 8192,
        runtime = 60
    shell:
        """
        set +o pipefail
        exec &> {log}
        mkdir -p {OUTPUT_DIR}/stats

        # Quality metrics for all variants (before PASS filter)
        bcftools query \
            -f '%TYPE\t%FILTER\t%INFO/QD\t%INFO/FS\t%INFO/MQ\t%INFO/DP\n' \
            {input.raw_vcf} \
            > {output.metrics_tsv}

        # Count PASS and Filtered variants per type.
        # Joint_all.vcf.gz is pre-filter so its FILTER column is always '.'.
        # Instead, derive counts by comparing raw vs Final_PASS VCF.
        for entry in "snps SNP" "indels INDEL"; do
            bcftools_type=$(echo $entry | cut -d' ' -f1)
            display_type=$(echo $entry | cut -d' ' -f2)
            raw_count=$(bcftools view -H -v "$bcftools_type" {input.raw_vcf} | wc -l)
            pass_count=$(bcftools view -H -v "$bcftools_type" {input.pass_vcf} | wc -l)
            filtered=$((raw_count - pass_count))
            printf '%s\tPASS\t%d\n%s\tFiltered\t%d\n' \
                "$display_type" "$pass_count" "$display_type" "$filtered"
        done > {output.filter_tsv}

        # Per-chromosome variant counts for clean (post all filters) variants
        bcftools query -f '%CHROM\t%TYPE\n' {input.pass_vcf} \
            | sort | uniq -c \
            | awk '{{print $2"\t"$3"\t"$1}}' \
            > {output.chrom_tsv}

        # Per-sample SNP (ts+tv) and INDEL counts from bcftools stats PSC block
        bcftools stats -s - {input.pass_vcf} \
            | awk '/^PSC/ {{print $3"\t"($7+$8)"\t"$9}}' \
            > {output.sample_tsv}

        # Ts/Tv ratio before and after filtering
        # Format: stage  ts  tv  ratio
        echo -e "stage\tts\ttv\tratio" > {output.tstv_tsv}
        bcftools stats {input.raw_vcf} \
            | awk '/^TSTV/ {{print "pre_filter\t"$3"\t"$4"\t"$5; exit}}' \
            >> {output.tstv_tsv}
        bcftools stats {input.pass_vcf} \
            | awk '/^TSTV/ {{print "post_filter\t"$3"\t"$4"\t"$5; exit}}' \
            >> {output.tstv_tsv}
        """


# ── Rule 17: Generate publication-quality figures ─────────────────────────────
rule visualize_variants:
    """
    Reads the TSV files produced by variant_statistics and writes six
    publication-quality PNG figures to results/plots/:

        variant_type_summary.png   – SNP / INDEL counts, PASS vs filtered
        filter_summary.png         – horizontal bar of filter flag frequencies
        quality_distributions.png  – QD, FS, MQ, DP histograms per variant type
        chrom_distribution.png     – per-chromosome variant counts (PASS only)
        sample_counts.png          – per-sample SNP / INDEL counts (PASS only)
        ts_tv_comparison.png       – Ts/Tv counts and ratio before vs after filtering
    """
    input:
        metrics_tsv = f"{OUTPUT_DIR}/stats/quality_metrics.tsv",
        filter_tsv  = f"{OUTPUT_DIR}/stats/filter_summary.tsv",
        chrom_tsv   = f"{OUTPUT_DIR}/stats/chrom_counts.tsv",
        sample_tsv  = f"{OUTPUT_DIR}/stats/sample_counts.tsv",
        tstv_tsv    = f"{OUTPUT_DIR}/stats/ts_tv.tsv"
    output:
        directory(f"{OUTPUT_DIR}/plots")
    singularity: PYTHON_VIZ_CONTAINER
    log: f"{LOG_DIR}/visualize_variants.log"
    params:
        script = "workflow/scripts/visualize_variants.py"
    resources:
        mem_mb  = 8192,
        runtime = 30
    shell:
        """
        exec &> {log}
        mkdir -p {output}
        python {params.script} \
            --metrics  {input.metrics_tsv} \
            --filters  {input.filter_tsv} \
            --chroms   {input.chrom_tsv} \
            --samples  {input.sample_tsv} \
            --tstv     {input.tstv_tsv} \
            --outdir   {output}
        """
