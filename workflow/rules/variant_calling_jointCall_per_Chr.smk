# GATK RNA-seq Variant Calling Pipeline - Snakemake Workflow

# Cluster wall-time cap (minutes) for the attempt-scaled rules. This is a
# per-site value and lives in config/config.yaml (`partition_max_runtime`) so
# readers never edit pipeline code — set it to your SLURM partition's MaxTime
# (`scontrol show partition <name>`). Default 7200 = 5 days (Snellius `genoa`).
# The attempt-scaled rules cap their retry runtime at this value so the third
# retry cannot exceed the partition limit and get rejected at submission.
PARTITION_MAX_RUNTIME = config.get("partition_max_runtime", 7200)

FASTP_CONTAINER = "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"

# Rule Pre: Adapter Trimming with fastp
# Trims adapters and low-quality bases before alignment.
# fastp auto-detects Illumina adapters for paired-end data.
# Output trimmed FASTQs are consumed by star_align (Rule 1).
if IS_PAIRED:
    rule fastp_trim:
        input:
            r1 = lambda wildcards: fq1_for(wildcards.sample),
            r2 = lambda wildcards: fq2_for(wildcards.sample)
        output:
            r1   = f"{OUTPUT_DIR}/trimmed/{{sample}}_1.trimmed.fastq",
            r2   = f"{OUTPUT_DIR}/trimmed/{{sample}}_2.trimmed.fastq",
            html = f"{OUTPUT_DIR}/trimmed/{{sample}}_fastp.html",
            json = f"{OUTPUT_DIR}/trimmed/{{sample}}_fastp.json"
        singularity: FASTP_CONTAINER
        log: f"{LOG_DIR}/{{sample}}.fastp.log"
        threads: 8
        resources:
            mem_mb=8192,
            runtime=120
        shell:
            """
            exec &>> {log}
            mkdir -p {OUTPUT_DIR}/trimmed
            fastp \
                --in1 {input.r1} \
                --in2 {input.r2} \
                --out1 {output.r1} \
                --out2 {output.r2} \
                --html {output.html} \
                --json {output.json} \
                --thread {threads} \
                --detect_adapter_for_pe \
                --trim_poly_g \
                --length_required 36
            """
else:
    rule fastp_trim:
        input:
            r1 = lambda wildcards: fq1_for(wildcards.sample)
        output:
            r1   = f"{OUTPUT_DIR}/trimmed/{{sample}}.trimmed.fastq",
            html = f"{OUTPUT_DIR}/trimmed/{{sample}}_fastp.html",
            json = f"{OUTPUT_DIR}/trimmed/{{sample}}_fastp.json"
        singularity: FASTP_CONTAINER
        log: f"{LOG_DIR}/{{sample}}.fastp.log"
        threads: 8
        resources:
            mem_mb=8192,
            runtime=120
        shell:
            """
            exec &>> {log}
            mkdir -p {OUTPUT_DIR}/trimmed
            fastp \
                --in1 {input.r1} \
                --out1 {output.r1} \
                --html {output.html} \
                --json {output.json} \
                --thread {threads} \
                --trim_poly_g \
                --length_required 36
            """


# Rule 0a: Prepare FAI (Uses Samtools)
rule prepare_fasta_index:
    input: REF_FASTA
    output: REF_FAI
    singularity: SAMTOOLS_CONTAINER
    log: f"{LOG_DIR}/prepare_fai.log"
    resources:
        mem_mb=8192,
        runtime=60
    shell:
        "samtools faidx {input}"


# Rule 0b: Prepare Dictionary (Uses GATK)
rule prepare_dict:
    input: REF_FASTA
    output: REF_DICT
    singularity: GATK_CONTAINER
    log: f"{LOG_DIR}/prepare_dict.log"
    resources:
        mem_mb=8192,
        runtime=60
    shell:
        "gatk CreateSequenceDictionary -R {input}"


# Rule 0c: Build STAR Index
# This runs automatically if the "star_hg38_index" folder is missing.
rule build_star_index:
    input:
        fasta = REF_FASTA,
        gtf = config["reference"]["gtf"]
    output:
        # We use directory() to tell Snakemake this creates a folder
        directory(config["reference"]["star_index_dir"])
    singularity: STAR_CONTAINER
    log: f"{LOG_DIR}/build_star_index.log"
    threads: 16
    resources:
        # Auto-escalating on retry: 64 -> 128 -> 192 GB. STAR hg38 indexing needs
        # ~30GB+; the higher attempts cover RAM-constrained nodes on other HPCs.
        mem_mb  = lambda wildcards, attempt: 64000 * attempt,
        runtime = lambda wildcards, attempt: min(3600 * attempt, PARTITION_MAX_RUNTIME)
    shell:
        """
        exec &>> {log}
        
        # Create the output directory if it doesn't exist
        mkdir -p {output}

        # Run STAR in genomeGenerate mode
        # --sjdbOverhang 100 is standard for read lengths of 100bp+
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 149
        """


# Rule 1: Align FASTQ to Reference and Sort
# Supports both paired-end (IS_PAIRED=True) and single-end (IS_PAIRED=False).
# Paired-end  → reads = [{sample}_1.fastq, {sample}_2.fastq]
# Single-end  → reads = [{sample}.fastq]
# STAR accepts --readFilesIn R1 [R2], so a space-separated list works for both.
rule star_align:
    input:
        reads = lambda wildcards: (
            [os.path.join(OUTPUT_DIR, "trimmed", f"{wildcards.sample}_1.trimmed.fastq"),
             os.path.join(OUTPUT_DIR, "trimmed", f"{wildcards.sample}_2.trimmed.fastq")]
            if IS_PAIRED else
            [os.path.join(OUTPUT_DIR, "trimmed", f"{wildcards.sample}.trimmed.fastq")]
        ),
        star_index = config["reference"]["star_index_dir"]
    output:
        bam      = f"{OUTPUT_DIR}/{{sample}}.sorted.bam",
        star_log = f"{OUTPUT_DIR}/{{sample}}.Log.final.out"
    singularity: STAR_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.star_align.log"
    params:
        sample_id  = "{sample}",
        out_prefix = f"{OUTPUT_DIR}/{{sample}}."
    threads: 24
    resources:
        # 48GB base — STAR + hg38 needs >32GB to load the SA index plus per-thread
        # buffers and sort headroom. Auto-escalating on retry: 48 -> 96 -> 144 GB.
        mem_mb  = lambda wildcards, attempt: 49152 * attempt,
        runtime = lambda wildcards, attempt: min(600 * attempt, PARTITION_MAX_RUNTIME)
    shell:
        """
        exec &>> {log}
        # Generates coordinate-sorted BAM and adds Read Groups required for GATK.
        # {input.reads} expands to "R1 R2" (paired) or "R1" (single-end).
        #
        # STAR's BAM-sort writes ~20-40 GB of intermediate bin files to its
        # tmp dir. Default location is <outFileNamePrefix>_STARtmp/ next to the
        # output BAM — i.e. on the project filesystem. On HPC that filesystem is
        # typically quota-bound, and a mid-sort write truncation manifests as
        # "BAM bin size does not agree with size on disk". We redirect the tmp
        # dir to Snakemake's `tmpdir` resource, which resolves to node-local
        # /scratch-local/<user>.<jobid> under SLURM and /tmp locally. STAR
        # requires --outTmpDir to not already exist.
        STAR_TMP="{resources.tmpdir}/{wildcards.sample}_STARtmp"
        rm -rf "$STAR_TMP"
        STAR --runThreadN {threads} \
             --genomeDir {input.star_index} \
             --readFilesIn {input.reads} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattrRGline ID:{params.sample_id} SM:{params.sample_id} LB:lib1 PL:ILLUMINA PU:unit1 \
             --outFileNamePrefix {params.out_prefix} \
             --outTmpDir "$STAR_TMP"

        # Rename STAR's default output name to match Snakemake's output
        mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam}
        mv {params.out_prefix}Log.final.out {output.star_log}
        """


# Rule 1b: Index the Aligned BAM (Uses Samtools Container)
# This rule picks up where Rule 1 left off.
rule index_sorted_bam:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.sorted.bam"
    output:
        bai=f"{OUTPUT_DIR}/{{sample}}.sorted.bam.bai"
    singularity: SAMTOOLS_CONTAINER
    log: f"{LOG_DIR}/{{sample}}.index_star_bam.log"
    threads: 4
    resources:
        mem_mb=4096,
        runtime=60
    shell:
        """
        exec &>> {log}
        samtools index -@ {threads} {input.bam} {output.bai}
        """


# Rule 2: Mark Duplicates
rule mark_duplicates:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.sorted.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.sorted.bam.bai"
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.markdup.bam.bai",
        metrics=f"{OUTPUT_DIR}/{{sample}}.markdup_metrics.txt"
    singularity:
        GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.mark_duplicates.log"
    threads: 4
    resources:
        # Auto-escalating on retry: attempt 1 = 80GB, attempt 2 = 160GB, attempt 3 = 240GB.
        # MarkDuplicates holds duplicate indices in heap during the mark-pass; on >300M-record
        # BAMs this can blow past 60GB. The retry pattern lets the first attempt cover typical
        # libraries cheaply while outliers auto-recover. Java heap = mem_mb-2GB (see -Xmx below).
        mem_mb  = lambda wildcards, attempt: 81920 * attempt,
        runtime = lambda wildcards, attempt: min(1800 * attempt, PARTITION_MAX_RUNTIME)
    shell:
        """
        exec &>> {log}
        echo "===== $(date '+%F %T') | mark_duplicates sample={wildcards.sample} | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT
        if [ -f "{output.bam}".bai ] && [ ! -f {output.bai} ]; then
            mv "{output.bam}".bai {output.bai}
        # Handle the case where GATK strips the extension (e.g. file.bai instead of file.bam.bai)
        elif [ -f "{OUTPUT_DIR}/{wildcards.sample}.markdup.bai" ] && [ ! -f {output.bai} ]; then
             mv "{OUTPUT_DIR}/{wildcards.sample}.markdup.bai" {output.bai}
        fi
        """



#Rule 3: Filter non-standard contigs
rule filter_standard_contigs:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.bam",
        bai = f"{OUTPUT_DIR}/{{sample}}.markdup.bam.bai"
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam"
    singularity: SAMTOOLS_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.filter_contigs.log"
    resources:
        # samtools-view + awk streams a 16-19GB BAM end-to-end (low-mem, but the
        # node still needs headroom). Auto-escalating on retry: 32 -> 64 -> 96 GB.
        mem_mb  = lambda wildcards, attempt: 32768 * attempt,
        runtime = lambda wildcards, attempt: min(600 * attempt, PARTITION_MAX_RUNTIME)
    shell:
        """
        exec &>> {log}
        samtools view -h {input.bam} \
        | awk 'BEGIN{{OFS="\t"}} /^@/ || $3 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X|Y|M)$/ {{print}}' \
        | samtools view -b -o {output.bam} -
        """




#Rule 4:Index "filtered bam files"
rule index_filtered_bam:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam"
    output:
        bai=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam.bai"
    singularity: SAMTOOLS_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.index_filtered_bam.log"
    resources:
        mem_mb=4096,
        runtime=120
    shell:
        """
        exec &>> {log}
        samtools index {input.bam} {output.bai}
        """



#Rule 5: SplitNCigar reads
rule split_n_cigar_reads:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam.bai",
        fasta=REF_FASTA,
        fai=REF_FAI,
        dict=REF_DICT
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.split.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.split.bai"
    singularity:
        GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.split_n_cigar_reads.log"
    resources:
        # Auto-escalating on retry: 40 -> 80 -> 120 GB and 60 -> 120 -> 180 min.
        mem_mb  = lambda wildcards, attempt: 40960 * attempt,
        runtime = lambda wildcards, attempt: min(3600 * attempt, PARTITION_MAX_RUNTIME)
    shell:
        """
        exec &>> {log}
        echo "===== $(date '+%F %T') | split_n_cigar_reads sample={wildcards.sample} | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" SplitNCigarReads \
            -R {config[reference][fasta]} \
            -I {input.bam} \
            -O {output.bam}
        """



#Rule 6: Calculate BQSR (Base Quality Score Recalibration) table
rule base_recalibrator:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.split.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.split.bai"
    output:
        table=f"{OUTPUT_DIR}/{{sample}}.recal_data.table"
    singularity:
        GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.base_recalibrator.log"
    resources:
        # Auto-escalating on retry: 40 -> 80 -> 120 GB and 60 -> 120 -> 180 min.
        mem_mb  = lambda wildcards, attempt: 40960 * attempt,
        runtime = lambda wildcards, attempt: min(3600 * attempt, PARTITION_MAX_RUNTIME)
    params:
        reference=config["reference"]["fasta"],
        dbsnp=config["reference"]["dbsnp"],
        mills=config["reference"]["mills_indels"]
    shell:
        """
        exec &>> {log}
        echo "===== $(date '+%F %T') | base_recalibrator sample={wildcards.sample} | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" BaseRecalibrator \
            -R {params.reference} \
            -I {input.bam} \
            --known-sites {params.dbsnp} \
            --known-sites {params.mills} \
            -O {output.table}
        """



# Rule 7: Apply  BQSR 
rule apply_BQSR:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.split.bam",
        table=f"{OUTPUT_DIR}/{{sample}}.recal_data.table"
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.BQSR.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.BQSR.bai"
    singularity: GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.gather_bqsr_reports.log"
    resources:
        # Auto-escalating on retry: 40 -> 80 -> 120 GB and 40 -> 80 -> 120 min.
        mem_mb  = lambda wildcards, attempt: 40960 * attempt,
        runtime = lambda wildcards, attempt: min(2400 * attempt, PARTITION_MAX_RUNTIME)
    params:
        reference=config["reference"]["fasta"]
    shell:
        """
        exec &>> {log}
        echo "===== $(date '+%F %T') | apply_BQSR sample={wildcards.sample} | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" ApplyBQSR \
            -R {params.reference} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output.bam}
        """



#Rule 8: Variant Calling 
rule haplotype_caller:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.BQSR.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.BQSR.bai",
        fasta=REF_FASTA,
        fai=REF_FAI,
        dict=REF_DICT
    output:
        vcf=f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz.tbi"
    singularity: GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.variant_calling.log"
    resources:
        # Auto-escalating on retry: 40 -> 80 -> 120 GB and 120 -> 240 -> 360 min.
        # HaplotypeCaller on full-genome RNA can take 1-3 hours per sample.
        mem_mb  = lambda wildcards, attempt: 40960 * attempt,
        runtime = lambda wildcards, attempt: min(7200 * attempt, PARTITION_MAX_RUNTIME)
    params:
        reference=config["reference"]["fasta"],
        intervals=" ".join([f"-L {c}" for c in config["chromosomes"]])
    shell:
        """
        exec &>> {log}
        echo "===== $(date '+%F %T') | haplotype_caller sample={wildcards.sample} | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" HaplotypeCaller \
            -ERC GVCF \
            -R {params.reference} \
            -I {input.bam} \
            {params.intervals} \
	    --dont-use-soft-clipped-bases true \
	    --standard-min-confidence-threshold-for-calling 20.0 \
            -O {output.vcf}
        """



#Rule 10: generate sample map for next step (create genomicsDB)
rule generate_sample_map:
    input: vcfs=expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz",sample=SAMPLES)
    output:
        map_file=f"{OUTPUT_DIR}/sample_map_file"
    run:
        with open(output.map_file, 'w') as f:
            for vcf_path in input.vcfs:
                sample_name = os.path.basename(vcf_path).replace('.g.vcf.gz', '')
                f.write(f"{sample_name}\t{vcf_path}\n")




# ==============================================================================
# OPTIMIZATION: FULL SCATTER-GATHER (DB + GENOTYPING)
# ==============================================================================
# Rule 11: GenomicsDB Import - SCATTERED
# Creates one DB folder per chromosome
rule GenomicsDB_scatter:
    input:
        sample_map = f"{OUTPUT_DIR}/sample_map_file",
        # We need the GVCFs to be ready, but the tool reads them from the sample_map
        gvcfs = expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz", sample=SAMPLES),
        tbis = expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES)
    output:
        # [CHANGED] Use directory() to mark the folder as output
        db_dir = directory(f"{OUTPUT_DIR}/genomicsdb_per_chrom/genomicsdb_{{chrom}}"),
        flag = touch(f"{OUTPUT_DIR}/genomicsdb_per_chrom/{{chrom}}.done")
    singularity: GATK_CONTAINER
    log:
        f"{LOG_DIR}/GenomicsDBImport_{{chrom}}.log"
    params:
        # [CHANGED] Removed tmp_dir
        interval = "{chrom}" 
    threads: 4
    resources:
        # GenomicsDBImport is memory-hungry. Auto-escalating: 24 -> 48 -> 72 GB.
        # Java heap = mem_mb-2GB (see -Xmx below).
        mem_mb  = lambda wildcards, attempt: 24576 * attempt,
        runtime = lambda wildcards, attempt: min(300 * attempt, PARTITION_MAX_RUNTIME)
    shell:
        """
        exec &>> {log}

        # [CHANGED] Removed mkdir tmp and -Djava.io.tmpdir
        # Ensure clean start by removing existing DB dir if retrying
        rm -rf {output.db_dir}

        echo "===== $(date '+%F %T') | GenomicsDB_scatter chrom={wildcards.chrom} | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" GenomicsDBImport \
            --genomicsdb-workspace-path {output.db_dir} \
            --sample-name-map {input.sample_map} \
            --batch-size 5 \
            -L {params.interval} \
            --reader-threads {threads} \
            --merge-input-intervals \
            --consolidate
        """

# Rule 12a: GenotypeGVCFs - SCATTERED
# Joins genotyping per chromosome
rule join_genotyping_scatter:
    input:
        # Triggered by the flag file from the previous step
        db_flag = f"{OUTPUT_DIR}/genomicsdb_per_chrom/{{chrom}}.done",
        db_dir = f"{OUTPUT_DIR}/genomicsdb_per_chrom/genomicsdb_{{chrom}}"
    output:
        vcf = temp(f"{OUTPUT_DIR}/split_vcfs/Joint_{{chrom}}.vcf.gz"),
        tbi = temp(f"{OUTPUT_DIR}/split_vcfs/Joint_{{chrom}}.vcf.gz.tbi")
    log:
        f"{LOG_DIR}/genotypeGVCFs_{{chrom}}.log"
    singularity: GATK_CONTAINER
    params:
        reference = config["reference"]["fasta"],
        # [CHANGED] Removed tmp_dir
        chrom = "{chrom}"
    resources:
        # 24GB base — GenotypeGVCFs needs >16GB heap for hg38 joint calling.
        # Auto-escalating on retry: 24 -> 48 -> 72 GB (heap = mem_mb-2GB).
        mem_mb  = lambda wildcards, attempt: 24576 * attempt,
        runtime = lambda wildcards, attempt: min(600 * attempt, PARTITION_MAX_RUNTIME)
    shell:
        """
        exec &>> {log}
        
        # [CHANGED] Removed mkdir tmp and -Djava.io.tmpdir
        
        echo "===== $(date '+%F %T') | join_genotyping_scatter chrom={wildcards.chrom} | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" GenotypeGVCFs \
            -R {params.reference} \
            -V gendb://{input.db_dir} \
            -L {params.chrom} \
            -O {output.vcf}
        """

# Rule 12b: Gather VCFs
# Merges the scattered chromosome VCFs into one
rule merge_joint_vcfs:
    input:
        vcfs = expand(f"{OUTPUT_DIR}/split_vcfs/Joint_{{chrom}}.vcf.gz", chrom=CHROMOSOMES),
        tbis = expand(f"{OUTPUT_DIR}/split_vcfs/Joint_{{chrom}}.vcf.gz.tbi", chrom=CHROMOSOMES)
    output:
        vcf = f"{OUTPUT_DIR}/Joint_all.vcf.gz",
        tbi = f"{OUTPUT_DIR}/Joint_all.vcf.gz.tbi"
    singularity: GATK_CONTAINER
    log:
        f"{LOG_DIR}/merge_vcfs.log"
    resources:
        mem_mb = 16384,
        runtime = 120
    params:
        input_args = lambda wildcards, input: " -I ".join(input.vcfs)
    shell:
        """
        exec &>> {log}
        echo "===== $(date '+%F %T') | merge_joint_vcfs (all chroms) | host=$(hostname) job=${{SLURM_JOB_ID:-local}} mem={resources.mem_mb}MB rt={resources.runtime}min ====="
        gatk --java-options "-Xmx$(( {resources.mem_mb} - 2048 ))M" GatherVcfs \
            -I {params.input_args} \
            -O {output.vcf}
        
        gatk IndexFeatureFile -I {output.vcf}

        """
