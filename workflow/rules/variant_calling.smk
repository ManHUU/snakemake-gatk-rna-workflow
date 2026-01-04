# GATK RNA-seq Variant Calling Pipeline - Snakemake Workflow


# Rule 0a: Prepare FAI (Uses Samtools)
rule prepare_fasta_index:
    input: REF_FASTA
    output: REF_FAI
    container: SAMTOOLS_CONTAINER
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
    container: GATK_CONTAINER
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
    container: STAR_CONTAINER
    log: f"{LOG_DIR}/build_star_index.log"
    threads: 16
    resources:
        mem_mb=64000,       # STAR indexing requires significant RAM (~30GB+)
        runtime=300
    shell:
        """
        exec &> {log}
        
        # Create the output directory if it doesn't exist
        mkdir -p {output}

        # Run STAR in genomeGenerate mode
        # --sjdbOverhang 100 is standard for read lengths of 100bp+
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 100
        """


# Rule 1: Align FASTQ to Reference and Sort
rule star_align:
    input:
        # Assuming single-end based on your file list; 
        # if paired-end, you'd use r1 and r2 here.
        fastq=lambda wildcards: os.path.join(FASTQ_DIR, f"{wildcards.sample}.fastq"),
        star_index = config["reference"]["star_index_dir"]
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.sorted.bam",
    container: STAR_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.star_align.log"
    params:
        sample_id = "{sample}",
        out_prefix = f"{OUTPUT_DIR}/{{sample}}."
    threads: 24
    resources:
        mem_mb=32768,        # 32GB
        runtime=600
    shell:
        """
        exec &> {log}
        # Generates coordinate-sorted BAM and adds Read Groups required for GATK
        STAR --runThreadN {threads} \
             --genomeDir {input.star_index} \
             --readFilesIn {input.fastq} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattrRGline ID:{params.sample_id} SM:{params.sample_id} LB:lib1 PL:ILLUMINA PU:unit1 \
             --outFileNamePrefix {params.out_prefix}

        # Rename STAR's default output name to match Snakemake's output
        mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam}
        """


# Rule 1b: Index the Aligned BAM (Uses Samtools Container)
# This rule picks up where Rule 1 left off.
rule index_sorted_bam:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.sorted.bam"
    output:
        bai=f"{OUTPUT_DIR}/{{sample}}.sorted.bam.bai"
    container: SAMTOOLS_CONTAINER
    log: f"{LOG_DIR}/{{sample}}.index_star_bam.log"
    threads: 4
    resources:
        mem_mb=4096,
        runtime=60
    shell:
        """
        exec &> {log}
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
    container:
        GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.mark_duplicates.log"
    threads: 4
    resources:
        mem_mb=32768,
        runtime=300
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}M" MarkDuplicates \
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
    container: SAMTOOLS_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.filter_contigs.log"
    resources:
        mem_mb=32768,
        runtime=300
    shell:
        """
        exec &> {log}
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
    container: SAMTOOLS_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.index_filtered_bam.log"
    shell:
        """
        exec &> {log}
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
    container:
        GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.split_n_cigar_reads.log"
    resources:
        mem_mb=32768,
        runtime=600
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}M" SplitNCigarReads \
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
    container:
        GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.base_recalibrator.log"
    resources:
        mem_mb=32768,
        runtime=300
    params:
        reference=config["reference"]["fasta"],
        dbsnp=config["reference"]["dbsnp"],
        mills=config["reference"]["mills_indels"]
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}M" BaseRecalibrator \
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
    container: GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.gather_bqsr_reports.log"
    resources:
        mem_mb=32768,
        runtime=300
    params:
        reference=config["reference"]["fasta"]
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}M" ApplyBQSR \
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
    container: GATK_CONTAINER
    log:
        f"{LOG_DIR}/{{sample}}.variant_calling.log"
    resources:
        mem_mb=32768,
        runtime=600
    params:
        reference=config["reference"]["fasta"],
        intervals=" ".join([f"-L {c}" for c in config["chromosomes"]])
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}M" HaplotypeCaller \
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



#Rule 11: generate genomicsDB
rule GenomicsDB:
    input:
        sample_map=f"{OUTPUT_DIR}/sample_map_file",
        gvcfs=expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz", sample=SAMPLES),
        tbis=expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES)
    output:
        touch(f"{OUTPUT_DIR}/GenomicsDB.done")
    container: GATK_CONTAINER
    params:
        workspace=f"{OUTPUT_DIR}/genomicsdb",
        chromosomes=" ".join([f"-L {chrom}" for chrom in CHROMOSOMES])
    log:
        f"{LOG_DIR}/genomicsDB.log"
    threads: 4
    resources:
        mem_mb=32768,
        runtime=600
    shell:
        """
        exec &> {log}
        rm -rf {params.workspace}
        gatk --java-options "-Xmx{resources.mem_mb}M" GenomicsDBImport \
            --genomicsdb-workspace-path {params.workspace} \
            --sample-name-map {input.sample_map} \
	    --batch-size 3 \
	    {params.chromosomes} \
	    --reader-threads {threads} \
            && touch {output}
        """




#Rule 12: perform joint genotyping 
rule join_genotyping:
    input:
        f"{OUTPUT_DIR}/GenomicsDB.done"
    output:
        f"{OUTPUT_DIR}/Joint_all.vcf.gz"
    log:
        f"{LOG_DIR}/genotypeGVCFs.log"
    container: GATK_CONTAINER
    params:
        reference=config["reference"]["fasta"],
	database=f"{OUTPUT_DIR}/genomicsdb"
    threads: config["resources"]["join_genotyping"]["threads"]
    resources:
        mem_mb=config["resources"]["join_genotyping"]["memory"] * 1024
    shell:
        """
        exec &> {log}
        
        echo "Starting joint genotyping..."
        echo "Date: $(date)"
        gatk --java-options "-Xmx{resources.mem_mb}M" GenotypeGVCFs \
           -R {params.reference} \
           -V gendb://{params.database} \
           -O {output}
        """



#Rule 13: Annotate variant filter 
## This step won't filter out any reads, the number of reads of output will be the same as input
rule annotate_variant_filter:
    input:
        f"{OUTPUT_DIR}/Joint_all.vcf.gz"
    output:
        f"{OUTPUT_DIR}/Filtered.vcf.gz"
    container:
        GATK_CONTAINER
    log:
        f"{LOG_DIR}/variant_filtration.log"
    params:
        reference=config["reference"]["fasta"],
        filter_expression=config["variant_filters"]["expression"],
        filter_name=config["variant_filters"]["name"]
    resources:
        mem_mb= 65536
    shell:
        """
        exec &> {log}
        echo "variant_filtration"
        gatk --java-options "-Xmx{resources.mem_mb}M" VariantFiltration \
            -R {params.reference} \
            -V {input} \
            -O {output} \
            --filter-name "{params.filter_name}" \
            --filter-expression "{params.filter_expression}"
       """
