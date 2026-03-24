# Snakemake GATK RNA-seq Variant Calling Workflow

A reproducible, containerized Snakemake pipeline for RNA-seq variant calling
following GATK Best Practices. Supports both local workstations and HPC/Slurm
environments. All tools run inside Singularity/Apptainer containers — no manual
tool installation required.

---

## Features

- Paired-end and single-end FASTQ support (set one flag in config)
- Adapter trimming with fastp
- Raw read quality control with FastQC + MultiQC
- STAR 2-pass alignment with read group tagging
- GATK Best Practices: MarkDuplicates → SplitNCigarReads → BQSR → HaplotypeCaller
- Scalable joint genotyping: GenomicsDBImport + GenotypeGVCFs scattered per chromosome
- Type-specific hard filtering (SNP and INDEL, RNA-seq tuned thresholds)
- RNA editing site removal (REDIportal hg38)
- Dense variant cluster filtering
- Functional annotation with SnpEff (hg38)
- CDS-only filter: retains coding variants (missense, frameshift, stop, etc.)
- rsID annotation from dbSNP138
- Publication-quality figures (6 plots, PNG + PDF)
- Runs locally or on Slurm HPC with no code changes

---

## Pipeline Workflow

The pipeline is split into three rule modules, executed in order.

```mermaid
flowchart TD
    FASTQ[/"FASTQ files\n{sample}_1.fastq / {sample}_2.fastq"/]

    subgraph QC["Module 1 — qc.smk"]
        direction TB
        QC1["FastQC\nper-sample quality report"]
        QC2["MultiQC\naggregated report\n(FastQC + STAR logs + MarkDup metrics)"]
        QC1 --> QC2
    end

    subgraph VC["Module 2 — variant_calling_jointCall_per_Chr.smk"]
        direction TB
        FASTP["fastp\nadapter trimming, poly-G removal\nmin length 36 bp"]
        STAR["STAR align\ncoordinate-sorted BAM\n+ read group tags"]
        MKDUP["GATK MarkDuplicates\nflag PCR duplicates"]
        FILT["samtools view\nkeep standard chromosomes only"]
        SPLIT["GATK SplitNCigarReads\nsplit reads spanning splice junctions"]
        BQSR["GATK BaseRecalibrator\n+ ApplyBQSR\nbase quality score recalibration"]
        HC["GATK HaplotypeCaller\nper-sample GVCF\n(ERC GVCF mode)"]
        MAP["create_sample_map\nGVCF path list for joint calling"]
        DB["GenomicsDBImport\nscattered per chromosome\n(25 parallel jobs)"]
        GT["GenotypeGVCFs\nscattered per chromosome"]
        MERGE["MergeVcfs\nJoint_all.vcf.gz"]

        FASTP --> STAR --> MKDUP --> FILT --> SPLIT --> BQSR --> HC --> MAP --> DB --> GT --> MERGE
    end

    subgraph FV["Module 3 — filter_and_visualize.smk"]
        direction TB
        FILTERING["Filtering Strategy\n(see diagram below)"]
        STATS["variant_statistics\nstats/*.tsv"]
        VIZ["visualize_variants\nresults/plots/"]
        FILTERING --> STATS --> VIZ
    end

    FASTQ --> QC1
    FASTQ --> FASTP
    STAR --> QC2
    MKDUP --> QC2
    MERGE --> FILTERING
```

> **Parallelism**: HaplotypeCaller runs per sample in parallel. GenomicsDBImport
> and GenotypeGVCFs each run as 25 independent jobs (one per chromosome),
> then MergeVcfs collects the results.

---

## Filtering Strategy

```mermaid
flowchart TD
    IN[/"Joint_all.vcf.gz\nraw joint-called variants"/]

    subgraph HARDFILTER["Hard Filtering"]
        direction LR
        SNP["filter_snps\nGATK SelectVariants + VariantFiltration\n─────────────────────\nQD < 2.0\nFS > 30.0\nMQ < 40.0\nMQRankSum < -12.5\nReadPosRankSum < -8.0\nDP < 10"]
        INDEL["filter_indels\nGATK SelectVariants + VariantFiltration\n─────────────────────\nQD < 2.0\nFS > 200.0\nReadPosRankSum < -20.0\nDP < 10"]
    end

    PASS["select_pass_variants\nMergeVcfs + SelectVariants\nkeep PASS-only variants\n→ Final_PASS.vcf.gz"]

    EDIT["remove_rna_editing\nbcftools view -T\nexclude REDIportal hg38 known A-to-I sites\n→ Final_PASS_noEdit.vcf.gz"]

    CLUSTER["filter_clusters\nGATK VariantFiltration\nexclude dense clusters\nwindow=35 bp, min=3 variants\n→ Final_clean.vcf.gz"]

    subgraph ANNOTATION["Annotation"]
        direction TB
        ANN["annotate_variants\nSnpEff hg38\nadd ANN fields\n→ Final_annotated.vcf.gz"]
        ANNREP[["Final_annotated_healthCheck.html"]]
        CDS["filter_cds\nSnpSift filter\nkeep coding effects only:\nmissense · synonymous · stop_gained\nstop_lost · start_lost · frameshift\ninframe_ins · inframe_del\n→ Final_CDS.vcf.gz"]
        RSID["annotate_rsid\nbcftools annotate\nadd dbSNP138 rsIDs\n→ Final_CDS_rsID.vcf.gz"]
        CDSREP[["Final_CDS_healthCheck.html"]]
        ANN --> ANNREP
        ANN --> CDS --> RSID --> CDSREP
    end

    OUT[/"Final_CDS_rsID.vcf.gz\nFinal output — coding variants with rsIDs"/]

    IN --> SNP & INDEL --> PASS --> EDIT --> CLUSTER --> ANN
    RSID --> OUT
```

---

## Repository Structure

```
.
├── Snakefile                          # Main entry point
├── config/
│   └── config.yaml                    # All user-facing settings
├── workflow/
│   ├── rules/
│   │   ├── qc.smk                     # Module 1: FastQC + MultiQC
│   │   ├── variant_calling_jointCall_per_Chr.smk   # Module 2: fastp → STAR → GATK
│   │   └── filter_and_visualize.smk   # Module 3: filtering → annotation → plots
│   ├── profiles/
│   │   └── local/config.yaml          # Snakemake profile for local runs
│   └── scripts/
│       ├── run_local.sh               # Run pipeline locally
│       ├── download_refs.slurm        # Download all reference files
│       ├── download_GSE256519.sh      # Download example test data (SRA)
│       └── visualize_variants.py      # Variant visualization script
├── resources/                         # Reference files (not tracked by git)
└── data/                              # Input FASTQ files (not tracked by git)
```

---

## Requirements

- **Conda** (Miniconda or Mamba) — only needed for Snakemake itself
- **Apptainer / Singularity** — all tools run in containers
- No GATK, STAR, or samtools installation needed

Install Snakemake:
```bash
conda create -n snakemake_env -c bioconda -c conda-forge snakemake
conda activate snakemake_env
```

Install SnpEff/SnpSift (required for annotation, runs outside containers):
```bash
conda create -n vcf_annotation -c bioconda snpeff snpsift
```

---

## Reference Files

Download all required references with the provided script:
```bash
bash workflow/scripts/download_refs.slurm
# or on HPC:
sbatch workflow/scripts/download_refs.slurm
```

This downloads:

| File | Used by |
|---|---|
| GRCh38 primary assembly FASTA | STAR index, GATK |
| Gencode v44 GTF | STAR index |
| Mills & 1000G gold standard indels (hg38) | BQSR |
| dbSNP138 hg38 VCF | BQSR, rsID annotation |
| REDIportal v2.0 hg38 BED | RNA editing filter |
| SnpEff hg38 database | Functional annotation |

After downloading, update absolute paths in `config/config.yaml` if you placed
files outside the default `resources/` directory.

---

## Configuration

All settings are in `config/config.yaml`:

```yaml
fastq_dir: "data/GSE256519"      # Directory containing input FASTQ files
output_dir: "results"
log_dir: "logs"

sequencing_type: "paired"        # "paired" or "single"
                                 # paired  → {sample}_1.fastq + {sample}_2.fastq
                                 # single  → {sample}.fastq
```

Sample names are discovered automatically by globbing `fastq_dir` at runtime —
no sample list needed.

---

## Usage Example

The following walks through a complete run from a fresh clone to final results,
using the GSE256519 test dataset (human heart RNA-seq, paired-end 2×150 bp).

### Step 1 — Clone the repository

```bash
git clone https://github.com/ManHUU/snakemake-gatk-rna-workflow.git
cd snakemake-gatk-rna-workflow
```

### Step 2 — Install Snakemake

```bash
conda create -n snakemake_env -c bioconda -c conda-forge snakemake
conda activate snakemake_env
```

### Step 3 — Download reference files

```bash
bash workflow/scripts/download_refs.slurm
```

This places all references under `resources/` and the paths in `config/config.yaml`
will match without any edits.

### Step 4 — Download input FASTQ data

```bash
bash workflow/scripts/download_GSE256519.sh
```

> **Note:** This downloads the GSE256519 test dataset. To use your own data,
> place your FASTQ files directly in `data/` (or any directory of your choice)
> and update `fastq_dir` in `config/config.yaml` to point to it.

After this you should have:
```
data/GSE256519/
├── SRR28074879_1.fastq
├── SRR28074879_2.fastq
├── SRR28074880_1.fastq
└── SRR28074880_2.fastq
```

Samples are discovered automatically — no need to edit the config.

### Step 5 — Dry run (verify the DAG)

```bash
bash workflow/scripts/run_local.sh --dry-run
```

You should see all rules listed for both samples across 25 chromosomes with no errors.

### Step 6 — Run the pipeline

```bash
# Run in a screen session so it survives disconnection
screen -S gatk_run
bash workflow/scripts/run_local.sh
# Detach with Ctrl+A then D
```

### Step 7 — Check outputs

```bash
# Final variant file
ls -lh results/Final_CDS_rsID.vcf.gz

# QC report
xdg-open results/qc/multiqc_report.html

# SnpEff health check
xdg-open results/Final_CDS_healthCheck.html

# Figures
ls results/plots/
```

### Expected output summary

```
results/
├── qc/
│   └── multiqc_report.html          # FastQC + STAR + MarkDuplicates QC
├── Final_CDS_rsID.vcf.gz            # Final output: coding variants with rsIDs
├── Final_annotated_healthCheck.html # SnpEff report (all variants)
├── Final_CDS_healthCheck.html       # SnpEff report (coding variants only)
├── stats/
│   ├── quality_metrics.tsv
│   ├── filter_summary.tsv
│   ├── chrom_counts.tsv
│   ├── sample_counts.tsv
│   └── ts_tv.tsv
└── plots/
    ├── variant_type_summary.png / .pdf
    ├── filter_summary.png / .pdf
    ├── quality_distributions.png / .pdf
    ├── chrom_distribution.png / .pdf
    ├── sample_counts.png / .pdf
    └── ts_tv_comparison.png / .pdf
```

---

## Running the Pipeline

### Local workstation

```bash
bash workflow/scripts/run_local.sh            # full run
bash workflow/scripts/run_local.sh --dry-run  # preview steps only
```

The local profile uses 8 cores by default. Edit
`workflow/profiles/local/config.yaml` to change.

### HPC / Slurm

```bash
# Dry run first
snakemake --snakefile Snakefile --profile workflow/profiles/slurm --dry-run

# Full run
sbatch workflow/scripts/run_pipeline.sh
```

---

## Outputs

| Path | Description |
|---|---|
| `results/qc/multiqc_report.html` | Aggregated QC (FastQC + STAR + MarkDuplicates) |
| `results/Joint_all.vcf.gz` | Raw joint-called variants (pre-filter) |
| `results/Final_PASS.vcf.gz` | Hard-filtered PASS variants |
| `results/Final_PASS_noEdit.vcf.gz` | RNA editing sites removed |
| `results/Final_clean.vcf.gz` | Dense clusters removed |
| `results/Final_annotated.vcf.gz` | SnpEff functional annotation |
| `results/Final_CDS.vcf.gz` | Coding variants only |
| `results/Final_CDS_rsID.vcf.gz` | **Final output** — coding variants with rsIDs |
| `results/Final_annotated_healthCheck.html` | SnpEff summary (all variants) |
| `results/Final_CDS_healthCheck.html` | SnpEff summary (coding variants) |
| `results/stats/quality_metrics.tsv` | Per-variant QD, FS, MQ, DP (pre-filter) |
| `results/stats/ts_tv.tsv` | Ts/Tv ratio before and after filtering |
| `results/plots/` | 6 publication-quality figures (PNG + PDF) |

### Figures generated

| Figure | Description |
|---|---|
| `variant_type_summary` | SNP / INDEL counts, PASS vs filtered |
| `filter_summary` | Total PASS vs filtered (log scale) |
| `quality_distributions` | QD and DP histograms per variant type |
| `chrom_distribution` | Per-chromosome variant counts |
| `sample_counts` | Per-sample coding variant counts |
| `ts_tv_comparison` | Ts/Tv counts and ratio before vs after filtering |

---

## Containers

All containers are pulled automatically on first run and cached in
`resources/containers/`.

| Tool | Container | Version |
|---|---|---|
| STAR | quay.io/biocontainers/star | 2.7.11b |
| samtools | quay.io/biocontainers/samtools | 1.19 |
| GATK | broadinstitute/gatk | 4.6.1.0 |
| fastp | quay.io/biocontainers/fastp | 0.23.4 |
| FastQC | quay.io/biocontainers/fastqc | 0.12.1 |
| MultiQC | quay.io/biocontainers/multiqc | 1.21 |
| bcftools | quay.io/biocontainers/bcftools | 1.21 |
| matplotlib | quay.io/biocontainers/matplotlib | 3.5.1 |

---

## Test Data

The pipeline was developed using human heart RNA-seq data from GEO dataset
GSE256519 (paired-end 2×150 bp). To download the test samples:

```bash
bash workflow/scripts/download_GSE256519.sh
```

This downloads SRR28074879 and SRR28074880 into `data/GSE256519/`.

---

## License

See LICENSE file.
