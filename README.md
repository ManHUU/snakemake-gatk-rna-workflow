# Snakemake Workflow: RNA-Seq Variant Calling (GATK)

This repository contains a reproducible Snakemake workflow for RNA-Seq variant calling using the GATK Best Practices pipeline. It is designed to run on HPC clusters using **Slurm** for job scheduling and **Apptainer (Singularity)** for containerized environment management.

## 📌 Features
* **Reproducible:** Uses Docker/Singularity containers for all major tools (STAR, GATK, Samtools).
* **Scalable:** Optimized for Slurm HPC environments with automatic resource scaling.
* **Conflict-Free:** Separate containers for separate rules to avoid dependency hell.
* **GATK Best Practices:** Implements SplitNCigarReads, BaseRecalibrator, and HaplotypeCaller.

## 🛠️ Prerequisites
You do not need to install GATK or STAR. You only need:
1.  **HPC Cluster** with Slurm scheduler.
2.  **Apptainer (or Singularity)** installed on the cluster.
3.  **Conda** (Miniconda/Mamba) to run Snakemake.

## 📂 Repository Structure
```text
.
├── data/                   # Test data
├── config/                 # Configuration files
│   └── config.yaml         # MAIN CONFIG: Paths, memory, threads
├── resources/              # Place references here (not tracked by git)
├── workflow/               # Snakemake logic
│   ├── rules/              # .smk rule definitions
│   ├── profiles/           # Slurm configuration
│   └── scripts/            # Submission scripts
├── Snakefile               # Main entry point
└── README.md
```

## 🚀 Quick Start
### 1. Clone the repository

```Bash
git clone [https://github.com/ManHUU/snakemake-gatk-rna-workflow.git](https://github.com/ManHUU/snakemake-gatk-rna-workflow.git)
cd snakemake-gatk-rna-workflow
```

### 2. Prepare Reference Files
Place your reference genome (.fa), GTF, and known sites (.vcf) in the resources/ folder.
Update config/config.yaml if your filenames differ from the defaults.

### 3. Install Snakemake
Create a simple environment for the workflow manager:
```Bash
conda create -n snakemake_env -c bioconda -c conda-forge snakemake=7.32.4
conda activate snakemake_env
```

### 4. Configure & Run
Edit config/config.yaml to point to your FASTQ directory. Then, submit the pipeline to Slurm:
```Bash
sbatch workflow/scripts/run_pipeline.sh
```


### 📦 Containers Used
The workflow automatically pulls these containers (no manual installation required):
```text
STAR: quay.io/biocontainers/star:2.7.11b--h43eeafb_3

Samtools: quay.io/biocontainers/samtools:1.19--h50ea8bc_0

GATK: broadinstitute/gatk:4.6.1.0
```



