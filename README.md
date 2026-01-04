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
├── config/                 # Configuration files
│   └── config.yaml         # MAIN CONFIG: Paths, memory, threads
├── resources/              # Place references here (not tracked by git)
├── workflow/               # Snakemake logic
│   ├── rules/              # .smk rule definitions
│   ├── profiles/           # Slurm configuration
│   └── scripts/            # Submission scripts
├── Snakefile               # Main entry point
└── README.md
