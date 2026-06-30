# Snakemake GATK RNA-seq Variant Calling Workflow

A reproducible Snakemake pipeline for RNA-seq variant calling following GATK
Best Practices. Supports both local workstations and HPC/Slurm environments.
Every scientific tool (STAR, GATK, fastp, FastQC, MultiQC, bcftools, SnpEff,
SnpSift, SRA toolkit) runs inside a pinned Singularity/Apptainer container —
the only software you install on the host is Snakemake itself.

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
        DB["GenomicsDBImport\nscattered per chromosome\n(24 parallel jobs)"]
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
> and GenotypeGVCFs each run as 24 independent jobs (one per chromosome — chr1–22, chrX, chrY),
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
├── run.sh                             # Unified runner — local OR SLURM (auto-detect)
├── config/
│   ├── config.yaml                    # All user-facing settings
│   └── units.tsv                      # Sample sheet (sample / unit / fq1 / fq2)
├── workflow/
│   ├── envs/
│   │   └── snakemake.yaml             # Snakemake 9 + slurm executor plugin (only host-side env)
│   ├── profiles/
│   │   ├── local/config.yaml          # Snakemake profile for local runs
│   │   └── slurm/config.yaml          # Snakemake profile for SLURM HPC runs
│   ├── schemas/
│   │   ├── config.schema.yaml         # JSON schema validating config.yaml
│   │   └── units.schema.yaml          # JSON schema validating the sample sheet
│   ├── rules/
│   │   ├── qc.smk                     # Module 1: FastQC + MultiQC
│   │   ├── variant_calling_jointCall_per_Chr.smk   # Module 2: fastp → STAR → GATK
│   │   └── filter_and_visualize.smk   # Module 3: filtering → annotation → plots
│   └── scripts/
│       ├── download_refs.sh           # Download all reference databases (checksum + Zenodo fallback)
│       ├── download_refs.slurm        # Thin SLURM wrapper around download_refs.sh
│       ├── make_units.sh              # Generate config/units.tsv from a FASTQ dir (run once)
│       ├── download_GSE256519.sh      # Download example test data (SRA)
│       ├── download_GSE256519.slurm   # SLURM wrapper for the test-data download
│       └── visualize_variants.py      # Variant visualization script
├── resources/                         # Reference files (not tracked by git)
└── data/                              # Input FASTQ files (not tracked by git)
```

> **Files an HPC reader must touch:** none. Set the two env vars
> `SLURM_ACCOUNT` and `SLURM_PARTITION` for your cluster (see Usage Example).
> Local readers edit nothing.

---

## Requirements

- **Snakemake** — the only host-side software you must install.
  Install however you prefer:
  ```bash
  conda env create -f workflow/envs/snakemake.yaml      # or: mamba / micromamba
  conda activate snakemake_env
  ```
  Other valid options: `pip install snakemake snakemake-executor-plugin-slurm`
  in a venv, or `module load snakemake` on HPC systems with environment modules.
- **Apptainer or Singularity** — required. Every scientific tool runs from a
  pinned biocontainer pulled on first use.
- No manual installation of GATK, STAR, samtools, SnpEff, SnpSift, SRA toolkit,
  bcftools, etc. — they all run from containers.

Run `bash workflow/scripts/check_prerequisites.sh` to verify your host before
starting.

---

## Reference Files

Download all required references with the provided script. Output paths
already match `config/config.yaml`, so no edits are required afterwards.
Step 5/5 (SnpEff hg38 database) runs `snpEff` from the same biocontainer
the pipeline uses — apptainer/singularity must be on PATH.

```bash
# Local (or HPC interactive node):
bash workflow/scripts/download_refs.sh

# HPC via SLURM — uses the same env vars as run.sh (set them once per shell):
export SLURM_ACCOUNT=<your_slurm_account>
export SLURM_PARTITION=<your_slurm_partition>
bash workflow/scripts/download_refs.slurm   # self-submits via sbatch
```

This downloads:

| File | Source | Used by |
|---|---|---|
| GRCh38 primary assembly FASTA | EBI / GENCODE | STAR index, GATK |
| GENCODE v44 GTF | EBI / GENCODE | STAR index |
| Mills & 1000G gold standard indels (hg38) | Broad Institute | BQSR |
| dbSNP138 hg38 VCF | Broad Institute | BQSR, rsID annotation |
| REDIportal v2.0 hg38 BED | rediportal.cloud.ba.infn.it | RNA editing filter |
| SnpEff hg38 database | snpeff.blob.core.windows.net | Functional annotation |

### Integrity checking and fallback mirror

Every download is verified against a pinned **sha256** checksum immediately after
fetching (see the `SHA_*` values and the `fetch_with_fallback` / `verify_sha256`
helpers at the top of `download_refs.sh`). A corrupted, truncated, or silently
changed file fails loudly instead of poisoning the run, and re-running the script
skips anything already present and verified.

The two link-rot-prone files — the REDIportal BED and the SnpEff hg38 database —
are also archived on Zenodo as a permanent fallback. If the upstream host is
unreachable, the script automatically falls back to the Zenodo copy and verifies
the same checksum:

> **Zenodo archive (DOI):** https://zenodo.org/records/21030946

---

## Disk space and runtime expectations

This is a heavyweight pipeline. Intermediate BAMs accumulate quickly and the
references alone are substantial. Plan your storage before launching —
especially on HPC, where the project filesystem is quota-bound.

**One-time downloads (`download_refs.sh` + test data):**

| Item | Size |
|---|---|
| STAR hg38 index | ~30 GB |
| GRCh38 FASTA + GTF + dbSNP + Mills + REDIportal + SnpEff db | ~10 GB |
| GSE256519 test FASTQs (2 samples, paired-end) | ~50 GB |

**Per-sample intermediate files (peak, before cleanup):**

| Step | Output written | Approx size | Approx wall time @ 24 cores |
|---|---|---|---|
| `fastp_trim` | trimmed FASTQ | ~ input size | 5–15 min |
| `star_align` | sorted BAM + scratch `_STARtmp/` | ~15–20 GB BAM + ~20–40 GB scratch | 15–30 min |
| `mark_duplicates` | markdup BAM | ~15–20 GB | 30–90 min |
| `filter_standard_contigs` | filtered BAM | ~15–20 GB | 15–30 min |
| `split_n_cigar_reads` | split BAM | ~15–20 GB | 1–3 hr |
| `base_recalibrator` + `apply_BQSR` | recal table + BQSR BAM | ~15–20 GB | 1–3 hr |
| `haplotype_caller` | per-sample GVCF | ~1–3 GB | 1–3 hr |
| `GenomicsDB_scatter` (× 24 chrom) | DB workspace | ~5–10 GB total | 5–10 min each |
| `join_genotyping_scatter` + `merge_joint_vcfs` | joint VCF | < 1 GB | 10–30 min |

Per sample this leaves ~80–120 GB of intermediate BAMs sitting in `results/`
until you clean them up. The final VCFs and QC reports total <2 GB. For a
two-sample test run, budget **~300 GB of free project space** including
references; for a typical cohort scale linearly per sample.

### HPC tip — STAR scratch must go to node-local SSD

STAR's BAM sort streams ~20–40 GB of intermediate bin files to a scratch
directory. The rule routes this to Snakemake's `tmpdir` resource, which under
SLURM resolves to node-local `/scratch-local/<user>.<jobid>` — **never** to
your project filesystem. If your `_STARtmp/` ever lands on the project FS
(e.g. because you bypassed `run.sh` and lost the `--bind /scratch-local`
arg), STAR fails mid-sort with:

```
EXITING because of FATAL ERROR: number of bytes expected from the BAM bin
does not agree with the actual size on disk
```

which usually means the project quota or per-file cap got hit. Fix: launch
via `run.sh` (which binds `/scratch-local`), and check `myquota` if you see
this error.

### Cleaning up after a run

Once `Final_CDS_rsID.vcf.gz` and the QC reports are in place, the intermediate
BAMs are safe to delete:

```bash
rm results/*.sorted.bam results/*.sorted.bam.bai
rm results/*.markdup.bam results/*.markdup.bam.bai results/*.markdup_metrics.txt
rm results/*.markdup.filtered.bam results/*.markdup.filtered.bam.bai
rm results/*.split.bam results/*.split.bai
rm results/*.BQSR.bam results/*.BQSR.bai
rm results/*.recal_data.table
```

Keep the GVCFs (`*.g.vcf.gz`) and the genomicsdb workspace if you want to
add samples later without re-running per-sample steps.

---

## Configuration

All settings are in `config/config.yaml`:

```yaml
fastq_dir: "data/GSE256519"      # only used by make_units.sh, not at runtime
output_dir: "results"
log_dir: "logs"

units: "config/units.tsv"        # sample sheet (sample / unit / fq1 / fq2)

sequencing_type: "paired"        # "paired" or "single"
                                 # paired  → {sample}_1.fastq + {sample}_2.fastq
                                 # single  → {sample}.fastq
```

`config/config.yaml` and the sample sheet are validated against JSON schemas in
`workflow/schemas/` at DAG-build time, so a mistyped field or missing path fails
immediately with a clear message instead of mid-run.

### Sample sheet (`config/units.tsv`)

The set of samples to analyse is specified explicitly in a small tab-separated sample
sheet, `config/units.tsv`, with one row per sample (columns: `sample`, `unit`, `fq1`,
`fq2`, where `fq1`/`fq2` give the paths to the FASTQ files). The pipeline reads this
sheet directly and does **not** infer samples by scanning `fastq_dir` at run time.
Because the sample sheet is a plain text file tracked alongside the code in version
control, the exact set of analysed samples — and any change to it — is permanently
recorded, so a run can be reproduced and inspected later. This also prevents a common
failure mode in which a stray, hidden, or partially downloaded file in the input
directory silently alters which samples enter the analysis.

```tsv
sample	unit	fq1	fq2
SRR28074879	1	data/GSE256519/SRR28074879_1.fastq	data/GSE256519/SRR28074879_2.fastq
SRR28074880	1	data/GSE256519/SRR28074880_1.fastq	data/GSE256519/SRR28074880_2.fastq
```

One row per sample (the `unit` column is kept for snakemake-workflows template
compatibility; multi-lane merging is not wired in — pre-merge lanes if needed).
For single-end data, leave `fq2` empty.

Regenerate the sheet from a FASTQ directory (run once, not part of the pipeline):

```bash
bash workflow/scripts/make_units.sh data/GSE256519 paired   # writes config/units.tsv
```

---

## Usage Example

A complete run from a fresh clone to final results, using the GSE256519 test
dataset (human heart RNA-seq, paired-end 2×150 bp). The same steps work for your
own data — only Step 5 changes (see the note there).

> ⚠️ **On HPC, run all of this from a login node (e.g. `int*`), not a compute
> node.** `run.sh` and the download wrappers **self-submit their own SLURM jobs** —
> you only launch them. Local-workstation users can ignore every "HPC" note below.

### Step 1 — Clone the repository

```bash
git clone https://github.com/ManHUU/snakemake-gatk-rna-workflow.git
cd snakemake-gatk-rna-workflow
```

On HPC, cloning onto scratch (e.g. `/scratch-shared/$USER/`) is recommended — fast
disk and large quota for the hundreds of GB of intermediates. Mind the scratch
purge & backup caveat under [Running the Pipeline](#running-the-pipeline).

### Step 2 — Install Snakemake (pinned)

```bash
conda env create -f workflow/envs/snakemake.yaml   # first time only
conda activate snakemake_env
```

That is the only thing you install. Every scientific tool runs from a container,
and `apptainer`/`singularity` itself is already on `PATH` on most HPCs (including
Snellius) — if a container runtime is ever missing, `run.sh` stops with a clear
message telling you to `module load apptainer`.

### Step 3 — Set HPC settings (HPC only — local users skip)

Everything that varies between clusters is set here as environment variables; you
never edit the pipeline code. Set them once per shell (or in your `~/.bashrc`):

| Variable | Required? | What it is |
|---|---|---|
| `SLURM_ACCOUNT` | **Required** | Your SLURM account / billing budget. |
| `SLURM_PARTITION` | **Required** | The partition jobs run on (e.g. `genoa`). |
| `HPC_SCRATCH_DIR` | Recommended | Fast scratch where apptainer builds its multi-GB images. If set it always wins; if unset, `run.sh` auto-detects (`$SCRATCH`, `/scratch-shared/$USER`, `/scratch/$USER`) and prints the resolved path. |
| `EXTRA_BIND_PATHS` | Optional | Extra host dirs made visible *inside* containers — only if `fastq_dir`/`output_dir` point outside the repo (the repo is auto-bound). |

```bash
export SLURM_ACCOUNT=<your_account>
export SLURM_PARTITION=<your_partition>           # e.g. genoa
export HPC_SCRATCH_DIR=/scratch-shared/$USER      # recommended; Snellius example
```

The one per-site value that is **not** an env var is `partition_max_runtime` in
`config/config.yaml` (your partition's MaxTime in minutes — find it with
`scontrol show partition <name> | grep MaxTime`). It lives in config because
Snakemake, not `run.sh`, consumes it. The default (7200 = 5 days) suits Snellius
`genoa`; lower it only if your partition's limit is lower.

### Step 4 — Download reference files

```bash
bash workflow/scripts/download_refs.sh        # local
bash workflow/scripts/download_refs.slurm     # HPC (self-submits via sbatch)
```

References land under `resources/`; the paths in `config/config.yaml` already match —
no edits.

### Step 5 — Download input FASTQ data

The wrapper auto-detects local vs SLURM execution. The SRA Toolkit runs from a
pinned biocontainer pulled on first use — no extra setup.

```bash
bash workflow/scripts/download_GSE256519.slurm
```

After this you should have:
```
data/GSE256519/
├── SRR28074879_1.fastq
├── SRR28074879_2.fastq
├── SRR28074880_1.fastq
└── SRR28074880_2.fastq
```

`config/units.tsv` is pre-filled for this dataset — no edits for the walkthrough.

> **Using your own data instead?** Two changes, both in `config/` — never in code:
> 1. Regenerate the sample sheet: `bash workflow/scripts/make_units.sh <your_dir> <paired|single>`
>    For example, paired-end FASTQs in `data/my_run/` named
>    `tumorA_1.fastq` / `tumorA_2.fastq`, `tumorB_1.fastq` / `tumorB_2.fastq`:
>    ```bash
>    bash workflow/scripts/make_units.sh data/my_run paired
>    ```
>    writes `config/units.tsv`:
>    ```
>    sample  unit  fq1                          fq2
>    tumorA  1     data/my_run/tumorA_1.fastq   data/my_run/tumorA_2.fastq
>    tumorB  1     data/my_run/tumorB_1.fastq   data/my_run/tumorB_2.fastq
>    ```
>    (columns are tab-separated; single-end leaves `fq2` empty). The naming
>    convention is `{sample}_1.fastq` / `{sample}_2.fastq` for paired and
>    `{sample}.fastq` for single-end.
> 2. Set `sequencing_type: paired|single` in `config/config.yaml` to match.
>
> Review `config/units.tsv` before running — see [Configuration § Sample sheet](#sample-sheet-configunitstsv).

### Step 6 — Dry run (verify the DAG)

```bash
bash run.sh --dry-run
```

You should see all rules listed for both samples across 24 chromosomes, no errors.
On HPC this plans on the login node and submits nothing.

### Step 7 — Run the pipeline

```bash
bash run.sh
```

Locally this runs on your machine; on HPC it self-submits an orchestrator job that
scatters the per-rule jobs. To **resume** a previous run and reuse work already on
disk (instead of recomputing rules whose code changed), add:

```bash
bash run.sh --rerun-triggers mtime
```

### Step 8 — Monitor (HPC)

```bash
squeue -u $USER
tail -f logs/gatk_pipeline_<jobid>.log
```

### Step 9 — Check outputs

```bash
ls -lh results/Final_CDS_rsID.vcf.gz         # final variant file
xdg-open results/qc/multiqc_report.html      # QC report
xdg-open results/Final_CDS_healthCheck.html  # SnpEff health check
ls results/plots/                            # figures
```

### Step 10 (optional) — Generate a provenance report

After a successful run, build a self-contained `report.html` with the workflow DAG,
per-rule runtimes, and software/provenance metadata — handy as supplementary material:

```bash
bash run.sh --report report.html
```

Like `--unlock`, this runs locally in seconds (no sbatch/apptainer needed).

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

`run.sh` is the single entrypoint and **auto-detects how to run** — you call it
the same way everywhere. The HPC settings and the full step-by-step live in
[Usage Example](#usage-example); this section covers the execution modes and a
couple of details.

- **Local workstation** — runs directly on your machine, using 8 cores by default
  (change in `workflow/profiles/local/config.yaml`). For long runs, launch from a
  `screen` (or `tmux`) session so it survives disconnection:
  ```bash
  screen -S gatk_run
  bash run.sh
  # Detach with Ctrl+A then D; reattach later with `screen -r gatk_run`
  ```
- **HPC with SLURM** — `run.sh` self-submits an orchestrator job that scatters the
  per-rule jobs, binds the working directory into every container, and forwards
  your `SLURM_ACCOUNT` / `SLURM_PARTITION` (set in
  [Usage Example § Step 3](#step-3--set-hpc-settings-hpc-only--local-users-skip)).
  `--dry-run` always plans locally and submits nothing; `sbatch run.sh` is the
  equivalent explicit form.

> **Running on scratch (purge & backup caveat).** Cloning and running the repo
> directly on scratch (e.g. `/scratch-shared/$USER/...`) is fully supported and
> recommended — all paths are repo-relative, scratch is bound into the
> containers, and it has the fast disk and large quota this pipeline's
> hundreds of GB of intermediates need. The caveat is a property of scratch
> itself, not the workflow: **scratch is auto-purged and not backed up.** Most
> sites delete files untouched for a retention window (on Snellius
> `/scratch-shared` it is roughly 14 days — check your site's policy). For this
> pipeline that means:
> - Don't let a run idle past the purge window. If you must pause for days,
>   keep it running or move the working tree to project/home storage —
>   otherwise references, the apptainer image cache, or intermediate files can
>   be deleted mid-run (a *partially* purged file can even look present but be
>   truncated, causing confusing failures).
> - Keep a durable copy of the **references** elsewhere if you reuse them across
>   weeks; re-downloading is slow.
> - **Copy `results/` off scratch once the run completes** — scratch is not
>   backed up.

---

## Recovering from an interrupted run

Long HPC runs occasionally get interrupted — a SLURM wall-time limit hits, a
node fails, `scancel` arrives, or the local terminal hosting the driver
disconnects. Snakemake re-uses every intermediate file already on disk and
picks up exactly where the previous run left off, so recovery is usually two
commands.

### 1. Clear the stale lock

A clean Snakemake exit removes its workflow lock; an un-graceful exit leaves
it behind, and the next launch fails with:

```
LockException:
Error: Directory cannot be locked. Please make sure that no other Snakemake
process is trying to create the same files in the following directory:
<your-repo-path>
```

First check no other run is actually live:

```bash
squeue -u $USER                            # any pipeline jobs still pending/running?
ps -u $USER -o pid,cmd | grep -i snakemake # any driver process still attached?
```

If both come back empty, the lock is stale — clear it with:

```bash
bash run.sh --unlock
```

`run.sh` short-circuits the sbatch path for `--unlock`, so this runs in
seconds wherever you launched it (head node, login shell, anywhere
`snakemake` is on PATH) — no job queueing.

### 2. Find the cause (optional, but recommended)

```bash
# Previous orchestrator log:
ls -lt logs/gatk_pipeline_*.log | head -5
tail -100 logs/gatk_pipeline_<PREV_JOBID>.log

# Per-rule SLURM logs — drill here if the orchestrator only reports
# "Reason: Unknown" (the real failure line usually lives in the per-rule log):
ls -lt .snakemake/slurm_logs/rule_*/*/*.log | head
```

Common culprits, in rough order of frequency:

- **SLURM wall-time hit.** `run.sh` requests `--time=72:00:00`; if your cohort
  is large or your partition has tighter limits, bump the `#SBATCH --time=`
  line at the top of `run.sh`.
- **Out-of-memory on an outlier sample.** The five GATK per-sample rules
  (`mark_duplicates`, `split_n_cigar_reads`, `base_recalibrator`, `apply_BQSR`,
  `haplotype_caller`) auto-escalate `mem_mb` and `runtime` over three retry
  attempts; extreme inputs can still exceed the third attempt's ceiling.
- **Node failure.** Snakemake re-queues failed jobs per `restart-times` in
  `workflow/profiles/slurm/config.yaml`.
- **Local terminal disconnected.** Wrap long local runs in `screen` or `tmux`
  (see [Running the Pipeline](#running-the-pipeline)).

### 3. Resume

```bash
bash run.sh
```

The DAG resolver skips everything still on disk and starts at the first
missing output — typically far past the original failure point.

---

## HPC tip: where apptainer builds its container images

When apptainer pulls a container image it writes several GB of temporary
files to a scratch directory. Putting this on cluster scratch (rather than
on your project filesystem or in RAM) makes builds faster and avoids
project-quota and out-of-memory issues. `run.sh` resolves the scratch path
in this order — first match wins:

1. **Explicit `APPTAINER_TMPDIR`** — power-user override, used as-is unless
   it points at a RAM-backed filesystem (see safety guard below).
2. **`HPC_SCRATCH_DIR`** — if you export it, the pipeline uses
   `${HPC_SCRATCH_DIR}/apptainer-tmp` and `${HPC_SCRATCH_DIR}/apptainer-cache`.
3. **Auto-detect** — otherwise the pipeline probes common HPC conventions
   (`$SCRATCH`, `/scratch-shared/$USER`, `/scratch/$USER`) and picks the
   first writable, disk-backed path.
4. **In-repo fallback** — if nothing else works, falls back to
   `resources/containers/{tmp,cache}` inside the repo. Always works, but
   counts against your project disk quota on HPC.

Every run prints the chosen path so it is visible in your SLURM log:

```
Apptainer scratch : /scratch-shared/<user>/apptainer-tmp
```

**Recommended on HPC: export `HPC_SCRATCH_DIR` explicitly.** Autodetect is
usually fine, but exporting it gives you explicit control and a clear
record in the job log of where containers live:

```bash
export HPC_SCRATCH_DIR=/scratch-shared/$USER     # Snellius / SURF
export HPC_SCRATCH_DIR=$SCRATCH                  # TACC and sites that set $SCRATCH
export HPC_SCRATCH_DIR=/scratch/$USER            # many university clusters
```

If unsure, check your cluster's documentation for "scratch" or "temporary
storage", or inspect your environment with `env | grep -iE 'scratch|tmp'`.

**Safety guard**: if `APPTAINER_TMPDIR` is already set in your environment
(e.g. by `module load apptainer`) and it points at a RAM-backed filesystem
like Snellius's `/tmp`, the pipeline rejects it with a warning and
re-resolves — large image builds on tmpfs would OOM-kill the job.

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
| SnpEff | quay.io/biocontainers/snpeff | 5.3.0a |
| SnpSift | quay.io/biocontainers/snpsift | 5.3.0a |
| SRA Toolkit | quay.io/biocontainers/sra-tools | 3.1.1 |
| matplotlib | quay.io/biocontainers/matplotlib | 3.5.1 |

---

## Test Data

The pipeline was developed using human heart RNA-seq data from GEO dataset
GSE256519 (paired-end 2×150 bp). The bundled wrapper downloads accessions
SRR28074879 and SRR28074880 into `data/GSE256519/` (see
[Usage Example § Step 5](#step-5--download-input-fastq-data) for the full
walkthrough):

```bash
bash workflow/scripts/download_GSE256519.slurm
```

---

## License

See LICENSE file.
