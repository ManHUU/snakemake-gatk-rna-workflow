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

The following walks through a complete run from a fresh clone to final results,
using the GSE256519 test dataset (human heart RNA-seq, paired-end 2×150 bp).

### Step 1 — Clone the repository

```bash
git clone https://github.com/ManHUU/snakemake-gatk-rna-workflow.git
cd snakemake-gatk-rna-workflow
```

### Step 2 — Install Snakemake (pinned)

```bash
conda env create -f workflow/envs/snakemake.yaml
conda activate snakemake_env
```

### Step 3 — Download reference files

On a local workstation:
```bash
bash workflow/scripts/download_refs.sh
```

On HPC, set the two SLURM env vars once (same ones used by `run.sh` later),
then launch the wrapper — it self-submits via sbatch:
```bash
export SLURM_ACCOUNT=<your_slurm_account>
export SLURM_PARTITION=<your_slurm_partition>
bash workflow/scripts/download_refs.slurm
```

This places all references under `resources/` and the paths in `config/config.yaml`
match without any edits.

### Step 4 — Download input FASTQ data

The wrapper auto-detects local vs SLURM execution and reuses the same
`SLURM_ACCOUNT` / `SLURM_PARTITION` env vars from Step 3. The SRA Toolkit
runs from a pinned biocontainer pulled on first use — no conda env required.

```bash
bash workflow/scripts/download_GSE256519.slurm
```

> **Note:** This downloads the GSE256519 test dataset. To use your own data,
> place your FASTQ files in a directory of your choice, then regenerate the
> sample sheet with `bash workflow/scripts/make_units.sh <your_dir> <paired|single>`.

After this you should have:
```
data/GSE256519/
├── SRR28074879_1.fastq
├── SRR28074879_2.fastq
├── SRR28074880_1.fastq
└── SRR28074880_2.fastq
```

The sample sheet `config/units.tsv` is pre-filled for this test dataset, so for the
GSE256519 walkthrough no edits are needed. For your own data, (re)generate it with
`make_units.sh` and review it before running — see
[Configuration § Sample sheet](#sample-sheet-configunitstsv).

### Step 5 — Dry run (verify the DAG)

```bash
bash run.sh --dry-run
```

You should see all rules listed for both samples across 24 chromosomes with no errors.

### Step 6 — Run the pipeline

```bash
bash run.sh
```

`run.sh` runs locally on workstations, or self-submits via `sbatch` on HPC
when `SLURM_ACCOUNT` / `SLURM_PARTITION` are exported. See
[Running the Pipeline](#running-the-pipeline) below for full details
(screen sessions, `EXTRA_BIND_PATHS`, etc.).

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

### Step 8 (optional) — Generate a provenance report

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

A single entrypoint, `run.sh`, handles both modes via auto-detection.

### Local workstation

```bash
bash run.sh             # full run
bash run.sh --dry-run   # preview steps only
```

For long runs, launch from a `screen` (or `tmux`) session so the run survives
disconnection:
```bash
screen -S gatk_run
bash run.sh
# Detach with Ctrl+A then D; reattach later with `screen -r gatk_run`
```

The local profile uses 8 cores by default. Edit
`workflow/profiles/local/config.yaml` to change.

### HPC / SLURM

```bash
export SLURM_ACCOUNT=<your_account>
export SLURM_PARTITION=<your_partition>

bash run.sh --dry-run   # plans the DAG on the head node, submits nothing
bash run.sh             # self-submits an orchestrator job that scatters per-rule jobs
sbatch run.sh           # equivalent explicit form
```

Optionally point apptainer at your cluster's scratch (recommended on HPC —
see [HPC tip: where apptainer builds its container images](#hpc-tip-where-apptainer-builds-its-container-images)
for details and autodetect behavior):
```bash
export HPC_SCRATCH_DIR=/scratch-shared/$USER     # Snellius example
```

Optionally make extra host directories visible **inside** the running
containers. This is unrelated to `HPC_SCRATCH_DIR` above — it controls
what the container can see, not where images are built. You only need it
if your `fastq_dir` or `output_dir` in `config/config.yaml` points outside
the repository (the repo directory is auto-bound):
```bash
export EXTRA_BIND_PATHS=/scratch-shared,/tmp
```

`run.sh` auto-detects whether you are on a SLURM-capable machine, dynamically
binds the current working directory into all containers, and forwards the
SLURM account/partition to per-rule submissions.

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
  (see [Running the Pipeline § Local workstation](#local-workstation)).

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
[Usage Example § Step 4](#step-4--download-input-fastq-data) for the full
walkthrough):

```bash
bash workflow/scripts/download_GSE256519.slurm
```

---

## License

See LICENSE file.
