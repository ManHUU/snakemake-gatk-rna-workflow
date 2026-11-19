import os
import sys
import csv
from snakemake.utils import validate
from snakemake.exceptions import WorkflowError


# 1. Load the configurations
# We point to the specific paths in your workflow directory
configfile: "config/config.yaml"
# Fail fast at DAG-build time on a mistyped/missing config field.
validate(config, "workflow/schemas/config.schema.yaml")

# Path-safety check on output_dir (reviewer comment C-06): reject empty, '/',
# the home directory, the repo root, parent traversal or symlinked paths BEFORE
# any job is submitted or any rm -rf cleanup runs. Complements the schema check.
PROJECT_ROOT = os.path.dirname(os.path.abspath(workflow.snakefile))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "workflow", "scripts"))
import pathsafe  # noqa: E402
try:
    _validated_output = pathsafe.resolve_output_dir(config["output_dir"], PROJECT_ROOT)
except ValueError as exc:
    raise WorkflowError(
        f"Refusing to start: unsafe output_dir={config['output_dir']!r} ({exc}). "
        "Set output_dir to a dedicated run directory (relative to the repo, or "
        "an absolute path on scratch); it may not be empty, '/', your home "
        "directory, or the repository root."
    )


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


# 3. Sequencing mode (set in config: sequencing_type: "paired" | "single")
IS_PAIRED = config.get("sequencing_type", "paired") == "paired"

# 4. Sample sheet – samples and FASTQ paths come from config/units.tsv
#   (columns: sample, unit, fq1, fq2). This replaces runtime globbing so the
#   exact set of analysed samples is explicit, version-controlled and
#   schema-validated. One row per sample (the `unit` column is kept for
#   snakemake-workflows template compatibility; multi-lane merging is not wired
#   in). Regenerate the sheet with workflow/scripts/make_units.sh.
UNITS_TSV = config.get("units", "config/units.tsv")

UNITS = {}
with open(UNITS_TSV, newline="") as _fh:
    for _row in csv.DictReader(_fh, delimiter="\t"):
        _row = {k: (v.strip() if isinstance(v, str) else v) for k, v in _row.items()}
        validate(_row, "workflow/schemas/units.schema.yaml")
        _sample = _row["sample"]
        if _sample in UNITS:
            raise ValueError(
                f"Duplicate sample '{_sample}' in {UNITS_TSV}: this pipeline "
                "expects one row per sample (multi-lane merging is not wired in)."
            )
        if IS_PAIRED and not _row.get("fq2"):
            raise ValueError(
                f"sequencing_type is 'paired' but sample '{_sample}' has no fq2 "
                f"in {UNITS_TSV}."
            )
        UNITS[_sample] = _row

SAMPLES = sorted(UNITS)


# FASTQ path lookups (sample sheet is the single source of truth).
def fq1_for(sample):
    return UNITS[sample]["fq1"]


def fq2_for(sample):
    return UNITS[sample]["fq2"]


# --- STAR read-length / sjdbOverhang resolution (reviewer comment C-09) -------
# The STAR index must be built with --sjdbOverhang ~ (max read length - 1).
# Rather than force every user to compute and hard-code this, we detect it from
# the actual input reads when config star.sjdb_overhang is 'auto' (the default)
# or null, and fall back to STAR's own default (100) if detection fails. An
# explicit integer is honoured as-is but cross-checked against the reads, so a
# stale or mistyped value surfaces as a warning instead of a silent wrong index.
def _detect_read_length(fastq_path, max_reads=10000):
    """Max read length over the first `max_reads` records of a FASTQ (plain or
    .gz), or None if it cannot be read."""
    if not fastq_path:
        return None
    import gzip
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as _fh:
            max_len = 0
            seen = 0
            for i, line in enumerate(_fh):
                if i % 4 == 1:  # sequence line of each 4-line FASTQ record
                    max_len = max(max_len, len(line.strip()))
                    seen += 1
                    if seen >= max_reads:
                        break
            return max_len or None
    except OSError:
        return None


def resolve_sjdb_overhang():
    """Resolve STAR --sjdbOverhang from config, detecting read length for
    'auto'/null and warning on an explicit value that disagrees with the reads."""
    configured = config.get("star", {}).get("sjdb_overhang")
    probe = fq1_for(SAMPLES[0]) if SAMPLES else None
    detected = _detect_read_length(probe)

    is_auto = configured is None or (
        isinstance(configured, str) and configured.strip().lower() == "auto"
    )
    if is_auto:
        if detected:
            print(f"[STAR] sjdb_overhang=auto: detected read length {detected} bp "
                  f"from {probe} -> --sjdbOverhang {detected - 1}", file=sys.stderr)
            return detected - 1
        print(f"[STAR] sjdb_overhang=auto but read length could not be detected "
              f"(probe={probe!r}); using STAR's default 100.", file=sys.stderr)
        return 100

    value = int(configured)
    if detected and (detected - 1) != value:
        print(f"[STAR] WARNING: config star.sjdb_overhang={value} but the first "
              f"sample's reads look like {detected} bp (expected {detected - 1}). "
              f"Using {value} as configured — if this is wrong, set "
              f"sjdb_overhang: auto (or {detected - 1}) and rebuild the index in a "
              f"fresh reference.star_index_dir.", file=sys.stderr)
    else:
        print(f"[STAR] sjdb_overhang={value} (explicit).", file=sys.stderr)
    return value


# 5. Target Rule
rule all:
    input:
        # Reference indices
        REF_FASTA + ".fai",
        # QC
        os.path.join(OUTPUT_DIR, "qc", "multiqc_report.html"),
        # Variant calling – final CDS-filtered VCF with rsIDs
        os.path.join(OUTPUT_DIR, "Final_CDS_rsID.vcf.gz"),
        # SnpEff health-check reports
        os.path.join(OUTPUT_DIR, "Final_annotated_healthCheck.html"),
        os.path.join(OUTPUT_DIR, "Final_CDS_healthCheck.html"),
        # Statistics
        os.path.join(OUTPUT_DIR, "stats", "quality_metrics.tsv"),
        os.path.join(OUTPUT_DIR, "stats", "ts_tv.tsv"),
        # Visualization
        os.path.join(OUTPUT_DIR, "plots")


# 6. Include QC rules
include: "workflow/rules/qc.smk"

# 7. Include the variant-calling rules
include: "workflow/rules/variant_calling_jointCall_per_Chr.smk"

# 8. Include filtering, statistics and visualization rules
include: "workflow/rules/filter_and_visualize.smk"


# Optional: On-start/On-error messaging
onsuccess:
    print(f"""
============================================================
✅  Workflow complete.

Key outputs (under {OUTPUT_DIR}/):
  • Final_CDS_rsID.vcf.gz            final CDS-filtered, rsID-annotated variants
  • Final_annotated_healthCheck.html / Final_CDS_healthCheck.html   SnpEff reports
  • qc/multiqc_report.html           read/alignment QC
  • stats/quality_metrics.tsv, stats/ts_tv.tsv, stats/chrom_counts.tsv …
  • stats/*_filter_counts.tsv        variants removed by each filter (before/after)
  • stats/figure3_sources.tsv        which VCF feeds each figure panel
  • plots/                           publication-quality figures

Suggested next steps:
  1. Review plots/ and stats/quality_metrics.tsv.
  2. Once you have confirmed the outputs, reclaim disk from the redundant
     intermediate BAMs/FASTQs (final BAMs and gVCFs are kept):
        python workflow/scripts/cleanup_intermediates.py            # preview
        python workflow/scripts/cleanup_intermediates.py --delete   # remove
============================================================
""")

onerror:
    print(f"""
============================================================
❌  Workflow failed. The actionable error is usually in the PER-RULE log,
    not in this orchestrator message:
        ls -t {LOG_DIR}/ | head          # newest log first
    Re-check host readiness (tools, disk, writable paths, SLURM):
        bash workflow/scripts/check_prerequisites.sh
============================================================
""")
