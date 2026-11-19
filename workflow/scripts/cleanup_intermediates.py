#!/usr/bin/env python3
"""
cleanup_intermediates.py  (reviewer comment C-13, uses the C-06 path guard)

A SEPARATE, EXPLICIT cleanup action - never called automatically by the
workflow. It deletes only an explicit, workflow-owned list of redundant
intermediate alignment files, inside the path-safety-validated output_dir, and
only AFTER you have confirmed the required outputs. Deletion is DISABLED by
default (preview/dry-run); pass --delete to actually remove the previewed files.

Eligible for cleanup (redundant earlier BAM stages, superseded once their
downstream outputs are validated):
    {sample}.sorted.bam(.bai)           -> superseded by markdup
    {sample}.markdup.bam(.bai)          -> superseded by markdup.filtered
    {sample}.markdup.filtered.bam(.bai) -> superseded by split
    {sample}.split.bam/.split.bai       -> superseded by BQSR
    trimmed/{sample}_*.trimmed.fastq    -> regenerable from the archived FASTQs

RETAINED by default (NOT touched here):
    {sample}.BQSR.bam/.bai             final BAM used for variant calling
    {sample}.g.vcf.gz(.tbi)            per-sample gVCFs for repeat joint calling
    Joint_all / Final_PASS / Final_PASS_noEdit / Final_clean / Final_annotated /
    Final_CDS / Final_CDS_rsID (+ indexes), stats/, plots/, qc/, logs, metrics.

Samples are read from the sample sheet (config: units, default config/units.tsv),
the single source of truth used by the workflow — this does NOT glob fastq_dir.

Usage:
    python -I workflow/scripts/cleanup_intermediates.py            # preview
    python -I workflow/scripts/cleanup_intermediates.py --delete   # delete
    python -I workflow/scripts/cleanup_intermediates.py --config config/config.yaml
"""

import argparse
import csv
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pathsafe  # noqa: E402


def discover_samples(units_tsv):
    """Read the `sample` column of the sample sheet (one row per sample)."""
    samples = []
    with open(units_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s = (row.get("sample") or "").strip()
            if s:
                samples.append(s)
    return sorted(set(samples))


def owned_intermediates(samples, paired):
    files = []
    for s in samples:
        files += [
            f"{s}.sorted.bam", f"{s}.sorted.bam.bai",
            f"{s}.markdup.bam", f"{s}.markdup.bam.bai",
            f"{s}.markdup.filtered.bam", f"{s}.markdup.filtered.bam.bai",
            f"{s}.split.bam", f"{s}.split.bai",
        ]
        if paired:
            files += [f"trimmed/{s}_1.trimmed.fastq", f"trimmed/{s}_2.trimmed.fastq"]
        else:
            files += [f"trimmed/{s}.trimmed.fastq"]
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--delete", action="store_true",
                    help="actually delete the previewed files (default: preview only)")
    args = ap.parse_args()

    project_root = os.getcwd()
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    output_dir = cfg["output_dir"]
    paired = cfg.get("sequencing_type", "paired") == "paired"
    units_tsv = cfg.get("units", "config/units.tsv")

    # Re-validate output_dir with the same guard the workflow uses at startup.
    out = pathsafe.resolve_output_dir(output_dir, project_root)
    allowed_root = str(out.parent)

    samples = discover_samples(units_tsv)
    if not samples:
        print(f"No samples found in {units_tsv}; nothing to do.")
        return
    owned = owned_intermediates(samples, paired)

    targets = pathsafe.cleanup_files(output_dir, allowed_root, project_root,
                                     owned, dry_run=not args.delete)
    action = "Deleted" if args.delete else "Would delete (preview)"
    if not targets:
        print("No eligible intermediate files are present.")
        return
    total = 0
    for t in targets:
        try:
            total += t.stat().st_size if not args.delete else 0
        except OSError:
            pass
        print(f"  {action}: {t}")
    if not args.delete:
        print(f"\n{len(targets)} file(s), ~{total/1e9:.2f} GB. "
              "Re-run with --delete to remove them.")
    else:
        print(f"\n{len(targets)} file(s) deleted.")


if __name__ == "__main__":
    main()
