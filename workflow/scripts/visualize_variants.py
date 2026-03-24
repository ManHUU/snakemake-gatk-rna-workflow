#!/usr/bin/env python3
"""
visualize_variants.py
---------------------
Generates publication-quality figures from GATK variant-calling statistics.

Input TSV files (produced by the variant_statistics Snakemake rule via bcftools):
    --metrics   quality_metrics.tsv  – TYPE, FILTER, QD, FS, MQ, DP per variant
    --filters   filter_summary.tsv   – variant count per filter flag
    --chroms    chrom_counts.tsv     – PASS variant counts per chromosome and type
    --samples   sample_counts.tsv    – per-sample SNP / INDEL counts (PASS only)
    --tstv      ts_tv.tsv            – Ts, Tv counts and ratio before/after filtering

Output PNGs written to --outdir:
    1. variant_type_summary.png   – SNP / INDEL counts, PASS vs filtered (log-scale stacked bar)
    2. filter_summary.png         – PASS vs Filtered counts (log-scale horizontal bar)
    3. quality_distributions.png  – QD and DP histograms by variant type (DP on log x-scale)
    4. chrom_distribution.png     – per-chromosome variant counts (stacked bar)
    5. sample_counts.png          – per-sample SNP / INDEL counts (stacked bar)
    6. ts_tv_comparison.png       – Ts/Tv counts (log scale) and ratio side-by-side

Usage:
    python visualize_variants.py \
        --metrics quality_metrics.tsv \
        --filters filter_summary.tsv \
        --chroms  chrom_counts.tsv \
        --samples sample_counts.tsv \
        --tstv    ts_tv.tsv \
        --outdir  results/plots
"""

import argparse
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


# ── Nature / Cell journal style ───────────────────────────────────────────────

# Okabe-Ito colorblind-safe palette (doi:10.1038/nmeth.1618)
PALETTE = {
    "SNP":      "#0072B2",   # blue
    "INDEL":    "#D55E00",   # vermilion
    "PASS":     "#009E73",   # bluish green
    "Filtered": "#999999",   # gray
    "Ts":       "#56B4E9",   # sky blue
    "Tv":       "#E69F00",   # orange
    "ratio":    "#000000",   # black
}

plt.rcParams.update({
    "font.family":           "sans-serif",
    "font.sans-serif":       ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":             10,
    "axes.spines.top":       False,
    "axes.spines.right":     False,
    "axes.linewidth":        0.8,
    "axes.labelsize":        10,
    "axes.titlesize":        11,
    "axes.titleweight":      "bold",
    "xtick.major.width":     0.8,
    "ytick.major.width":     0.8,
    "xtick.labelsize":       9,
    "ytick.labelsize":       9,
    "legend.fontsize":       8,
    "legend.frameon":        False,
    "figure.facecolor":      "white",
    "axes.facecolor":        "white",
    "savefig.facecolor":     "white",
})

FIG_DPI = 300
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def fmt_count(v, _):
    """Tick formatter: display large numbers with K/M suffixes."""
    if v >= 1_000_000:
        return f"{v/1_000_000:.1f}M"
    if v >= 1_000:
        return f"{v/1_000:.0f}K"
    return f"{int(v)}"


def save_fig(fig, outdir, filename):
    for ext in ("png", "pdf"):
        path = os.path.join(outdir, filename.replace(".png", f".{ext}"))
        fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
        print(f"  Saved: {path}")
    plt.close(fig)


# ── Plot 1: Variant type summary ──────────────────────────────────────────────

def plot_type_summary(filter_tsv: str, outdir: str):
    """Log-scale stacked bar: SNP and INDEL counts split by PASS vs Filtered."""
    counts = defaultdict(lambda: {"PASS": 0, "Filtered": 0})

    with open(filter_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            vtype, status = parts[0].upper(), parts[1]
            if vtype not in ("SNP", "INDEL") or status not in ("PASS", "Filtered"):
                continue
            try:
                counts[vtype][status] = int(parts[2])
            except ValueError:
                continue

    types  = ["SNP", "INDEL"]
    passes = [counts[t]["PASS"]     for t in types]
    fails  = [counts[t]["Filtered"] for t in types]

    fig, ax = plt.subplots(figsize=(5, 4))
    x = range(len(types))

    ax.bar(x, passes, color=PALETTE["PASS"],     label="PASS",     width=0.5)
    ax.bar(x, fails,  color=PALETTE["Filtered"], label="Filtered", width=0.5,
           bottom=passes)

    ax.set_yscale("log")
    ax.set_xticks(list(x))
    ax.set_xticklabels(types)
    ax.set_ylabel("Variant count (log scale)")
    ax.set_title("Variant counts by type and filter status")
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt_count))

    # Annotate PASS count inside bar and total on top
    for i, (p, f) in enumerate(zip(passes, fails)):
        if p > 0:
            ax.text(i, p * 0.5, f"PASS\n{p:,}", ha="center", va="center",
                    fontsize=7, color="white", fontweight="bold")
        ax.text(i, p + f, f"{p + f:,}", ha="center", va="bottom", fontsize=8)

    save_fig(fig, outdir, "variant_type_summary.png")


# ── Plot 2: Filter summary ────────────────────────────────────────────────────

def plot_filter_summary(filter_tsv: str, outdir: str):
    """Log-scale horizontal bar: total PASS vs Filtered variants."""
    total_pass = 0
    total_filtered = 0

    with open(filter_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            try:
                status, count = parts[1], int(parts[2])
                if status == "PASS":
                    total_pass += count
                elif status == "Filtered":
                    total_filtered += count
            except ValueError:
                continue

    if total_pass == 0 and total_filtered == 0:
        print("  Warning: filter_summary.tsv has no data – skipping plot 2.")
        return

    labels = ["PASS", "Filtered"]
    values = [total_pass, total_filtered]
    colors = [PALETTE["PASS"], PALETTE["Filtered"]]

    fig, ax = plt.subplots(figsize=(7, 3))
    bars = ax.barh([0, 1], values, color=colors, height=0.5)
    ax.set_xscale("log")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(labels)
    ax.set_xlabel("Variant count (log scale)")
    ax.set_title("Filter outcome: PASS vs Filtered (all types combined)")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_count))

    for i, v in enumerate(values):
        ax.text(v * 1.05, i, f"{v:,}", va="center", fontsize=9)

    plt.tight_layout()
    save_fig(fig, outdir, "filter_summary.png")


# ── Plot 3: Quality metric distributions ─────────────────────────────────────

def plot_quality_distributions(metrics_tsv: str, outdir: str):
    """
    2×2 grid of histograms: rows = SNP / INDEL, cols = QD, DP.
    FS and MQ removed. DP uses log x-scale to handle the wide range.
    Red dashed lines mark the hard-filter thresholds.
    """
    METRIC_COL   = {"QD": 2, "DP": 5}
    METRIC_LABEL = {
        "QD": "Quality by Depth (QD)",
        "DP": "Read Depth (DP, log$_{10}$)",
    }
    THRESHOLDS = {
        "QD": {"SNP": 2.0,  "INDEL": 2.0},
        "DP": {"SNP": 10,   "INDEL": 10},
    }
    DP_LOG = True   # log x-scale for read depth

    data: dict = {
        "SNP":   {m: [] for m in METRIC_COL},
        "INDEL": {m: [] for m in METRIC_COL},
    }

    with open(metrics_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            vtype = parts[0].upper()
            if vtype not in data:
                continue
            for metric, col in METRIC_COL.items():
                raw = parts[col]
                if raw not in (".", ""):
                    try:
                        val = float(raw)
                        if val > 0:      # log-safe
                            data[vtype][metric].append(val)
                    except ValueError:
                        pass

    vtypes  = ["SNP", "INDEL"]
    metrics = list(METRIC_COL.keys())

    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    fig.suptitle("Quality metric distributions (pre-filter variants)", fontsize=12,
                 fontweight="bold")

    for row, vtype in enumerate(vtypes):
        for col, metric in enumerate(metrics):
            ax = axes[row][col]
            values = data[vtype][metric]

            if not values:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes, color="grey")
                ax.set_title(f"{vtype} – {metric}", fontsize=9)
                continue

            use_log_x = (metric == "DP" and DP_LOG)

            if use_log_x:
                log_vals = np.log10(values)
                bins = np.linspace(log_vals.min(), min(log_vals.max(), 6), 61)
                ax.hist(values, bins=10**bins,
                        color=PALETTE[vtype], alpha=0.85,
                        edgecolor="white", linewidth=0.3)
                ax.set_xscale("log")
                ax.set_xlim(left=1)
                ax.xaxis.set_major_formatter(ticker.FuncFormatter(fmt_count))
            else:
                ax.hist(values, bins=60, color=PALETTE[vtype],
                        alpha=0.85, edgecolor="white", linewidth=0.3)

            thresh = THRESHOLDS[metric][vtype]
            if thresh is not None:
                ax.axvline(thresh, color="#D55E00", linestyle="--",
                           linewidth=1.4, label=f"threshold = {thresh}")
                ax.legend(fontsize=7, loc="upper right")

            ax.set_title(f"{vtype} – {METRIC_LABEL[metric]}", fontsize=10)
            ax.set_xlabel("Read Depth (log$_{10}$ scale)" if use_log_x else "")
            ax.set_ylabel("Count" if col == 0 else "")
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt_count))

    plt.tight_layout()
    save_fig(fig, outdir, "quality_distributions.png")


# ── Plot 4: Per-chromosome distribution ──────────────────────────────────────

def plot_chrom_distribution(chrom_tsv: str, outdir: str):
    """Stacked bar of PASS SNP and INDEL counts per chromosome."""
    chrom_data: dict = defaultdict(lambda: defaultdict(int))

    with open(chrom_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            chrom, vtype, count = parts[0], parts[1].upper(), int(parts[2])
            # collapse OVERLAP subtypes (SNP,OVERLAP → SNP)
            vtype = vtype.split(",")[0]
            chrom_data[chrom][vtype] += count

    chroms = [c for c in CHROM_ORDER if c in chrom_data]
    if not chroms:
        print("  Warning: chrom_counts.tsv is empty – skipping plot 4.")
        return

    snp_counts   = [chrom_data[c].get("SNP",   0) for c in chroms]
    indel_counts = [chrom_data[c].get("INDEL", 0) for c in chroms]

    fig, ax = plt.subplots(figsize=(13, 4))
    x = range(len(chroms))
    ax.bar(x, snp_counts,   color=PALETTE["SNP"],   label="SNP",   width=0.7)
    ax.bar(x, indel_counts, color=PALETTE["INDEL"], label="INDEL", width=0.7,
           bottom=snp_counts)
    ax.set_xticks(list(x))
    ax.set_xticklabels([c.replace("chr", "") for c in chroms],
                       rotation=0, fontsize=9)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Variant count (PASS)")
    ax.set_title("Per-chromosome variant distribution")
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt_count))

    plt.tight_layout()
    save_fig(fig, outdir, "chrom_distribution.png")


# ── Plot 5: Per-sample variant counts ────────────────────────────────────────

def plot_sample_counts(sample_tsv: str, outdir: str):
    """Stacked bar of PASS SNP and INDEL counts per sample."""
    samples, snps, indels = [], [], []

    with open(sample_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                samples.append(parts[0])
                snps.append(int(parts[1]))
                indels.append(int(parts[2]))
            except ValueError:
                continue

    if not samples:
        print("  Warning: sample_counts.tsv is empty – skipping plot 5.")
        return

    MAX_LABEL = 25
    labels = [s if len(s) <= MAX_LABEL else s[:MAX_LABEL] + "…" for s in samples]

    x = range(len(samples))
    fig, ax = plt.subplots(figsize=(max(5, len(samples) * 1.6), 4))
    ax.bar(x, snps,   color=PALETTE["SNP"],   label="SNP",   width=0.5)
    ax.bar(x, indels, color=PALETTE["INDEL"], label="INDEL", width=0.5,
           bottom=snps)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Variant count (coding, PASS)")
    ax.set_title("Per-sample coding variant counts")
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt_count))

    for i, (s, d) in enumerate(zip(snps, indels)):
        ax.text(i, s + d, f"{s + d:,}", ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    save_fig(fig, outdir, "sample_counts.png")


# ── Plot 6: Ts/Tv comparison ──────────────────────────────────────────────────

def plot_ts_tv_comparison(tstv_tsv: str, outdir: str):
    """
    Two-panel figure:
      Left  — Grouped bars of Ts and Tv counts (log y-scale) pre- and post-filter.
      Right — Bar chart of Ts/Tv ratio with expected RNA-seq range annotated.

    Expected Ts/Tv ranges:
        Genome-wide SNPs:   ~2.1
        Coding / exonic:    ~2.8 – 3.3
        RNA-seq variants:   ~3.5 – 4.0
    """
    stages, ts_vals, tv_vals, ratios = [], [], [], []

    with open(tstv_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in ("stage", ""):
                continue
            try:
                stages.append(parts[0])
                ts_vals.append(int(parts[1]))
                tv_vals.append(int(parts[2]))
                ratios.append(float(parts[3]))
            except (ValueError, IndexError):
                continue

    if not stages:
        print("  Warning: ts_tv.tsv is empty – skipping plot 6.")
        return

    stage_labels = [s.replace("pre_filter", "Pre-filter")
                     .replace("post_filter", "Post-filter") for s in stages]

    fig, (ax_counts, ax_ratio) = plt.subplots(1, 2, figsize=(10, 4))

    # ── Left panel: Ts / Tv counts (log scale) ────────────────────────────────
    x     = np.arange(len(stages))
    width = 0.35

    ax_counts.bar(x - width / 2, ts_vals, width=width,
                  color=PALETTE["Ts"], label="Transitions (Ts)", alpha=0.9)
    ax_counts.bar(x + width / 2, tv_vals, width=width,
                  color=PALETTE["Tv"], label="Transversions (Tv)", alpha=0.9)

    ax_counts.set_yscale("log")
    ax_counts.set_xticks(x)
    ax_counts.set_xticklabels(stage_labels, fontsize=10)
    ax_counts.set_ylabel("Variant count (log scale)")
    ax_counts.set_title("Ts and Tv counts")
    ax_counts.legend()
    ax_counts.yaxis.set_major_formatter(ticker.FuncFormatter(fmt_count))

    # Annotate count values above bars
    for i, (ts, tv) in enumerate(zip(ts_vals, tv_vals)):
        ax_counts.text(i - width / 2, ts * 1.1, fmt_count(ts, None),
                       ha="center", va="bottom", fontsize=8)
        ax_counts.text(i + width / 2, tv * 1.1, fmt_count(tv, None),
                       ha="center", va="bottom", fontsize=8)

    # ── Right panel: Ts/Tv ratio ───────────────────────────────────────────────
    bar_colors = [PALETTE["Ts"] if i == 0 else PALETTE["Tv"] for i in range(len(stages))]
    ax_ratio.bar(x, ratios, color=bar_colors, width=0.45, alpha=0.9)

    # Expected range bands
    ax_ratio.axhspan(3.5, 4.0, color="#009E73", alpha=0.15,
                     label="Expected RNA-seq (3.5–4.0)")
    ax_ratio.axhspan(2.8, 3.3, color="#E69F00", alpha=0.15,
                     label="Expected coding (2.8–3.3)")

    ax_ratio.set_xticks(x)
    ax_ratio.set_xticklabels(stage_labels, fontsize=10)
    ax_ratio.set_ylabel("Ts/Tv ratio")
    ax_ratio.set_title("Ts/Tv ratio")
    ax_ratio.set_ylim(0, max(ratios) * 1.4)
    ax_ratio.legend(loc="upper right")

    # Annotate ratio values
    for i, r in enumerate(ratios):
        ax_ratio.text(i, r + max(ratios) * 0.04, f"{r:.2f}",
                      ha="center", va="bottom", fontsize=11, fontweight="bold")

    plt.tight_layout()
    save_fig(fig, outdir, "ts_tv_comparison.png")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--metrics",  required=True,
                        help="quality_metrics.tsv from bcftools query")
    parser.add_argument("--filters",  required=True,
                        help="filter_summary.tsv from bcftools query")
    parser.add_argument("--chroms",   required=True,
                        help="chrom_counts.tsv from bcftools query")
    parser.add_argument("--samples",  required=True,
                        help="sample_counts.tsv from bcftools stats")
    parser.add_argument("--tstv",     required=True,
                        help="ts_tv.tsv from bcftools stats (pre and post filter)")
    parser.add_argument("--outdir",   required=True,
                        help="Output directory for PNG files")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("Generating variant visualizations...")

    plot_type_summary(args.filters, args.outdir)
    plot_filter_summary(args.filters, args.outdir)
    plot_quality_distributions(args.metrics, args.outdir)
    plot_chrom_distribution(args.chroms, args.outdir)
    plot_sample_counts(args.samples, args.outdir)
    plot_ts_tv_comparison(args.tstv, args.outdir)

    print(f"Done. All plots written to: {args.outdir}")


if __name__ == "__main__":
    main()
