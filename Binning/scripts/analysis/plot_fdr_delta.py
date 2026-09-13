#!/usr/bin/env python3
"""
plot_fdr_delta_by_interval.py

Generate interval summaries and plots for classic-vs-binned FDR deltas.

Input
-----
classic_vs_binned_scan_peptide_fdr_wide.tsv from
compare_classic_vs_binned_fdr_summary.py.

Outputs
-------
fdr_delta_interval_summary.tsv
fdr_delta_by_base_interval.tsv
fdr_threshold_transition_summary.tsv
classic_vs_binned_fdr_wide_with_intervals.tsv
plots/hist_<delta>.png
plots/boxplot_<delta>_by_BaseFDRInterval.png
plots/boxplot_all_delta_stats_by_BaseFDRInterval.png

By default, both median and mean binned-FDR deltas are plotted automatically.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_edges(s: str) -> List[float]:
    vals = []
    for part in s.split(","):
        part = part.strip()
        if part:
            vals.append(float(part))
    if len(vals) < 2:
        raise ValueError("At least two interval edges are required")
    vals = sorted(set(vals))
    if vals[0] > 0:
        vals = [0.0] + vals
    if vals[-1] < 1:
        vals = vals + [1.0]
    return vals


def interval_labels(edges: List[float]) -> List[str]:
    return [f"{a:g}-{b:g}" for a, b in zip(edges[:-1], edges[1:])]


def add_interval(df: pd.DataFrame, value_col: str, edges: List[float], out_col: str) -> pd.DataFrame:
    df[out_col] = pd.cut(
        df[value_col],
        bins=edges,
        labels=interval_labels(edges),
        include_lowest=True,
        right=True,
    )
    return df


def safe_skew(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    return float(x.skew()) if len(x) >= 3 else np.nan


def safe_kurtosis(x: pd.Series) -> float:
    x = pd.to_numeric(x, errors="coerce").dropna()
    return float(x.kurtosis()) if len(x) >= 4 else np.nan


def normality_tests(x: pd.Series, max_n: int = 5000) -> dict:
    x = pd.to_numeric(x, errors="coerce").dropna()
    if len(x) < 8:
        return {
            "NormalityTest": "not_run_n_lt_8",
            "NormalityN": len(x),
            "NormalityStatistic": np.nan,
            "NormalityPValue": np.nan,
        }
    if len(x) > max_n:
        x = x.sample(max_n, random_state=1)
    try:
        from scipy import stats  # type: ignore
        if len(x) < 20:
            stat, p = stats.shapiro(x)
            test = "Shapiro_sampled"
        else:
            stat, p = stats.normaltest(x)
            test = "DAgostinoK2_sampled"
        return {
            "NormalityTest": test,
            "NormalityN": len(x),
            "NormalityStatistic": float(stat),
            "NormalityPValue": float(p),
        }
    except Exception as e:
        return {
            "NormalityTest": f"not_run_{type(e).__name__}",
            "NormalityN": len(x),
            "NormalityStatistic": np.nan,
            "NormalityPValue": np.nan,
        }


def summarize_delta(df: pd.DataFrame, group_col: str, delta_col: str, peptide_col: str) -> pd.DataFrame:
    rows = []
    for group_val, sub in df.groupby(group_col, dropna=False, observed=False):
        x = pd.to_numeric(sub[delta_col], errors="coerce").dropna()
        if x.empty:
            continue
        nt = normality_tests(x)
        rows.append({
            "IntervalBasis": group_col,
            "Interval": str(group_val),
            "N_ScanPeptideKeys": int(len(sub)),
            "N_UniquePeptides": int(sub[peptide_col].nunique(dropna=True)) if peptide_col in sub.columns else np.nan,
            "DeltaColumn": delta_col,
            "DeltaMean": float(x.mean()),
            "DeltaMedian": float(x.median()),
            "DeltaStd": float(x.std(ddof=1)) if len(x) > 1 else 0.0,
            "DeltaMin": float(x.min()),
            "DeltaQ01": float(x.quantile(0.01)),
            "DeltaQ05": float(x.quantile(0.05)),
            "DeltaQ25": float(x.quantile(0.25)),
            "DeltaQ75": float(x.quantile(0.75)),
            "DeltaQ95": float(x.quantile(0.95)),
            "DeltaQ99": float(x.quantile(0.99)),
            "DeltaMax": float(x.max()),
            "Skew": safe_skew(x),
            "ExcessKurtosis": safe_kurtosis(x),
            "FractionNegative": float((x < 0).mean()),
            "FractionPositive": float((x > 0).mean()),
            "FractionZero": float((x == 0).mean()),
            **nt,
        })
    return pd.DataFrame(rows)


def make_histogram(df: pd.DataFrame, delta_col: str, outpath: Path, title: str) -> None:
    x = pd.to_numeric(df[delta_col], errors="coerce").dropna()
    if x.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(x, bins=120)
    ax.axvline(0, linewidth=1)
    ax.axvline(x.mean(), linewidth=1, linestyle="--")
    ax.axvline(x.median(), linewidth=1, linestyle=":")
    ax.set_xlabel(delta_col)
    ax.set_ylabel("ScanNum/CorePeptide count")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def make_interval_boxplot(df: pd.DataFrame, interval_col: str, delta_col: str, outpath: Path, title: str) -> None:
    plot_df = df[[interval_col, delta_col]].copy()
    plot_df[delta_col] = pd.to_numeric(plot_df[delta_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[interval_col, delta_col])
    if plot_df.empty:
        return

    groups = []
    labels = []
    for label, sub in plot_df.groupby(interval_col, observed=False):
        vals = sub[delta_col].dropna().to_numpy()
        if len(vals):
            groups.append(vals)
            labels.append(str(label))
    if not groups:
        return

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.1), 5))
    ax.boxplot(groups, labels=labels, showfliers=False)
    ax.axhline(0, linewidth=1)
    ax.set_xlabel(interval_col)
    ax.set_ylabel(delta_col)
    ax.set_title(title)
    plt.xticks(rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def make_combined_delta_boxplot(df: pd.DataFrame, interval_col: str, delta_cols: List[str], outpath: Path) -> None:
    frames = []
    for col in delta_cols:
        if col not in df.columns:
            continue
        tmp = df[[interval_col, col]].copy()
        tmp = tmp.rename(columns={col: "DeltaFDR"})
        tmp["DeltaStatistic"] = col.replace("DeltaBinned", "").replace("MinusBase", "")
        tmp["DeltaFDR"] = pd.to_numeric(tmp["DeltaFDR"], errors="coerce")
        tmp = tmp.dropna(subset=[interval_col, "DeltaFDR"])
        frames.append(tmp)
    if not frames:
        return
    plot_df = pd.concat(frames, ignore_index=True)
    if plot_df.empty:
        return

    interval_values = [str(x) for x in plot_df[interval_col].cat.categories] if hasattr(plot_df[interval_col], "cat") else sorted(plot_df[interval_col].dropna().astype(str).unique())
    stat_values = list(dict.fromkeys(plot_df["DeltaStatistic"].tolist()))

    groups = []
    positions = []
    labels = []
    pos = 1
    gap = 1
    for interval in interval_values:
        interval_positions = []
        for stat in stat_values:
            vals = plot_df[(plot_df[interval_col].astype(str) == interval) & (plot_df["DeltaStatistic"] == stat)]["DeltaFDR"].dropna().to_numpy()
            if len(vals):
                groups.append(vals)
                positions.append(pos)
                interval_positions.append(pos)
            pos += 1
        if interval_positions:
            labels.append((sum(interval_positions) / len(interval_positions), interval))
        pos += gap
    if not groups:
        return

    fig, ax = plt.subplots(figsize=(max(10, len(positions) * 0.55), 5.5))
    ax.boxplot(groups, positions=positions, showfliers=False)
    ax.axhline(0, linewidth=1)
    ax.set_xticks([x for x, _ in labels])
    ax.set_xticklabels([lab for _, lab in labels], rotation=45, ha="right")
    ax.set_xlabel(interval_col)
    ax.set_ylabel("Delta FDR")
    ax.set_title("Binned-minus-classic FDR deltas by base FDR interval")

    # Add a small text legend rather than relying on colors.
    legend_text = "Statistics per interval: " + ", ".join(stat_values)
    ax.text(0.01, 0.99, legend_text, transform=ax.transAxes, va="top", ha="left", fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def threshold_transition_summary(df: pd.DataFrame, thresholds: List[float], binned_cols: List[str], peptide_col: str) -> pd.DataFrame:
    rows = []
    for t in thresholds:
        base_pass = df["BaseFDR"] <= t
        for bcol in binned_cols:
            if bcol not in df.columns:
                continue
            binned_pass = df[bcol] <= t
            rows.append({
                "Threshold": t,
                "BinnedStatistic": bcol,
                "TotalKeys": int(len(df)),
                "BasePassKeys": int(base_pass.sum()),
                "BinnedPassKeys": int(binned_pass.sum()),
                "PassBothKeys": int((base_pass & binned_pass).sum()),
                "BasePass_BinnedFailKeys": int((base_pass & ~binned_pass).sum()),
                "BaseFail_BinnedPassKeys": int((~base_pass & binned_pass).sum()),
                "UniquePeptides_BasePass": int(df.loc[base_pass, peptide_col].nunique(dropna=True)) if peptide_col in df.columns else np.nan,
                "UniquePeptides_BinnedPass": int(df.loc[binned_pass, peptide_col].nunique(dropna=True)) if peptide_col in df.columns else np.nan,
                "UniquePeptides_BasePass_BinnedFail": int(df.loc[base_pass & ~binned_pass, peptide_col].nunique(dropna=True)) if peptide_col in df.columns else np.nan,
                "UniquePeptides_BaseFail_BinnedPass": int(df.loc[~base_pass & binned_pass, peptide_col].nunique(dropna=True)) if peptide_col in df.columns else np.nan,
            })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Analyze classic-vs-binned FDR deltas by FDR intervals."
    )
    ap.add_argument("--wide_table", required=True, type=Path,
                    help="classic_vs_binned_scan_peptide_fdr_wide.tsv")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--edges", default="0,0.01,0.05,0.10,0.20,0.50,0.75,1.0",
                    help="Comma-separated FDR interval edges. Default: 0,0.01,0.05,0.10,0.20,0.50,0.75,1.0")
    ap.add_argument("--delta_cols", default="DeltaBinnedMedianMinusBase,DeltaBinnedMeanMinusBase",
                    help="Comma-separated delta columns to summarize/plot. Default includes median and mean.")
    ap.add_argument("--binned_cols", default="BinnedFDRMedian,BinnedFDRMean,BinnedFDRMin,BinnedFDRMax",
                    help="Comma-separated binned FDR columns for threshold summaries.")
    ap.add_argument("--thresholds", default="0.01,0.05",
                    help="Comma-separated thresholds for pass/fail transition report.")
    ap.add_argument("--shared_only", action="store_true",
                    help="Analyze only rows with ComparisonStatus == shared.")
    ap.add_argument("--no_plots", action="store_true")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    edges = parse_edges(args.edges)
    delta_cols = [x.strip() for x in args.delta_cols.split(",") if x.strip()]
    binned_cols = [x.strip() for x in args.binned_cols.split(",") if x.strip()]
    thresholds = [float(x.strip()) for x in args.thresholds.split(",") if x.strip()]

    print("Loading wide comparison table...")
    df = pd.read_csv(args.wide_table, sep="\t", dtype=str, low_memory=False)

    if "BaseFDR" not in df.columns:
        raise SystemExit("ERROR: missing required column: BaseFDR")

    if args.shared_only and "ComparisonStatus" in df.columns:
        df = df[df["ComparisonStatus"] == "shared"].copy()

    numeric_cols = ["BaseFDR"] + binned_cols + delta_cols
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Generate expected deltas if the wide table does not already contain them.
    for stat in ["Median", "Mean", "Min", "Max"]:
        dcol = f"DeltaBinned{stat}MinusBase"
        bcol = f"BinnedFDR{stat}"
        if dcol not in df.columns and {bcol, "BaseFDR"}.issubset(df.columns):
            df[dcol] = df[bcol] - df["BaseFDR"]

    usable_delta_cols = [c for c in delta_cols if c in df.columns]
    if not usable_delta_cols:
        raise SystemExit(f"ERROR: none of the requested delta columns exist: {delta_cols}")

    df = df.dropna(subset=["BaseFDR"] + usable_delta_cols, how="any").copy()
    peptide_col = "UnmodifiedPeptide" if "UnmodifiedPeptide" in df.columns else ("CorePeptide" if "CorePeptide" in df.columns else "Peptide")

    print(f"Rows available for interval analysis: {len(df):,}")
    if peptide_col in df.columns:
        print(f"Unique peptide sequences: {df[peptide_col].nunique(dropna=True):,}")
    print()

    add_interval(df, "BaseFDR", edges, "BaseFDRInterval")
    for bcol in binned_cols:
        if bcol in df.columns:
            add_interval(df, bcol, edges, f"{bcol}Interval")

    interval_tables = []
    for dcol in usable_delta_cols:
        interval_tables.append(summarize_delta(df, "BaseFDRInterval", dcol, peptide_col))
        for bcol in binned_cols:
            interval_col = f"{bcol}Interval"
            if interval_col in df.columns:
                interval_tables.append(summarize_delta(df, interval_col, dcol, peptide_col))

    interval_summary = pd.concat(interval_tables, ignore_index=True) if interval_tables else pd.DataFrame()
    interval_summary_out = args.outdir / "fdr_delta_interval_summary.tsv"
    interval_summary.to_csv(interval_summary_out, sep="\t", index=False)

    compact_rows = []
    for dcol in usable_delta_cols:
        tmp = summarize_delta(df, "BaseFDRInterval", dcol, peptide_col)
        tmp.insert(0, "View", f"BaseFDR intervals for {dcol}")
        compact_rows.append(tmp)
    compact = pd.concat(compact_rows, ignore_index=True) if compact_rows else pd.DataFrame()
    compact_out = args.outdir / "fdr_delta_by_base_interval.tsv"
    compact.to_csv(compact_out, sep="\t", index=False)

    transition = threshold_transition_summary(df, thresholds, binned_cols, peptide_col)
    transition_out = args.outdir / "fdr_threshold_transition_summary.tsv"
    transition.to_csv(transition_out, sep="\t", index=False)

    annotated_out = args.outdir / "classic_vs_binned_fdr_wide_with_intervals.tsv"
    df.to_csv(annotated_out, sep="\t", index=False)

    if not args.no_plots:
        plot_dir = args.outdir / "plots"
        plot_dir.mkdir(exist_ok=True)
        for dcol in usable_delta_cols:
            make_histogram(df, dcol, plot_dir / f"hist_{dcol}.png", f"Distribution of {dcol}")
            make_interval_boxplot(
                df,
                "BaseFDRInterval",
                dcol,
                plot_dir / f"boxplot_{dcol}_by_BaseFDRInterval.png",
                f"{dcol} by base FDR interval",
            )
        if len(usable_delta_cols) > 1:
            make_combined_delta_boxplot(
                df,
                "BaseFDRInterval",
                usable_delta_cols,
                plot_dir / "boxplot_all_delta_stats_by_BaseFDRInterval.png",
            )

    print("Wrote interval summary:       ", interval_summary_out)
    print("Wrote base-interval summary:  ", compact_out)
    print("Wrote threshold summary:      ", transition_out)
    print("Wrote interval-annotated rows:", annotated_out)
    if not args.no_plots:
        print("Wrote plots in:               ", args.outdir / "plots")
    print()

    print("Delta by BaseFDR interval")
    print("-" * 80)
    show_cols = [
        "DeltaColumn", "Interval", "N_ScanPeptideKeys", "N_UniquePeptides",
        "DeltaMean", "DeltaMedian", "DeltaStd", "Skew", "ExcessKurtosis",
        "FractionNegative", "NormalityTest", "NormalityPValue",
    ]
    if not compact.empty:
        print(compact[[c for c in show_cols if c in compact.columns]].to_string(index=False))

    print()
    print("Threshold transition summary")
    print("-" * 80)
    if not transition.empty:
        print(transition.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
