#!/usr/bin/env python3
"""
compare_classic_vs_binned_fdr_summary.py

Build a transparent classic-vs-binned FDR comparison table using the generic
outputs from collect_pair_bin_fdr_values.py.

Inputs
------
1. Classical/full-run withsyn TSV, usually:
       <classic_run>/fdr_esti/<sample>_withsyn.tsv

2. Core-level binned FDR summary from collect_pair_bin_fdr_values.py:
       <binning_run>/Binning_Downstream/scan_core_peptide_fdr_summary.tsv

The script matches at ScanNum + CorePeptide. By default, modification notation
is retained in CorePeptide, while UnmodifiedPeptide is also carried for peptide-
level reporting.

Outputs
-------
classic_vs_binned_scan_peptide_fdr_wide.tsv
    One row per ScanNum + CorePeptide observed in either source, with explicit
    BaseFDR and BinnedFDR* columns plus delta columns.

classic_vs_binned_threshold_crossing_summary.tsv
    Pass/fail transition summaries for user-specified thresholds.

classic_vs_binned_<status>_scan_peptides.tsv
    Convenience subsets for shared/base_only/binning_only keys.

summary.txt
    Compact human-readable console-style summary.

Design choice
-------------
No generic final 'FDR' column is created. Binned statistics remain explicitly
named: BinnedFDRMean, BinnedFDRMedian, BinnedFDRMin, BinnedFDRMax, etc.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd


def normalize_core_peptide(value: object, strip_mods: bool = False) -> str:
    """Strip PHRP/MS-GF+ flanking residue notation; optionally strip PTM symbols."""
    if pd.isna(value):
        return ""
    s = str(value).strip()
    if not s:
        return ""
    if s.count(".") >= 2:
        first = s.find(".")
        last = s.rfind(".")
        if first < last:
            s = s[first + 1:last]
    s = s.strip()
    if strip_mods:
        s = re.sub(r"[^A-Z]", "", s)
    return s


def boolish_any(series: pd.Series) -> bool:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y", "t"}).any()


def require_columns(df: pd.DataFrame, required: set[str], path: Path) -> None:
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: {path} missing required columns: {sorted(missing)}")


def join_unique(values: Iterable[object], max_items: int = 100) -> str:
    vals = sorted({str(v) for v in values if pd.notna(v) and str(v) != ""})
    if len(vals) > max_items:
        return ";".join(vals[:max_items]) + f";...(+{len(vals) - max_items} more)"
    return ";".join(vals)


def load_base_withsyn(
    path: Path,
    scan_col: str,
    peptide_col: str,
    protein_col: str,
    fdr_col: str,
) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    require_columns(df, {scan_col, peptide_col, fdr_col}, path)

    df = df.copy()
    df[scan_col] = pd.to_numeric(df[scan_col], errors="coerce")
    df = df[df[scan_col].notna()].copy()
    df[scan_col] = df[scan_col].astype("int64")

    df["CorePeptide"] = df[peptide_col].map(lambda x: normalize_core_peptide(x, strip_mods=False))
    df["UnmodifiedPeptide"] = df[peptide_col].map(lambda x: normalize_core_peptide(x, strip_mods=True))
    df["BaseFDR_num"] = pd.to_numeric(df[fdr_col], errors="coerce")
    df = df[df["CorePeptide"] != ""].copy()

    group_cols = [scan_col, "CorePeptide"]
    agg_dict = {
        "BaseFDR": ("BaseFDR_num", "median"),
        "BaseFDRMean": ("BaseFDR_num", "mean"),
        "BaseFDRMin": ("BaseFDR_num", "min"),
        "BaseFDRMax": ("BaseFDR_num", "max"),
        "BaseRows": (fdr_col, "size"),
        "BasePeptideForms": (peptide_col, join_unique),
        "UnmodifiedPeptide": ("UnmodifiedPeptide", lambda x: sorted(set(x))[0] if len(set(x)) else ""),
    }
    if protein_col in df.columns:
        agg_dict.update({
            "BaseProteinCount": (protein_col, "nunique"),
            "BaseProteins": (protein_col, join_unique),
        })
    if "IsDecoy" in df.columns:
        agg_dict["BaseIsDecoy"] = ("IsDecoy", boolish_any)

    out = df.groupby(group_cols, dropna=False).agg(**agg_dict).reset_index()
    out = out.rename(columns={scan_col: "ScanNum"})
    if "BaseIsDecoy" not in out.columns:
        out["BaseIsDecoy"] = False
    if "BaseProteinCount" not in out.columns:
        out["BaseProteinCount"] = np.nan
        out["BaseProteins"] = ""
    return out


def load_binned_summary(path: Path, scan_col: str, peptide_col: str) -> pd.DataFrame:
    """
    Load an already-finalized one-row-per-ScanNum+CorePeptide binned summary.

    The collector owns all pair-level aggregation. This function validates
    uniqueness and never performs a second median-of-medians collapse.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    required = {scan_col, "CorePeptide", "BinnedFDRMean", "BinnedFDRMedian"}
    require_columns(df, required, path)

    df = df.copy()
    df[scan_col] = pd.to_numeric(df[scan_col], errors="coerce")
    df = df[df[scan_col].notna()].copy()
    df[scan_col] = df[scan_col].astype("int64")
    df["CorePeptide"] = df["CorePeptide"].map(
        lambda x: normalize_core_peptide(x, strip_mods=False)
    )
    df = df[df["CorePeptide"] != ""].copy()

    if "UnmodifiedPeptide" in df.columns:
        df["UnmodifiedPeptide_binned"] = df["UnmodifiedPeptide"].astype(str)
    else:
        df["UnmodifiedPeptide_binned"] = df["CorePeptide"].map(
            lambda x: normalize_core_peptide(x, strip_mods=True)
        )

    numeric_cols = [
        "BinnedFDRMean", "BinnedFDRMedian", "BinnedFDRMin", "BinnedFDRMax",
        "BinnedFDRStd", "BinnedFDRRange", "BinnedFDRCount", "RowsAcrossBins",
        "PairBinsRepresented", "MSGFScoreMin", "MSGFScoreMax", "MSGFScoreRange",
        "BinnedPeptideFormCount", "BinnedProteinCount",
        "MaxExactFormsWithinPair", "PairsWithMultipleExactForms",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    key_cols = [scan_col, "CorePeptide"]
    duplicate_mask = df.duplicated(key_cols, keep=False)
    if duplicate_mask.any():
        preview_cols = [c for c in key_cols + ["BinnedFDRMean", "BinnedFDRMedian"] if c in df.columns]
        preview = df.loc[duplicate_mask, preview_cols].head(20).to_string(index=False)
        raise SystemExit(
            "ERROR: binned core summary is not unique by ScanNum + CorePeptide. "
            "The comparison script no longer performs a second aggregation. "
            f"Example duplicate rows:\n{preview}"
        )

    if "BinnedPeptideForms" not in df.columns:
        if peptide_col in df.columns:
            df["BinnedPeptideForms"] = df[peptide_col].astype(str)
        else:
            df["BinnedPeptideForms"] = df["CorePeptide"]

    df = df.rename(columns={scan_col: "ScanNum"})
    keep = [
        "ScanNum", "CorePeptide", "UnmodifiedPeptide_binned",
        "BinnedPeptideForms", "BinnedPeptideFormCount",
        "BinnedFDRMean", "BinnedFDRMedian", "BinnedFDRMin", "BinnedFDRMax",
        "BinnedFDRStd", "BinnedFDRRange", "BinnedFDRCount",
        "RowsAcrossBins", "PairBinsRepresented", "BinnedRunFolders",
        "MSGFScoreMin", "MSGFScoreMax", "MSGFScoreRange",
        "BinnedProteinCount", "BinnedProteins",
        "MaxExactFormsWithinPair", "PairsWithMultipleExactForms",
    ]
    return df[[c for c in keep if c in df.columns]].copy()


def classify_status(base_present: bool, bin_present: bool) -> str:
    if base_present and bin_present:
        return "shared"
    if base_present:
        return "base_only"
    if bin_present:
        return "binning_only"
    return "neither"


def add_delta_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for stat in ["Mean", "Median", "Min", "Max"]:
        bcol = f"BinnedFDR{stat}"
        dcol = f"DeltaBinned{stat}MinusBase"
        if bcol in out.columns:
            out[dcol] = out[bcol] - out["BaseFDR"]
    return out


def add_threshold_columns(df: pd.DataFrame, thresholds: List[float]) -> pd.DataFrame:
    out = df.copy()
    for t in thresholds:
        tag = str(t).replace(".", "p")
        out[f"BasePass_{tag}"] = out["BaseFDR"].le(t)
        for stat in ["Mean", "Median", "Min", "Max"]:
            bcol = f"BinnedFDR{stat}"
            if bcol in out.columns:
                out[f"Binned{stat}Pass_{tag}"] = out[bcol].le(t)
                out[f"Delta{stat}MinusBase_{tag}"] = out[bcol] - out["BaseFDR"]
    return out


def threshold_report(df: pd.DataFrame, thresholds: List[float]) -> pd.DataFrame:
    rows = []
    for t in thresholds:
        base_pass = df["BaseFDR"].le(t).fillna(False)
        shared = df["ComparisonStatus"].eq("shared")
        for stat_col, label in [
            ("BinnedFDRMean", "mean"),
            ("BinnedFDRMedian", "median"),
            ("BinnedFDRMin", "min"),
            ("BinnedFDRMax", "max"),
        ]:
            if stat_col not in df.columns:
                continue
            bin_pass = df[stat_col].le(t).fillna(False)
            for group_name, mask in [
                ("All ScanNum/CorePeptide keys", pd.Series(True, index=df.index)),
                ("Shared keys only", shared),
            ]:
                sub = df[mask].copy()
                bp = base_pass[mask]
                sp = bin_pass[mask]
                rows.append({
                    "Threshold": t,
                    "BinnedStatistic": label,
                    "Scope": group_name,
                    "TotalKeys": int(len(sub)),
                    "BasePassKeys": int(bp.sum()),
                    "BinnedPassKeys": int(sp.sum()),
                    "PassBothKeys": int((bp & sp).sum()),
                    "BasePass_BinnedFailKeys": int((bp & ~sp).sum()),
                    "BaseFail_BinnedPassKeys": int((~bp & sp).sum()),
                    "UniquePeptides_Total": int(sub["UnmodifiedPeptide"].nunique(dropna=True)),
                    "UniquePeptides_BasePass": int(sub.loc[bp, "UnmodifiedPeptide"].nunique(dropna=True)),
                    "UniquePeptides_BinnedPass": int(sub.loc[sp, "UnmodifiedPeptide"].nunique(dropna=True)),
                    "UniquePeptides_PassBoth": int(sub.loc[bp & sp, "UnmodifiedPeptide"].nunique(dropna=True)),
                    "UniquePeptides_BasePass_BinnedFail": int(sub.loc[bp & ~sp, "UnmodifiedPeptide"].nunique(dropna=True)),
                    "UniquePeptides_BaseFail_BinnedPass": int(sub.loc[~bp & sp, "UnmodifiedPeptide"].nunique(dropna=True)),
                })
    return pd.DataFrame(rows)


def write_text_summary(path: Path, final: pd.DataFrame, report: pd.DataFrame) -> None:
    lines = []
    lines.append("Comparison status")
    lines.append("-" * 80)
    lines.append(final["ComparisonStatus"].value_counts().to_string())
    lines.append("")
    lines.append(f"Total ScanNum/CorePeptide keys: {len(final):,}")
    lines.append(f"Unique unmodified peptides:     {final['UnmodifiedPeptide'].nunique(dropna=True):,}")
    lines.append(f"Shared keys:                    {(final['ComparisonStatus'] == 'shared').sum():,}")
    lines.append("")
    lines.append("FDR delta summary for shared keys")
    lines.append("-" * 80)
    shared = final[final["ComparisonStatus"] == "shared"].copy()
    if shared.empty:
        lines.append("No shared keys.")
    else:
        cols = [
            "BaseFDR", "BinnedFDRMean", "BinnedFDRMedian", "BinnedFDRMin", "BinnedFDRMax",
            "DeltaBinnedMeanMinusBase", "DeltaBinnedMedianMinusBase",
        ]
        lines.append(shared[[c for c in cols if c in shared.columns]].describe().to_string())
    lines.append("")
    lines.append("Threshold crossing report")
    lines.append("-" * 80)
    lines.append(report.to_string(index=False))
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Compare classical full-run FDR with binned FDR summary outputs."
    )
    ap.add_argument("--base_withsyn", required=True, type=Path,
                    help="Classical/full-run withsyn TSV, e.g. <classic>/fdr_esti/<sample>_withsyn.tsv")
    ap.add_argument("--binned_summary", required=True, type=Path,
                    help="Core-level binned summary TSV: Binning_Downstream/scan_core_peptide_fdr_summary.tsv")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--scan_col", default="ScanNum")
    ap.add_argument("--peptide_col", default="Peptide")
    ap.add_argument("--protein_col", default="Protein")
    ap.add_argument("--base_fdr_col", default="FDR")
    ap.add_argument("--thresholds", default="0.01,0.05",
                    help="Comma-separated FDR thresholds for pass reports. Default: 0.01,0.05")
    args = ap.parse_args()

    base_withsyn = args.base_withsyn.expanduser().resolve()
    binned_summary_path = args.binned_summary.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    thresholds = [float(x.strip()) for x in args.thresholds.split(",") if x.strip()]

    print("Loading classical/base withsyn...")
    base = load_base_withsyn(base_withsyn, args.scan_col, args.peptide_col, args.protein_col, args.base_fdr_col)
    print(f"Base ScanNum/CorePeptide keys: {len(base):,}")

    print("Loading binned FDR summary...")
    binned = load_binned_summary(binned_summary_path, args.scan_col, args.peptide_col)
    print(f"Binned ScanNum/CorePeptide keys: {len(binned):,}")

    final = base.merge(binned, on=["ScanNum", "CorePeptide"], how="outer", suffixes=("", "_bin"))

    if "UnmodifiedPeptide_binned" in final.columns:
        final["UnmodifiedPeptide"] = final["UnmodifiedPeptide"].fillna(final["UnmodifiedPeptide_binned"])
        final = final.drop(columns=["UnmodifiedPeptide_binned"], errors="ignore")

    final["InBase"] = final["BaseFDR"].notna()
    final["InBinning"] = final["BinnedFDRMedian"].notna() if "BinnedFDRMedian" in final.columns else False
    final["ComparisonStatus"] = [classify_status(b, n) for b, n in zip(final["InBase"], final["InBinning"])]
    final = add_delta_columns(final)
    final = add_threshold_columns(final, thresholds)

    leading = [
        "ScanNum", "CorePeptide", "UnmodifiedPeptide", "ComparisonStatus", "InBase", "InBinning",
        "BaseFDR", "BaseFDRMean", "BaseFDRMin", "BaseFDRMax",
        "BinnedFDRMean", "BinnedFDRMedian", "BinnedFDRMin", "BinnedFDRMax", "BinnedFDRStd",
        "BinnedFDRRange", "BinnedFDRCount", "PairBinsRepresented", "RowsAcrossBins",
        "DeltaBinnedMeanMinusBase", "DeltaBinnedMedianMinusBase", "DeltaBinnedMinMinusBase", "DeltaBinnedMaxMinusBase",
        "BaseRows", "BaseProteinCount", "BinnedProteinCount", "BaseIsDecoy",
        "BasePeptideForms", "BinnedPeptideForms", "BaseProteins", "BinnedProteins",
        "BinnedRunFolders", "BinnedPairLabels",
    ]
    ordered = [c for c in leading if c in final.columns] + [c for c in final.columns if c not in leading]
    final = final[ordered].sort_values(["ComparisonStatus", "ScanNum", "CorePeptide"], kind="mergesort")

    main_out = outdir / "classic_vs_binned_scan_peptide_fdr_wide.tsv"
    final.to_csv(main_out, sep="\t", index=False)

    report = threshold_report(final, thresholds)
    report_out = outdir / "classic_vs_binned_threshold_crossing_summary.tsv"
    report.to_csv(report_out, sep="\t", index=False)

    for status in ["shared", "base_only", "binning_only"]:
        sub = final[final["ComparisonStatus"] == status].copy()
        sub.to_csv(outdir / f"classic_vs_binned_{status}_scan_peptides.tsv", sep="\t", index=False)

    summary_out = outdir / "summary.txt"
    write_text_summary(summary_out, final, report)

    print(f"\nWrote main wide table:       {main_out}")
    print(f"Wrote threshold summary:     {report_out}")
    print(f"Wrote status-specific files: {outdir}/classic_vs_binned_<status>_scan_peptides.tsv")
    print(f"Wrote summary:               {summary_out}")
    print()
    print(summary_out.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
