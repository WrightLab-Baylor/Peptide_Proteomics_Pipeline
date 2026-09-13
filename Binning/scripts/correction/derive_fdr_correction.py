#!/usr/bin/env python3
"""
derive_fdr_correction.py

Derive a condition-specific empirical FDR correction curve from scan-peptide
pairs observed in both a classical search and one partitioned search condition.

The model is intentionally transparent:

    raw_delta = partitioned_fdr - classical_fdr
    raw_corrected_fdr = clip(
        classical_fdr + linear_interpolation(classical_fdr; raw_delta_anchors),
        0,
        1,
    )
    corrected_fdr = clip(
        linear_interpolation(
            classical_fdr; fdr_anchors, monotonic_corrected_fdr_anchors
        ),
        0,
        1,
    )

Key design choices
------------------
- One formula is derived for one tested partition condition.
- ScanNum + peptide is the default stable input observation key.
- Rolling windows operate on distinct classical FDR values, preventing tied
  FDR blocks from being split across neighboring windows.
- The full 0-1 FDR spectrum is modeled by default.
- Quantile bins are built independently inside user-defined FDR segments, so
  the dense high-FDR tail cannot reduce resolution in the 0-0.05 or 0-0.2
  regions.
- Sparse adjacent bins are merged within each segment until they meet the
  requested minimum support when possible.
- Median delta is the default robust statistic; mean delta is optional.
- Explicit anchors at FDR=0 and FDR=1 eliminate out-of-domain endpoint
  extrapolation for valid FDR values.
- All shared rows are used by default; an optional deterministic hash split remains available for diagnostics.

Expected input
--------------
A TSV containing, at minimum, the default columns:
    ScanNum
    CorePeptide
    BaseFDR
    BinnedFDRMedian

The optional default delta column is:
    DeltaBinnedMedianMinusBase

If the delta column is absent, it is calculated as:
    BinnedFDRMedian - BaseFDR

All column names can be overridden with the corresponding *_col arguments.

Primary outputs
---------------
correction_formula.json
correction_anchors.tsv
correction_interval_summary.tsv
formula_metadata.tsv
validation_overall.tsv
validation_by_fdr_interval.tsv
cutoff_behavior.tsv
prediction_sample.tsv

This script derives and validates the formula. It does not modify pipeline FDR
files; application is handled by a separate script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

SCRIPT_VERSION = "1.4.0"
DEFAULT_SEGMENTS = (
    "0,0.001,5;"
    "0.001,0.01,10;"
    "0.01,0.05,20;"
    "0.05,0.1,10;"
    "0.1,0.2,10;"
    "0.2,0.5,15;"
    "0.5,1.0,20"
)
DEFAULT_VALIDATION_WINDOWS = (
    "0,0.001;"
    "0.001,0.005;"
    "0.005,0.01;"
    "0.01,0.02;"
    "0.02,0.05;"
    "0.05,0.1;"
    "0.1,0.2;"
    "0.2,0.5;"
    "0.5,1.0;"
    "0,0.05;"
    "0,0.1;"
    "0,0.2;"
    "0,1.0"
)
DEFAULT_CUTOFFS = "0.01,0.05,0.1,0.2"


@dataclass(frozen=True)
class Segment:
    low: float
    high: float
    requested_bins: int
    index: int

    @property
    def segment_id(self) -> str:
        return f"segment_{self.index:02d}_{self.low:g}_{self.high:g}"


@dataclass(frozen=True)
class Window:
    low: float
    high: float
    index: int

    @property
    def label(self) -> str:
        return f"{self.low:g}-{self.high:g}"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Derive a condition-specific, full-spectrum linear-interpolation "
            "FDR correction curve from shared classical/partitioned scan-peptide pairs."
        )
    )
    p.add_argument("--shared", required=True, type=Path,
                   help="Shared classical-versus-partitioned scan-peptide transition TSV.")
    p.add_argument("--outdir", required=True, type=Path,
                   help="Output directory for the formula package and validation files.")
    p.add_argument("--condition_name", required=True,
                   help="Stable name for this tested condition, e.g. MyDataset_pairwise.")
    p.add_argument("--statistic", choices=["median", "mean"], default="median",
                   help="Local delta statistic used as interpolation anchors. Default: median.")
    p.add_argument("--segments", default=DEFAULT_SEGMENTS,
                   help=(
                       "Semicolon-separated low,high,n_quantile_bins segments. "
                       "Segments are fitted independently. Default covers 0-1 with "
                       "extra resolution below 0.2."
                   ))
    p.add_argument("--validation_windows", default=DEFAULT_VALIDATION_WINDOWS,
                   help="Semicolon-separated low,high windows for validation summaries.")
    p.add_argument("--cutoffs", default=DEFAULT_CUTOFFS,
                   help="Comma-separated FDR cutoffs for pass/fail diagnostics.")
    p.add_argument("--curve_method", choices=["disjoint_bins", "rolling_median"],
                   default="rolling_median",
                   help=("Curve-construction method. rolling_median builds dense overlapping "
                         "local windows; disjoint_bins preserves the original quantile-bin method. "
                         "Default: rolling_median."))
    p.add_argument("--min_bin_n", type=int, default=100,
                   help="Minimum training rows per final disjoint interval. Used only by --curve_method disjoint_bins. Default: 100.")
    p.add_argument("--rolling_window_n", type=int, default=100,
                   help=("Distinct classical FDR values per overlapping rolling window. "
                         "Tied scan-peptide rows are summarized before windowing. Default: 100."))
    p.add_argument("--rolling_stride_n", type=int, default=10,
                   help=("Distinct classical FDR values advanced between rolling windows. "
                         "Default: 10."))
    p.add_argument("--no_merge_sparse_bins", action="store_true",
                   help="Do not merge adjacent sparse bins within a segment.")
    p.add_argument("--monotonic_method", choices=["isotonic", "none"], default="isotonic",
                   help=("Constraint applied to corrected-FDR anchors. isotonic uses weighted "
                         "pool-adjacent-violators regression to enforce a nondecreasing production "
                         "mapping while retaining raw rolling anchors for diagnostics. Default: isotonic."))
    p.add_argument("--test_fraction", type=float, default=0.0,
                   help=("Optional deterministic held-out validation fraction. "
                         "Use 0 to fit the production curve with all rows. Default: 0."))
    p.add_argument("--sample_rows", type=int, default=2000,
                   help="Maximum prediction rows written for inspection. Default: 2000.")
    p.add_argument("--random_seed", type=int, default=1,
                   help="Seed used only for prediction-sample selection. Default: 1.")

    p.add_argument("--scan_col", default="ScanNum")
    p.add_argument("--peptide_col", default="CorePeptide")
    p.add_argument("--classic_fdr_col", default="BaseFDR")
    p.add_argument("--partitioned_fdr_col", default="BinnedFDRMedian")
    p.add_argument("--delta_col", default="DeltaBinnedMedianMinusBase")
    return p.parse_args()


def parse_segments(spec: str) -> list[Segment]:
    segments: list[Segment] = []
    for index, raw in enumerate(spec.split(";"), start=1):
        raw = raw.strip()
        if not raw:
            continue
        parts = [x.strip() for x in raw.split(",")]
        if len(parts) != 3:
            raise ValueError(f"Invalid segment '{raw}'. Expected low,high,n_bins.")
        low, high = float(parts[0]), float(parts[1])
        n_bins = int(parts[2])
        if not (0.0 <= low < high <= 1.0):
            raise ValueError(f"Invalid segment '{raw}'. Require 0 <= low < high <= 1.")
        if n_bins < 1:
            raise ValueError(f"Invalid segment '{raw}'. n_bins must be >= 1.")
        segments.append(Segment(low, high, n_bins, index))

    if not segments:
        raise ValueError("No segments were supplied.")
    segments.sort(key=lambda s: (s.low, s.high))

    tol = 1e-12
    if abs(segments[0].low - 0.0) > tol or abs(segments[-1].high - 1.0) > tol:
        raise ValueError("Segments must cover the full FDR domain from 0 to 1.")
    for left, right in zip(segments[:-1], segments[1:]):
        if abs(left.high - right.low) > tol:
            raise ValueError(
                f"Segments must be contiguous and non-overlapping: "
                f"{left.low:g}-{left.high:g} followed by {right.low:g}-{right.high:g}."
            )
    return segments


def parse_windows(spec: str) -> list[Window]:
    windows: list[Window] = []
    for index, raw in enumerate(spec.split(";"), start=1):
        raw = raw.strip()
        if not raw:
            continue
        parts = [x.strip() for x in raw.split(",")]
        if len(parts) != 2:
            raise ValueError(f"Invalid validation window '{raw}'. Expected low,high.")
        low, high = float(parts[0]), float(parts[1])
        if not (0.0 <= low < high <= 1.0):
            raise ValueError(f"Invalid validation window '{raw}'.")
        windows.append(Window(low, high, index))
    return windows


def parse_cutoffs(spec: str) -> list[float]:
    cutoffs = sorted({float(x.strip()) for x in spec.split(",") if x.strip()})
    if not cutoffs:
        raise ValueError("At least one cutoff is required.")
    if any(x < 0 or x > 1 for x in cutoffs):
        raise ValueError("All cutoffs must be between 0 and 1.")
    return cutoffs


def stable_test_flag(scan: object, peptide: object, test_fraction: float) -> bool:
    key = f"{scan}|{peptide}".encode("utf-8", errors="replace")
    digest = hashlib.sha256(key).hexdigest()
    value = int(digest[:16], 16) / float(0xFFFFFFFFFFFFFFFF)
    return value < test_fraction


def interval_mask(values: pd.Series, low: float, high: float, include_high: bool) -> pd.Series:
    if include_high:
        return (values >= low) & (values <= high)
    return (values >= low) & (values < high)


def validate_columns(path: Path, columns: Sequence[str]) -> None:
    available = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    missing = [c for c in columns if c not in available]
    if missing:
        raise ValueError(
            f"Missing required columns in {path}: {missing}. "
            f"Available columns include: {available[:80]}"
        )


def load_input(args: argparse.Namespace) -> tuple[pd.DataFrame, bool]:
    required = [
        args.scan_col,
        args.peptide_col,
        args.classic_fdr_col,
        args.partitioned_fdr_col,
    ]
    validate_columns(args.shared, required)

    header = list(pd.read_csv(args.shared, sep="\t", nrows=0).columns)
    delta_present = args.delta_col in header
    usecols = set(required)
    if delta_present:
        usecols.add(args.delta_col)

    df = pd.read_csv(args.shared, sep="\t", usecols=lambda c: c in usecols, low_memory=False)
    df = df.rename(columns={
        args.scan_col: "ScanNum",
        args.peptide_col: "PeptideKey",
        args.classic_fdr_col: "ClassicFDR",
        args.partitioned_fdr_col: "PartitionedFDR",
    })
    if delta_present:
        df = df.rename(columns={args.delta_col: "DeltaFDR"})

    df["ScanNum"] = pd.to_numeric(df["ScanNum"], errors="coerce")
    df["ClassicFDR"] = pd.to_numeric(df["ClassicFDR"], errors="coerce")
    df["PartitionedFDR"] = pd.to_numeric(df["PartitionedFDR"], errors="coerce")
    if delta_present:
        df["DeltaFDR"] = pd.to_numeric(df["DeltaFDR"], errors="coerce")

    df["PeptideKey"] = df["PeptideKey"].astype("string")
    df = df.dropna(subset=["ScanNum", "PeptideKey", "ClassicFDR", "PartitionedFDR"]).copy()
    df = df[df["PeptideKey"].str.len().fillna(0) > 0].copy()
    df["ScanNum"] = df["ScanNum"].astype("int64")
    df = df[
        df["ClassicFDR"].between(0.0, 1.0, inclusive="both")
        & df["PartitionedFDR"].between(0.0, 1.0, inclusive="both")
    ].copy()

    calculated_delta = df["PartitionedFDR"] - df["ClassicFDR"]
    if delta_present:
        bad_delta = df["DeltaFDR"].isna() | ~np.isclose(
            df["DeltaFDR"].to_numpy(dtype=float),
            calculated_delta.to_numpy(dtype=float),
            rtol=1e-9,
            atol=1e-12,
        )
        if bad_delta.any():
            # The two FDR columns are authoritative; normalize inconsistent delta values.
            df.loc[bad_delta, "DeltaFDR"] = calculated_delta.loc[bad_delta]
    else:
        df["DeltaFDR"] = calculated_delta

    # Enforce one observation per canonical key. Exact duplicates are harmless;
    # conflicting duplicates are collapsed conservatively using best FDR values.
    df = (
        df.groupby(["ScanNum", "PeptideKey"], as_index=False, dropna=False)
        .agg(
            ClassicFDR=("ClassicFDR", "min"),
            PartitionedFDR=("PartitionedFDR", "min"),
            SourceRowCount=("ClassicFDR", "size"),
        )
    )
    df["DeltaFDR"] = df["PartitionedFDR"] - df["ClassicFDR"]
    return df, delta_present


def build_initial_bins(
    train: pd.DataFrame,
    segments: Sequence[Segment],
    min_bin_n: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    bin_counter = 0
    values = train["ClassicFDR"]

    for seg_index, segment in enumerate(segments):
        include_high = seg_index == len(segments) - 1
        sub = values.loc[interval_mask(values, segment.low, segment.high, include_high)].dropna()
        n_rows = len(sub)
        if n_rows == 0:
            # Preserve full-domain coverage even when a condition has no rows in a segment.
            actual_bins = 1
            edges = np.asarray([segment.low, segment.high], dtype=float)
        else:
            max_supported_bins = max(1, n_rows // min_bin_n) if min_bin_n > 0 else segment.requested_bins
            actual_bins = max(1, min(segment.requested_bins, max_supported_bins))
            if actual_bins == 1:
                edges = np.asarray([segment.low, segment.high], dtype=float)
            else:
                quantiles = np.linspace(0.0, 1.0, actual_bins + 1)
                edges = np.quantile(sub.to_numpy(dtype=float), quantiles)
                edges[0] = segment.low
                edges[-1] = segment.high
                edges = np.unique(edges)
                if len(edges) < 2:
                    edges = np.asarray([segment.low, segment.high], dtype=float)
                actual_bins = len(edges) - 1

        for j in range(actual_bins):
            bin_counter += 1
            rows.append({
                "bin_id": f"bin_{bin_counter:03d}",
                "segment_id": segment.segment_id,
                "segment_index": segment.index,
                "segment_low": segment.low,
                "segment_high": segment.high,
                "requested_bins_in_segment": segment.requested_bins,
                "initial_bins_in_segment": actual_bins,
                "classic_fdr_min_edge": float(edges[j]),
                "classic_fdr_max_edge": float(edges[j + 1]),
            })
    return pd.DataFrame(rows)


def assign_bin_ids(values: pd.Series, bins: pd.DataFrame) -> pd.Series:
    assigned = pd.Series(pd.NA, index=values.index, dtype="object")
    for i, row in bins.reset_index(drop=True).iterrows():
        include_high = i == len(bins) - 1
        mask = interval_mask(
            values,
            float(row["classic_fdr_min_edge"]),
            float(row["classic_fdr_max_edge"]),
            include_high,
        )
        assigned.loc[mask] = row["bin_id"]
    return assigned


def add_bin_counts(train: pd.DataFrame, bins: pd.DataFrame) -> pd.DataFrame:
    out = bins.copy().reset_index(drop=True)
    assigned = assign_bin_ids(train["ClassicFDR"], out)
    counts = assigned.value_counts(dropna=False)
    out["initial_n_train"] = out["bin_id"].map(counts).fillna(0).astype(int)
    return out


def merge_sparse_bins(train: pd.DataFrame, bins: pd.DataFrame, min_bin_n: int) -> pd.DataFrame:
    work = add_bin_counts(train, bins)
    if min_bin_n <= 0:
        work["merged_from_initial_bin_count"] = 1
        work["merged_from_initial_bin_ids"] = work["bin_id"]
        work["post_merge_enabled"] = False
        work["segment_underpowered_for_min_bin_n"] = False
        work["meets_min_bin_n_after_merge"] = True
        return work

    merged_rows: list[dict[str, object]] = []
    final_counter = 0

    for segment_id, seg in work.groupby("segment_id", sort=False):
        seg = seg.sort_values("classic_fdr_min_edge").reset_index(drop=True)
        segment_total = int(seg["initial_n_train"].sum())
        underpowered = segment_total < min_bin_n
        pending: list[pd.Series] = []
        pending_n = 0

        def flush_pending() -> None:
            nonlocal pending, pending_n, final_counter
            if not pending:
                return
            first = pending[0]
            last = pending[-1]
            final_counter += 1
            merged_rows.append({
                "bin_id": f"bin_{final_counter:03d}",
                "segment_id": segment_id,
                "segment_index": int(first["segment_index"]),
                "segment_low": float(first["segment_low"]),
                "segment_high": float(first["segment_high"]),
                "requested_bins_in_segment": int(first["requested_bins_in_segment"]),
                "initial_bins_in_segment": int(first["initial_bins_in_segment"]),
                "classic_fdr_min_edge": float(first["classic_fdr_min_edge"]),
                "classic_fdr_max_edge": float(last["classic_fdr_max_edge"]),
                "initial_n_train": int(pending_n),
                "segment_total_n_train": segment_total,
                "segment_underpowered_for_min_bin_n": bool(underpowered),
                "merged_from_initial_bin_count": len(pending),
                "merged_from_initial_bin_ids": ";".join(str(x["bin_id"]) for x in pending),
                "post_merge_enabled": True,
                "meets_min_bin_n_after_merge": bool(pending_n >= min_bin_n or underpowered),
            })
            pending = []
            pending_n = 0

        for _, row in seg.iterrows():
            pending.append(row)
            pending_n += int(row["initial_n_train"])
            if pending_n >= min_bin_n:
                flush_pending()

        if pending:
            # Attach a trailing sparse remainder to the previous bin in the same
            # segment where possible; otherwise retain it as an underpowered bin.
            if merged_rows and merged_rows[-1]["segment_id"] == segment_id and not underpowered:
                prev = merged_rows[-1]
                prev["classic_fdr_max_edge"] = float(pending[-1]["classic_fdr_max_edge"])
                prev["initial_n_train"] = int(prev["initial_n_train"]) + pending_n
                prev["merged_from_initial_bin_count"] = int(prev["merged_from_initial_bin_count"]) + len(pending)
                prev["merged_from_initial_bin_ids"] = (
                    str(prev["merged_from_initial_bin_ids"])
                    + ";"
                    + ";".join(str(x["bin_id"]) for x in pending)
                )
                prev["meets_min_bin_n_after_merge"] = int(prev["initial_n_train"]) >= min_bin_n
                pending = []
                pending_n = 0
            else:
                flush_pending()

    merged = pd.DataFrame(merged_rows)
    if not merged.empty:
        merged["final_bins_in_segment"] = merged.groupby("segment_id")["bin_id"].transform("count")
    return merged


def summarize_intervals(train: pd.DataFrame, bins: pd.DataFrame) -> pd.DataFrame:
    assigned = assign_bin_ids(train["ClassicFDR"], bins)
    work = train.copy()
    work["bin_id"] = assigned
    if work["bin_id"].isna().any():
        raise RuntimeError("Some training rows were not assigned to a full-spectrum interval.")

    stats = (
        work.groupby("bin_id", observed=True)
        .agg(
            n_train=("ClassicFDR", "size"),
            classic_fdr_min_observed=("ClassicFDR", "min"),
            classic_fdr_max_observed=("ClassicFDR", "max"),
            median_classic_fdr=("ClassicFDR", "median"),
            mean_classic_fdr=("ClassicFDR", "mean"),
            median_partitioned_fdr=("PartitionedFDR", "median"),
            mean_partitioned_fdr=("PartitionedFDR", "mean"),
            median_delta=("DeltaFDR", "median"),
            mean_delta=("DeltaFDR", "mean"),
            sd_delta=("DeltaFDR", "std"),
            q10_delta=("DeltaFDR", lambda x: x.quantile(0.10)),
            q25_delta=("DeltaFDR", lambda x: x.quantile(0.25)),
            q75_delta=("DeltaFDR", lambda x: x.quantile(0.75)),
            q90_delta=("DeltaFDR", lambda x: x.quantile(0.90)),
            fraction_partitioned_lower=("DeltaFDR", lambda x: float((x < 0).mean())),
            fraction_partitioned_equal=("DeltaFDR", lambda x: float((x == 0).mean())),
            fraction_partitioned_higher=("DeltaFDR", lambda x: float((x > 0).mean())),
        )
        .reset_index()
    )
    out = bins.merge(stats, on="bin_id", how="left")
    out["sem_delta"] = out["sd_delta"] / np.sqrt(out["n_train"])
    return out



def summarize_distinct_fdr_levels(train: pd.DataFrame) -> pd.DataFrame:
    """Collapse training rows to one equally weighted observation per ClassicFDR.

    Exact repeated classical FDR values are treated as one level. This prevents
    a tied FDR block from being split between adjacent rolling windows. Median
    and mean outcomes are retained so ``--statistic`` can select the production
    summary later without rebuilding the level table.
    """
    levels = (
        train.groupby("ClassicFDR", as_index=False, sort=True)
        .agg(
            source_row_count=("ClassicFDR", "size"),
            source_scan_count=("ScanNum", "nunique"),
            source_peptide_count=("PeptideKey", "nunique"),
            level_median_partitioned_fdr=("PartitionedFDR", "median"),
            level_mean_partitioned_fdr=("PartitionedFDR", "mean"),
            level_median_delta=("DeltaFDR", "median"),
            level_mean_delta=("DeltaFDR", "mean"),
            level_sd_delta=("DeltaFDR", "std"),
        )
    )
    levels["level_min_partitioned_fdr"] = (
        train.groupby("ClassicFDR", sort=True)["PartitionedFDR"].min().to_numpy()
    )
    levels["level_max_partitioned_fdr"] = (
        train.groupby("ClassicFDR", sort=True)["PartitionedFDR"].max().to_numpy()
    )
    return levels


def summarize_rolling_windows(
    train: pd.DataFrame,
    segments: Sequence[Segment],
    window_n: int,
    stride_n: int,
) -> pd.DataFrame:
    """Build overlapping summaries across distinct classical FDR values.

    Before windowing, all scan-peptide rows sharing an exact ``ClassicFDR`` are
    collapsed to one FDR-level observation. Each window therefore spans up to
    ``window_n`` distinct FDR values and advances by ``stride_n`` distinct
    values. No repeated FDR block can be divided between neighboring windows.
    Each distinct FDR level receives equal weight in the window statistic.
    """
    rows: list[dict[str, object]] = []
    window_counter = 0
    levels = summarize_distinct_fdr_levels(train)

    for seg_i, segment in enumerate(segments):
        include_high = seg_i == len(segments) - 1
        sub = levels.loc[
            interval_mask(levels["ClassicFDR"], segment.low, segment.high, include_high)
        ].sort_values("ClassicFDR").reset_index(drop=True)
        n_levels = len(sub)
        if n_levels == 0:
            continue

        effective_window = min(window_n, n_levels)
        if n_levels <= effective_window:
            starts = [0]
        else:
            starts = list(range(0, n_levels - effective_window + 1, stride_n))
            final_start = n_levels - effective_window
            if starts[-1] != final_start:
                starts.append(final_start)

        for local_index, start in enumerate(starts, start=1):
            stop = start + effective_window
            w = sub.iloc[start:stop]
            window_counter += 1

            # Window statistics are calculated across one summary per distinct
            # classical FDR value, rather than across duplicated source rows.
            median_delta_values = w["level_median_delta"]
            mean_delta_values = w["level_mean_delta"]
            source_rows = int(w["source_row_count"].sum())

            rows.append({
                "bin_id": f"window_{window_counter:05d}",
                "segment_id": segment.segment_id,
                "segment_index": segment.index,
                "segment_low": segment.low,
                "segment_high": segment.high,
                "requested_bins_in_segment": segment.requested_bins,
                "initial_bins_in_segment": pd.NA,
                "classic_fdr_min_edge": float(w["ClassicFDR"].min()),
                "classic_fdr_max_edge": float(w["ClassicFDR"].max()),
                "initial_n_train": int(len(w)),
                "segment_total_n_train": int(n_levels),
                "segment_total_source_rows": int(sub["source_row_count"].sum()),
                "segment_underpowered_for_min_bin_n": bool(n_levels < window_n),
                "merged_from_initial_bin_count": pd.NA,
                "merged_from_initial_bin_ids": "",
                "post_merge_enabled": False,
                "meets_min_bin_n_after_merge": bool(len(w) >= window_n or n_levels < window_n),
                "final_bins_in_segment": int(len(starts)),
                "rolling_window_index": local_index,
                "rolling_start_distinct_fdr_index": int(start),
                "rolling_stop_distinct_fdr_index_exclusive": int(stop),
                "rolling_window_unit": "distinct_classical_fdr",
                "rolling_window_n_requested": int(window_n),
                "rolling_window_n_effective": int(len(w)),
                "rolling_stride_n": int(stride_n),
                "n_train": int(len(w)),
                "n_source_rows": source_rows,
                "n_source_scans": int(w["source_scan_count"].sum()),
                "n_source_peptides_summed_across_levels": int(w["source_peptide_count"].sum()),
                "classic_fdr_min_observed": float(w["ClassicFDR"].min()),
                "classic_fdr_max_observed": float(w["ClassicFDR"].max()),
                "median_classic_fdr": float(w["ClassicFDR"].median()),
                "mean_classic_fdr": float(w["ClassicFDR"].mean()),
                "median_partitioned_fdr": float(w["level_median_partitioned_fdr"].median()),
                "mean_partitioned_fdr": float(w["level_mean_partitioned_fdr"].mean()),
                "median_delta": float(median_delta_values.median()),
                "mean_delta": float(mean_delta_values.mean()),
                "sd_delta": float(median_delta_values.std()) if len(w) > 1 else np.nan,
                "q10_delta": float(median_delta_values.quantile(0.10)),
                "q25_delta": float(median_delta_values.quantile(0.25)),
                "q75_delta": float(median_delta_values.quantile(0.75)),
                "q90_delta": float(median_delta_values.quantile(0.90)),
                "fraction_partitioned_lower": float((median_delta_values < 0).mean()),
                "fraction_partitioned_equal": float((median_delta_values == 0).mean()),
                "fraction_partitioned_higher": float((median_delta_values > 0).mean()),
            })

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No populated distinct-FDR rolling windows could be derived.")
    out["sem_delta"] = out["sd_delta"] / np.sqrt(out["n_train"])
    return out

def build_anchors(intervals: pd.DataFrame, statistic: str) -> pd.DataFrame:
    delta_col = "median_delta" if statistic == "median" else "mean_delta"
    anchors = intervals.loc[
        intervals["n_train"].fillna(0) > 0,
        [
            "bin_id",
            "segment_id",
            "classic_fdr_min_edge",
            "classic_fdr_max_edge",
            "n_train",
            "median_classic_fdr",
            "mean_classic_fdr",
            "median_delta",
            "mean_delta",
        ],
    ].copy()
    anchors["delta_anchor"] = anchors[delta_col]
    anchors["fdr_anchor"] = anchors["median_classic_fdr"].astype(float)
    anchors["delta_anchor"] = pd.to_numeric(anchors["delta_anchor"], errors="coerce")
    anchors = anchors.dropna(subset=["fdr_anchor", "delta_anchor"]).copy()
    if anchors.empty:
        raise ValueError("No supported formula anchors could be derived.")

    anchors = anchors.sort_values(["fdr_anchor", "classic_fdr_min_edge"]).reset_index(drop=True)
    if anchors["fdr_anchor"].duplicated().any():
        anchors = (
            anchors.groupby("fdr_anchor", as_index=False)
            .agg(
                delta_anchor=("delta_anchor", "median" if statistic == "median" else "mean"),
                n_train=("n_train", "sum"),
                bin_id=("bin_id", lambda x: ";".join(map(str, x))),
                segment_id=("segment_id", lambda x: ";".join(sorted(set(map(str, x))))),
                classic_fdr_min_edge=("classic_fdr_min_edge", "min"),
                classic_fdr_max_edge=("classic_fdr_max_edge", "max"),
                median_classic_fdr=("median_classic_fdr", "median"),
                mean_classic_fdr=("mean_classic_fdr", "mean"),
                median_delta=("median_delta", "median"),
                mean_delta=("mean_delta", "mean"),
            )
        )

    # Add exact domain boundaries so every valid FDR in [0,1] is interpolated,
    # never extrapolated. Boundary deltas inherit the nearest observed interval.
    boundary_rows: list[dict[str, object]] = []
    if float(anchors["fdr_anchor"].iloc[0]) > 0.0:
        first = anchors.iloc[0]
        boundary_rows.append({
            **first.to_dict(),
            "bin_id": "boundary_0",
            "segment_id": "domain_boundary",
            "fdr_anchor": 0.0,
            "classic_fdr_min_edge": 0.0,
            "classic_fdr_max_edge": 0.0,
            "n_train": 0,
        })
    if float(anchors["fdr_anchor"].iloc[-1]) < 1.0:
        last = anchors.iloc[-1]
        boundary_rows.append({
            **last.to_dict(),
            "bin_id": "boundary_1",
            "segment_id": "domain_boundary",
            "fdr_anchor": 1.0,
            "classic_fdr_min_edge": 1.0,
            "classic_fdr_max_edge": 1.0,
            "n_train": 0,
        })
    if boundary_rows:
        anchors = pd.concat([anchors, pd.DataFrame(boundary_rows)], ignore_index=True)

    anchors = anchors.sort_values("fdr_anchor").reset_index(drop=True)
    anchors["anchor_index"] = np.arange(len(anchors), dtype=int)
    anchors["statistic"] = statistic
    anchors["anchor_role"] = np.where(
        anchors["bin_id"].astype(str).str.startswith("boundary_"),
        "domain_boundary",
        "observed_interval",
    )
    return anchors[
        [
            "anchor_index",
            "anchor_role",
            "statistic",
            "fdr_anchor",
            "delta_anchor",
            "bin_id",
            "segment_id",
            "classic_fdr_min_edge",
            "classic_fdr_max_edge",
            "n_train",
            "median_classic_fdr",
            "mean_classic_fdr",
            "median_delta",
            "mean_delta",
        ]
    ]



def weighted_isotonic_nondecreasing(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted pool-adjacent-violators algorithm for a nondecreasing fit."""
    y = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    if y.ndim != 1 or w.ndim != 1 or len(y) != len(w):
        raise ValueError("Isotonic values and weights must be equal-length one-dimensional arrays.")
    if len(y) == 0:
        return y.copy()
    if np.any(~np.isfinite(y)) or np.any(~np.isfinite(w)) or np.any(w <= 0):
        raise ValueError("Isotonic values must be finite and weights must be finite and positive.")

    levels: list[float] = []
    block_weights: list[float] = []
    starts: list[int] = []
    ends: list[int] = []

    for i, (yi, wi) in enumerate(zip(y, w)):
        levels.append(float(yi))
        block_weights.append(float(wi))
        starts.append(i)
        ends.append(i)
        while len(levels) >= 2 and levels[-2] > levels[-1]:
            new_weight = block_weights[-2] + block_weights[-1]
            new_level = (levels[-2] * block_weights[-2] + levels[-1] * block_weights[-1]) / new_weight
            levels[-2] = new_level
            block_weights[-2] = new_weight
            ends[-2] = ends[-1]
            levels.pop(); block_weights.pop(); starts.pop(); ends.pop()

    fitted = np.empty(len(y), dtype=float)
    for level, start, end in zip(levels, starts, ends):
        fitted[start:end + 1] = level
    return fitted


def apply_monotonic_constraint(anchors: pd.DataFrame, method: str) -> pd.DataFrame:
    """Add raw and production corrected-FDR anchors, optionally enforcing monotonicity."""
    out = anchors.copy()
    out["raw_delta_anchor"] = out["delta_anchor"].astype(float)
    out["raw_corrected_fdr_anchor"] = np.clip(
        out["fdr_anchor"].to_numpy(dtype=float) + out["raw_delta_anchor"].to_numpy(dtype=float),
        0.0,
        1.0,
    )
    if method == "isotonic":
        weights = np.maximum(out["n_train"].to_numpy(dtype=float), 1.0)
        production = weighted_isotonic_nondecreasing(
            out["raw_corrected_fdr_anchor"].to_numpy(dtype=float), weights
        )
    elif method == "none":
        production = out["raw_corrected_fdr_anchor"].to_numpy(dtype=float)
    else:
        raise ValueError(f"Unsupported monotonic method: {method}")

    out["corrected_fdr_anchor"] = np.clip(production, 0.0, 1.0)
    out["delta_anchor"] = out["corrected_fdr_anchor"] - out["fdr_anchor"].to_numpy(dtype=float)
    out["monotonic_adjustment"] = out["corrected_fdr_anchor"] - out["raw_corrected_fdr_anchor"]
    out["monotonic_method"] = method
    return out

def predict(values: Iterable[float], anchors: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    x = anchors["fdr_anchor"].to_numpy(dtype=float)
    y = anchors["corrected_fdr_anchor"].to_numpy(dtype=float)
    vals = np.asarray(list(values) if not isinstance(values, np.ndarray) else values, dtype=float)
    if len(x) < 2:
        raise ValueError("At least two anchors are required for full-domain interpolation.")
    if not (np.isclose(x[0], 0.0) and np.isclose(x[-1], 1.0)):
        raise ValueError("Formula anchors do not span the full 0-1 domain.")

    # Interpolate the monotonic production mapping directly. Clip once more as
    # a defensive guarantee that no numerical edge case can produce an FDR
    # below 0 or above 1.
    corrected = np.clip(np.interp(vals, x, y), 0.0, 1.0)
    delta_hat = corrected - vals
    return delta_hat, corrected


def error_summary(df: pd.DataFrame, dataset_label: str) -> dict[str, object]:
    error = df["PredictedFDR"] - df["PartitionedFDR"]
    abs_error = error.abs()
    return {
        "dataset": dataset_label,
        "n": int(len(df)),
        "mean_error_predicted_minus_partitioned": float(error.mean()) if len(df) else np.nan,
        "median_error_predicted_minus_partitioned": float(error.median()) if len(df) else np.nan,
        "mae": float(abs_error.mean()) if len(df) else np.nan,
        "median_abs_error": float(abs_error.median()) if len(df) else np.nan,
        "q75_abs_error": float(abs_error.quantile(0.75)) if len(df) else np.nan,
        "q90_abs_error": float(abs_error.quantile(0.90)) if len(df) else np.nan,
        "q95_abs_error": float(abs_error.quantile(0.95)) if len(df) else np.nan,
        "max_abs_error": float(abs_error.max()) if len(df) else np.nan,
        "classic_fdr_median": float(df["ClassicFDR"].median()) if len(df) else np.nan,
        "partitioned_fdr_median": float(df["PartitionedFDR"].median()) if len(df) else np.nan,
        "predicted_fdr_median": float(df["PredictedFDR"].median()) if len(df) else np.nan,
    }


def validation_by_windows(
    df: pd.DataFrame,
    dataset_label: str,
    windows: Sequence[Window],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for window in windows:
        include_high = np.isclose(window.high, 1.0)
        sub = df.loc[interval_mask(df["ClassicFDR"], window.low, window.high, include_high)]
        row = error_summary(sub, dataset_label)
        row.update({
            "window_index": window.index,
            "classic_fdr_window": window.label,
            "window_low": window.low,
            "window_high": window.high,
        })
        rows.append(row)
    return pd.DataFrame(rows)


def cutoff_behavior(df: pd.DataFrame, dataset_label: str, cutoffs: Sequence[float]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cutoff in cutoffs:
        classic_pass = df["ClassicFDR"] <= cutoff
        observed_pass = df["PartitionedFDR"] <= cutoff
        predicted_pass = df["PredictedFDR"] <= cutoff
        tp = int((predicted_pass & observed_pass).sum())
        tn = int((~predicted_pass & ~observed_pass).sum())
        fp = int((predicted_pass & ~observed_pass).sum())
        fn = int((~predicted_pass & observed_pass).sum())
        n = len(df)
        rows.append({
            "dataset": dataset_label,
            "cutoff": cutoff,
            "n": n,
            "classic_pass_count": int(classic_pass.sum()),
            "observed_partitioned_pass_count": int(observed_pass.sum()),
            "predicted_pass_count": int(predicted_pass.sum()),
            "true_positive_predicted_and_observed_pass": tp,
            "true_negative_predicted_and_observed_fail": tn,
            "false_positive_predicted_pass_observed_fail": fp,
            "false_negative_predicted_fail_observed_pass": fn,
            "accuracy": (tp + tn) / n if n else np.nan,
            "precision_predicted_pass": tp / (tp + fp) if (tp + fp) else np.nan,
            "recall_observed_partitioned_pass": tp / (tp + fn) if (tp + fn) else np.nan,
            "observed_gain_vs_classic": int((observed_pass & ~classic_pass).sum()),
            "observed_loss_vs_classic": int((~observed_pass & classic_pass).sum()),
            "predicted_gain_vs_classic": int((predicted_pass & ~classic_pass).sum()),
            "predicted_loss_vs_classic": int((~predicted_pass & classic_pass).sum()),
        })
    return pd.DataFrame(rows)


def write_metadata(
    path: Path,
    args: argparse.Namespace,
    df: pd.DataFrame,
    train: pd.DataFrame,
    test: pd.DataFrame,
    intervals: pd.DataFrame,
    anchors: pd.DataFrame,
    delta_present: bool,
) -> None:
    rows = [
        ("script", Path(__file__).name),
        ("script_version", SCRIPT_VERSION),
        ("generated_utc", datetime.now(timezone.utc).isoformat()),
        ("python_version", platform.python_version()),
        ("pandas_version", pd.__version__),
        ("numpy_version", np.__version__),
        ("condition_name", args.condition_name),
        ("source_shared_tsv", str(args.shared.resolve())),
        ("statistic", args.statistic),
        ("curve_method", args.curve_method),
        ("rolling_window_n", args.rolling_window_n),
        ("rolling_stride_n", args.rolling_stride_n),
        ("rolling_window_unit", "distinct_classical_fdr" if args.curve_method == "rolling_median" else "source_rows"),
        ("distinct_classical_fdr_values", df["ClassicFDR"].nunique()),
        ("monotonic_method", args.monotonic_method),
        ("raw_corrected_anchor_reversals", int((anchors["raw_corrected_fdr_anchor"].diff() < -1e-15).sum())),
        ("production_corrected_anchor_reversals", int((anchors["corrected_fdr_anchor"].diff() < -1e-15).sum())),
        ("max_abs_monotonic_adjustment", float(anchors["monotonic_adjustment"].abs().max())),
        ("formula", "corrected_fdr = linear_interp(classic_fdr; fdr_anchor, corrected_fdr_anchor)"),
        ("formula_domain", "0-1 inclusive"),
        ("endpoint_policy", "explicit anchors at 0 and 1; no valid-domain extrapolation"),
        ("segments", args.segments),
        ("min_bin_n", args.min_bin_n),
        ("sparse_bin_merging_enabled", int(not args.no_merge_sparse_bins)),
        ("test_fraction", args.test_fraction),
        ("fit_scope", "all_rows" if args.test_fraction == 0 else "training_subset"),
        ("input_delta_column_present", int(delta_present)),
        ("rows_after_validation_and_key_collapse", len(df)),
        ("train_rows", len(train)),
        ("test_rows", len(test)),
        ("unique_scans", df["ScanNum"].nunique()),
        ("unique_scan_peptide_pairs", len(df)),
        ("intervals_total", len(intervals)),
        ("intervals_with_training_support", int((intervals["n_train"].fillna(0) > 0).sum())),
        ("anchors_total_including_boundaries", len(anchors)),
        ("scan_column", args.scan_col),
        ("peptide_column", args.peptide_col),
        ("classic_fdr_column", args.classic_fdr_col),
        ("partitioned_fdr_column", args.partitioned_fdr_col),
        ("delta_column", args.delta_col),
    ]
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(path, sep="\t", index=False)


def main() -> int:
    args = parse_args()
    if not args.shared.exists():
        raise SystemExit(f"ERROR: shared input not found: {args.shared}")
    if not (0.0 <= args.test_fraction < 1.0):
        raise SystemExit("ERROR: --test_fraction must be at least 0 and less than 1.")
    if args.min_bin_n < 1:
        raise SystemExit("ERROR: --min_bin_n must be at least 1.")
    if args.rolling_window_n < 2:
        raise SystemExit("ERROR: --rolling_window_n must be at least 2.")
    if args.rolling_stride_n < 1:
        raise SystemExit("ERROR: --rolling_stride_n must be at least 1.")

    try:
        segments = parse_segments(args.segments)
        windows = parse_windows(args.validation_windows)
        cutoffs = parse_cutoffs(args.cutoffs)
        df, delta_present = load_input(args)
    except (ValueError, OSError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    if df.empty:
        raise SystemExit("ERROR: no valid shared scan-peptide rows remain after filtering.")

    if args.test_fraction == 0.0:
        df["IsTest"] = False
        train = df.copy()
        test = df.iloc[0:0].copy()
    else:
        df["IsTest"] = [
            stable_test_flag(scan, peptide, args.test_fraction)
            for scan, peptide in zip(df["ScanNum"], df["PeptideKey"])
        ]
        train = df.loc[~df["IsTest"]].copy()
        test = df.loc[df["IsTest"]].copy()
        if train.empty or test.empty:
            raise SystemExit(
                "ERROR: deterministic train/test split produced an empty partition. "
                "Adjust --test_fraction or verify the input size."
            )

    if args.curve_method == "rolling_median":
        intervals = summarize_rolling_windows(
            train=train,
            segments=segments,
            window_n=args.rolling_window_n,
            stride_n=args.rolling_stride_n,
        )
    else:
        initial_bins = build_initial_bins(train, segments, args.min_bin_n)
        if args.no_merge_sparse_bins:
            bins = add_bin_counts(train, initial_bins)
            bins["segment_total_n_train"] = bins.groupby("segment_id")["initial_n_train"].transform("sum")
            bins["segment_underpowered_for_min_bin_n"] = bins["segment_total_n_train"] < args.min_bin_n
            bins["merged_from_initial_bin_count"] = 1
            bins["merged_from_initial_bin_ids"] = bins["bin_id"]
            bins["post_merge_enabled"] = False
            bins["meets_min_bin_n_after_merge"] = bins["initial_n_train"] >= args.min_bin_n
            bins["final_bins_in_segment"] = bins.groupby("segment_id")["bin_id"].transform("count")
        else:
            bins = merge_sparse_bins(train, initial_bins, args.min_bin_n)
        intervals = summarize_intervals(train, bins)

    anchors = build_anchors(intervals, args.statistic)
    anchors = apply_monotonic_constraint(anchors, args.monotonic_method)

    prediction_sets = [train]
    if not test.empty:
        prediction_sets.append(test)
    for subset in prediction_sets:
        delta_hat, predicted = predict(subset["ClassicFDR"].to_numpy(dtype=float), anchors)
        subset["PredictedDelta"] = delta_hat
        subset["PredictedFDR"] = predicted
        subset["PredictionError"] = subset["PredictedFDR"] - subset["PartitionedFDR"]
        subset["AbsolutePredictionError"] = subset["PredictionError"].abs()

    args.outdir.mkdir(parents=True, exist_ok=True)

    intervals.to_csv(args.outdir / "correction_interval_summary.tsv", sep="\t", index=False)
    anchors.to_csv(args.outdir / "correction_anchors.tsv", sep="\t", index=False)

    formula = {
        "schema_version": 1,
        "script_version": SCRIPT_VERSION,
        "condition_name": args.condition_name,
        "model": "monotonic_linear_interpolation_of_local_delta",
        "curve_method": args.curve_method,
        "statistic": args.statistic,
        "formula": "corrected_fdr = interp(classic_fdr; fdr_anchor, corrected_fdr_anchor)",
        "raw_formula": "raw_corrected_fdr = clip(classic_fdr + interp(classic_fdr; fdr_anchor, raw_delta_anchor), 0, 1)",
        "monotonic_method": args.monotonic_method,
        "domain": {"minimum": 0.0, "maximum": 1.0, "inclusive": True},
        "endpoint_policy": "Explicit anchors at 0 and 1; no extrapolation for valid FDR values.",
        "training_key": [args.scan_col, args.peptide_col],
        "source_shared_tsv": str(args.shared.resolve()),
        "source_columns": {
            "scan": args.scan_col,
            "peptide": args.peptide_col,
            "classic_fdr": args.classic_fdr_col,
            "partitioned_fdr": args.partitioned_fdr_col,
            "delta": args.delta_col,
        },
        "segments": [
            {"low": s.low, "high": s.high, "requested_quantile_bins": s.requested_bins}
            for s in segments
        ],
        "min_bin_n": args.min_bin_n,
        "rolling_window_n": args.rolling_window_n,
        "rolling_stride_n": args.rolling_stride_n,
        "rolling_window_unit": "distinct_classical_fdr" if args.curve_method == "rolling_median" else "source_rows",
        "distinct_classical_fdr_values": int(df["ClassicFDR"].nunique()),
        "sparse_bin_merging_enabled": not args.no_merge_sparse_bins,
        "test_fraction": args.test_fraction,
        "fit_scope": "all_rows" if args.test_fraction == 0 else "training_subset",
        "fdr_anchor": anchors["fdr_anchor"].astype(float).tolist(),
        "delta_anchor": anchors["delta_anchor"].astype(float).tolist(),
        "corrected_fdr_anchor": anchors["corrected_fdr_anchor"].astype(float).tolist(),
        "raw_delta_anchor": anchors["raw_delta_anchor"].astype(float).tolist(),
        "raw_corrected_fdr_anchor": anchors["raw_corrected_fdr_anchor"].astype(float).tolist(),
        "monotonic_adjustment": anchors["monotonic_adjustment"].astype(float).tolist(),
        "anchor_role": anchors["anchor_role"].astype(str).tolist(),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
    }
    with (args.outdir / "correction_formula.json").open("w", encoding="utf-8") as handle:
        json.dump(formula, handle, indent=2)

    validation_rows = [error_summary(train, "all_data" if test.empty else "train")]
    if not test.empty:
        validation_rows.append(error_summary(test, "test"))
    validation_overall = pd.DataFrame(validation_rows)
    validation_overall.to_csv(args.outdir / "validation_overall.tsv", sep="\t", index=False)

    validation_frames = [
        validation_by_windows(train, "all_data" if test.empty else "train", windows)
    ]
    if not test.empty:
        validation_frames.append(validation_by_windows(test, "test", windows))
    validation_windows = pd.concat(validation_frames, ignore_index=True)
    validation_windows.to_csv(args.outdir / "validation_by_fdr_interval.tsv", sep="\t", index=False)

    cutoff_frames = [cutoff_behavior(train, "all_data" if test.empty else "train", cutoffs)]
    if not test.empty:
        cutoff_frames.append(cutoff_behavior(test, "test", cutoffs))
    cutoff = pd.concat(cutoff_frames, ignore_index=True)
    cutoff.to_csv(args.outdir / "cutoff_behavior.tsv", sep="\t", index=False)

    combined_frames = [train.assign(DatasetSplit="all_data" if test.empty else "train")]
    if not test.empty:
        combined_frames.append(test.assign(DatasetSplit="test"))
    combined = pd.concat(combined_frames, ignore_index=True)
    sample = combined.sample(
        n=min(args.sample_rows, len(combined)),
        random_state=args.random_seed,
    ).sort_values(["ClassicFDR", "ScanNum", "PeptideKey"])
    sample.to_csv(args.outdir / "prediction_sample.tsv", sep="\t", index=False)

    write_metadata(
        args.outdir / "formula_metadata.tsv",
        args,
        df,
        train,
        test,
        intervals,
        anchors,
        delta_present,
    )

    print("FDR correction formula derived successfully.")
    print(f"Condition:          {args.condition_name}")
    print(f"Statistic:          {args.statistic}")
    print(f"Curve method:       {args.curve_method}")
    print(f"Monotonic method:   {args.monotonic_method}")
    if args.curve_method == "rolling_median":
        print("Rolling unit:       distinct classical FDR values")
    if test.empty:
        print(f"Rows:               {len(df):,} (all rows used for fitting)")
    else:
        print(f"Rows:               {len(df):,} ({len(train):,} train; {len(test):,} test)")
    print(f"Supported intervals:{int((intervals['n_train'].fillna(0) > 0).sum()):,}")
    print(f"Formula anchors:    {len(anchors):,} (including domain boundaries)")
    print(f"Output directory:   {args.outdir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        raise SystemExit(130)
