#!/usr/bin/env python3
"""
collate_binning_results.py

Collate scan-recovery, peptide/protein gain, and optional FDR summary
statistics from Binning runs into a self-contained results folder for each
numbered run.

Expected layout
---------------
<root>/
├── Fecal_1/
│   ├── Fecal_Binning/
│   │   ├── Group_Downstream/
│   │   └── Pairwise_Downstream/
│   ├── Fecal_Groups_1.2M/
│   ├── Fecal_Groups_600k/
│   └── Fecal_Groups_<LABEL>/   # arbitrary additional size-series labels supported
├── Fecal_2/ ...
├── Ocean_1/ ...
├── Soil_1/ ...
├── Fecal_Full/
├── Ocean_Full/
├── Soil_Full/
└── results/

Run discovery is generic: any top-level directory named <condition>_<integer>
with a matching <condition>_Full baseline is eligible. Therefore a user can run
one sample (e.g. MyDataset_1) or many without editing this script.

For each numbered run, writes
------------------------------
results/<run>/
├── binning_results_summary.tsv
├── main_partition_scan_recovery_summary.tsv
├── group_size_scan_recovery_summary.tsv
├── main_partition_peptide_protein_gains.tsv
├── group_size_peptide_protein_gains.tsv
├── input_manifest.tsv
├── results_guide.tsv
└── recovery_inputs/
    ├── main_individual_recovery_summary.tsv
    ├── main_pairwise_recovery_summary.tsv
    ├── group_1.2M_recovery_summary.tsv
    ├── group_600k_recovery_summary.tsv
    ├── group_300k_recovery_summary.tsv
    ├── group_150k_recovery_summary.tsv
    └── group_75k_recovery_summary.tsv

The compact TSV outputs are intended for downstream inspection, plotting, and
reporting. When ``Matched_Formula_Output`` is present, the collator also writes
a consolidated matched-formula peptide/protein gain table at the results root.

Notes
-----
- Recovery summaries are validated before use and copied into results/ because
  they are very small and useful for provenance.
- Peptide counts use distinct ``Peptide`` sequences, with distinct
  ``PeptideFlanked`` counts retained as a parallel QC metric. Protein counts use
  distinct ``Protein`` identifiers. This supersedes the legacy row-count logic,
  which over-counted peptides expanded across multiple protein mappings.
- Missing/incomplete analyses warn and are skipped; other numbered runs continue.
- Existing collated TSVs are rewritten completely when the collator is rerun;
  copied recovery-summary files are overwritten with their current sources.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import sys

import numpy as np
import pandas as pd
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

SCRIPT_VERSION = "1.4.4"

PREFERRED_SIZE_LABEL_ORDER = ["1.2M", "600k", "300k", "150k", "75k"]
PARTITION_KIND_INDIVIDUAL = "individual"
PARTITION_KIND_PAIRWISE = "pairwise"

RECOVERY_GROUP_REL = Path(
    "Group_Downstream/classical_guided_culling/"
    "classical_guided_group_reference_recovery_summary.tsv"
)
RECOVERY_PAIR_REL = Path(
    "Pairwise_Downstream/classical_guided_culling/"
    "classical_guided_pair_reference_recovery_summary.tsv"
)
GROUP_PEPTIDE_REL = Path("Group_Downstream/protein_rollup/peptide_crosstab.tsv")
GROUP_PROTEIN_REL = Path("Group_Downstream/protein_rollup/protein_crosstab_annotated.tsv")
PAIR_PEPTIDE_REL = Path("Pairwise_Downstream/protein_rollup/peptide_crosstab.tsv")
PAIR_PROTEIN_REL = Path("Pairwise_Downstream/protein_rollup/protein_crosstab_annotated.tsv")

PEPTIDE_CROSSTAB_NAME = "peptide_crosstab_annotated.tsv"
PROTEIN_CROSSTAB_NAME = "protein_crosstab_annotated.tsv"
PEPTIDE_ID_COLUMN = "Peptide"
PEPTIDE_FLANKED_COLUMN = "PeptideFlanked"
PROTEIN_ID_COLUMN = "Protein"
MATCHED_OUTPUT_DIRNAME = "Matched_Formula_Output"

RUN_RE = re.compile(r"^(?P<condition>.+)_(?P<replicate>[1-9]\d*)$")

FDR_BINNED_STATISTIC = "Median"
FDR_INTERVAL_EDGES = [0.0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.2, 0.5, 1.0]
FDR_THRESHOLDS = [0.01, 0.05]
FDR_SHARED_ONLY = True


@dataclass(frozen=True)
class NumberedRun:
    name: str
    condition: str
    replicate: int
    path: Path
    full_path: Path


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def natural_key(text: str):
    parts = re.split(r"(\d+(?:\.\d+)?)", text)
    key = []
    for part in parts:
        if not part:
            continue
        try:
            key.append((0, float(part)))
        except ValueError:
            key.append((1, part.lower()))
    return key


def size_label_sort_key(label: str) -> tuple[int, object]:
    """Sort familiar paper labels first, then arbitrary labels naturally."""
    if label in PREFERRED_SIZE_LABEL_ORDER:
        return (0, PREFERRED_SIZE_LABEL_ORDER.index(label))
    return (1, natural_key(label))


def discover_size_series(run: NumberedRun) -> list[tuple[str, Path]]:
    """Discover <Condition>_Groups_<LABEL> datasets for one numbered run."""
    prefix = f"{run.condition}_Groups_"
    discovered: list[tuple[str, Path]] = []
    if not run.path.is_dir():
        return discovered
    for candidate in run.path.iterdir():
        if not candidate.is_dir() or not candidate.name.startswith(prefix):
            continue
        label = candidate.name[len(prefix):]
        if not label:
            continue
        discovered.append((label, candidate))
    return sorted(discovered, key=lambda item: size_label_sort_key(item[0]))


def project_relative(path: Path, project_root: Path) -> str:
    """Return a project-root-relative path for portable provenance output."""
    try:
        return str(path.resolve().relative_to(project_root.resolve()))
    except ValueError:
        return str(path.resolve())


def count_partition_directories(condition_root: Path, partition_kind: str) -> int:
    """Count actual Individual-bin or Pairwise-bin directories."""
    if not condition_root.is_dir():
        raise FileNotFoundError(condition_root)

    if partition_kind == PARTITION_KIND_INDIVIDUAL:
        pattern = re.compile(r".*_group_[1-9]\d*$")
    elif partition_kind == PARTITION_KIND_PAIRWISE:
        pattern = re.compile(r".*_pair_[1-9]\d*_[1-9]\d*$")
    else:
        raise ValueError(f"unknown partition kind: {partition_kind}")

    partitions = [
        p for p in condition_root.iterdir()
        if p.is_dir() and pattern.fullmatch(p.name)
    ]
    if not partitions:
        raise ValueError(
            f"{condition_root}: no {partition_kind} partition directories discovered"
        )
    return len(partitions)


def partition_structure(count: int, partition_kind: str) -> str:
    """Format a human-readable partition structure label."""
    if partition_kind == PARTITION_KIND_INDIVIDUAL:
        noun = "group" if count == 1 else "groups"
    elif partition_kind == PARTITION_KIND_PAIRWISE:
        noun = "pair" if count == 1 else "pairs"
    else:
        raise ValueError(f"unknown partition kind: {partition_kind}")
    return f"{count} {noun}"


def discover_numbered_runs(root: Path) -> list[NumberedRun]:
    runs: list[NumberedRun] = []
    for path in root.iterdir():
        if not path.is_dir() or path.name in {"results", "HPC_1"}:
            continue
        match = RUN_RE.match(path.name)
        if not match:
            continue
        condition = match.group("condition")
        replicate = int(match.group("replicate"))
        full_path = root / f"{condition}_Full"
        if not full_path.is_dir():
            warn(
                f"{path.name}: matching full-database baseline directory not found: "
                f"{full_path.name}; skipping"
            )
            continue
        runs.append(NumberedRun(path.name, condition, replicate, path, full_path))

    return sorted(runs, key=lambda r: (natural_key(r.condition), r.replicate))


def _distinct_nonempty(series: pd.Series) -> int:
    """Count distinct non-empty identifiers after conservative string cleanup."""
    values = series.astype("string").str.strip()
    values = values[values.notna() & values.ne("") & values.str.lower().ne("nan")]
    return int(values.nunique(dropna=True))


def count_peptide_entities(path: Path) -> dict[str, int]:
    """Count distinct peptide sequence and flanked-peptide identifiers.

    Annotated peptide crosstabs can contain many rows for one peptide because a
    sequence may map to multiple proteins. Counting rows therefore inflates the
    biological peptide total. Distinct ``Peptide`` is the production metric;
    distinct ``PeptideFlanked`` is retained as a QC quantity.
    """
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    header = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()
    required = [PEPTIDE_ID_COLUMN, PEPTIDE_FLANKED_COLUMN]
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(f"{path} is missing required peptide identifier columns: {missing}")
    df = pd.read_csv(path, sep="\t", usecols=required, low_memory=False)
    return {
        "Peptides": _distinct_nonempty(df[PEPTIDE_ID_COLUMN]),
        "PeptideFlanked": _distinct_nonempty(df[PEPTIDE_FLANKED_COLUMN]),
        "SourceRows": int(len(df)),
    }


def count_protein_entities(path: Path) -> dict[str, int]:
    """Count distinct protein identifiers in a protein crosstab."""
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    header = pd.read_csv(path, sep="\t", nrows=0).columns.astype(str).tolist()
    if PROTEIN_ID_COLUMN not in header:
        raise ValueError(f"{path} is missing required column {PROTEIN_ID_COLUMN!r}")
    df = pd.read_csv(path, sep="\t", usecols=[PROTEIN_ID_COLUMN], low_memory=False)
    return {
        "Proteins": _distinct_nonempty(df[PROTEIN_ID_COLUMN]),
        "SourceRows": int(len(df)),
    }

def read_metric_value_tsv(path: Path) -> tuple[dict[str, str], int]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"Metric", "Value"}.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must contain Metric and Value columns; "
                f"found {reader.fieldnames or []}"
            )
        metrics: dict[str, str] = {}
        rows = 0
        for row in reader:
            rows += 1
            metric = (row.get("Metric") or "").strip()
            value = (row.get("Value") or "").strip()
            if not metric:
                continue
            if metric in metrics:
                raise ValueError(f"duplicate metric {metric!r} in {path}")
            metrics[metric] = value
    return metrics, rows


def first_metric(metrics: dict[str, str], candidates: Iterable[str], path: Path) -> tuple[str, str]:
    for name in candidates:
        if name in metrics:
            return name, metrics[name]
    raise ValueError(f"{path} is missing required metric; tried: {', '.join(candidates)}")


def parse_int_metric(raw: str, name: str, path: Path) -> int:
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"{path}: metric {name} is not numeric: {raw!r}") from exc
    if not math.isfinite(value) or not value.is_integer():
        raise ValueError(f"{path}: metric {name} is not an integer count: {raw!r}")
    return int(value)


def load_recovery(path: Path, partition_kind: str) -> tuple[dict[str, object], dict[str, object]]:
    metrics, metric_rows = read_metric_value_tsv(path)

    total_name, total_raw = first_metric(
        metrics, ["TotalReferenceScans", "ReferenceScans", "TotalScans"], path
    )
    if partition_kind == PARTITION_KIND_INDIVIDUAL:
        recovered_names = [
            "RecoveredInAtLeastOneGroup",
            "RecoveredInAtLeastOnePartition",
            "RecoveredScans",
        ]
        unrecovered_names = [
            "UnrecoveredInAllGroups",
            "Unrecovered",
            "UnrecoveredScans",
        ]
    else:
        recovered_names = [
            "RecoveredInAtLeastOnePair",
            "RecoveredInAtLeastOnePartition",
            "RecoveredScans",
        ]
        unrecovered_names = [
            "UnrecoveredInAllPairs",
            "Unrecovered",
            "UnrecoveredScans",
        ]

    recovered_name, recovered_raw = first_metric(metrics, recovered_names, path)
    unrecovered_name, unrecovered_raw = first_metric(metrics, unrecovered_names, path)
    fraction_name, fraction_raw = first_metric(
        metrics, ["RecoveredFraction", "RecoveryFraction"], path
    )

    total = parse_int_metric(total_raw, total_name, path)
    recovered = parse_int_metric(recovered_raw, recovered_name, path)
    unrecovered = parse_int_metric(unrecovered_raw, unrecovered_name, path)
    try:
        recovered_fraction = float(fraction_raw)
    except ValueError as exc:
        raise ValueError(
            f"{path}: metric {fraction_name} is not numeric: {fraction_raw!r}"
        ) from exc

    if min(total, recovered, unrecovered) < 0:
        raise ValueError(f"{path}: negative recovery count")
    if recovered > total or unrecovered > total:
        raise ValueError(f"{path}: recovery count exceeds reference total")
    if total - recovered != unrecovered:
        raise ValueError(
            f"{path}: total - recovered = {total - recovered}, "
            f"stored unrecovered = {unrecovered}"
        )
    calculated = recovered / total if total else float("nan")
    if total and not math.isclose(
        recovered_fraction, calculated, rel_tol=1e-9, abs_tol=1e-12
    ):
        raise ValueError(
            f"{path}: stored RecoveredFraction={recovered_fraction:.12g}, "
            f"calculated={calculated:.12g}"
        )

    row = {
        "TotalReferenceScans": total,
        "RecoveredScans": recovered,
        "UnrecoveredScans": unrecovered,
        "RecoveredFraction": recovered_fraction,
        "UnrecoveredFraction": unrecovered / total if total else "",
        "RecoveryPercent": recovered_fraction * 100.0,
        "LossPercent": unrecovered / total * 100.0 if total else "",
    }
    manifest = {
        "InputPath": str(path),
        "TotalMetric": total_name,
        "RecoveredMetric": recovered_name,
        "UnrecoveredMetric": unrecovered_name,
        "FractionMetric": fraction_name,
        "MetricRows": metric_rows,
        "ValidationStatus": "PASS",
    }
    return row, manifest


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def baseline_counts_from_full(
    full_path: Path,
    project_root: Path,
) -> tuple[dict[str, int], dict[str, str]]:
    peptide = full_path / PEPTIDE_CROSSTAB_NAME
    protein = full_path / PROTEIN_CROSSTAB_NAME
    pep = count_peptide_entities(peptide)
    pro = count_protein_entities(protein)
    counts = {
        "BaselinePeptides": pep["Peptides"],
        "BaselinePeptideFlanked": pep["PeptideFlanked"],
        "BaselineProteins": pro["Proteins"],
        "BaselinePeptideSourceRows": pep["SourceRows"],
        "BaselineProteinSourceRows": pro["SourceRows"],
    }
    paths = {
        "BaselinePeptidePath": project_relative(peptide, project_root),
        "BaselineProteinPath": project_relative(protein, project_root),
    }
    return counts, paths


def baseline_counts(
    run: NumberedRun,
    project_root: Path,
) -> tuple[dict[str, int], dict[str, str]]:
    return baseline_counts_from_full(run.full_path, project_root)


def gain_row(
    run: NumberedRun,
    downstream_root: Path,
    partition_structure_label: str,
    baseline: dict[str, int],
    baseline_paths: dict[str, str],
    project_root: Path,
    bin_size: str | None = None,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    peptide_path = downstream_root / "protein_rollup" / PEPTIDE_CROSSTAB_NAME
    protein_path = downstream_root / "protein_rollup" / PROTEIN_CROSSTAB_NAME

    pep = count_peptide_entities(peptide_path)
    pro = count_protein_entities(protein_path)
    corrected_peptides = pep["Peptides"]
    corrected_flanked = pep["PeptideFlanked"]
    corrected_proteins = pro["Proteins"]

    row: dict[str, object] = {
        "Dataset": run.condition,
        "Run": run.name,
        "PartitionStructure": partition_structure_label,
        "DownstreamRoot": project_relative(downstream_root, project_root),
        "CorrectedPeptides": corrected_peptides,
        "CorrectedPeptideFlanked": corrected_flanked,
        "CorrectedProteins": corrected_proteins,
        "BaselinePeptides": baseline["BaselinePeptides"],
        "BaselinePeptideFlanked": baseline["BaselinePeptideFlanked"],
        "BaselineProteins": baseline["BaselineProteins"],
        "PeptideGain": corrected_peptides - baseline["BaselinePeptides"],
        "PeptideFlankedGain": corrected_flanked - baseline["BaselinePeptideFlanked"],
        "ProteinGain": corrected_proteins - baseline["BaselineProteins"],
        "PeptideSourceRows": pep["SourceRows"],
        "ProteinSourceRows": pro["SourceRows"],
        "PeptidePath": project_relative(peptide_path, project_root),
        "ProteinPath": project_relative(protein_path, project_root),
        "BaselinePeptidePath": baseline_paths["BaselinePeptidePath"],
        "BaselineProteinPath": baseline_paths["BaselineProteinPath"],
    }
    if bin_size is not None:
        row["BinSize"] = bin_size

    analysis_label = (
        partition_structure_label.replace(" ", "_")
        if bin_size is None
        else f"size_{bin_size}_{partition_structure_label.replace(' ', '_')}"
    )
    manifest_rows = [
        {
            "Run": run.name,
            "Category": "gain_input",
            "Analysis": analysis_label,
            "SourcePath": project_relative(peptide_path, project_root),
            "CopiedPath": "",
            "ValidationStatus": "PASS",
            "Notes": (
                f"source rows={pep['SourceRows']}; distinct Peptide={corrected_peptides}; "
                f"distinct PeptideFlanked={corrected_flanked}"
            ),
        },
        {
            "Run": run.name,
            "Category": "gain_input",
            "Analysis": analysis_label,
            "SourcePath": project_relative(protein_path, project_root),
            "CopiedPath": "",
            "ValidationStatus": "PASS",
            "Notes": f"source rows={pro['SourceRows']}; distinct Protein={corrected_proteins}",
        },
    ]
    return row, manifest_rows



def fdr_quantile(series: pd.Series, q: float) -> float:
    return float(series.quantile(q)) if len(series) else float("nan")


def load_fdr_rows(
    path: Path,
    dataset: str,
    partition_structure_label: str,
    bin_size: str | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)

    binned_col = f"BinnedFDR{FDR_BINNED_STATISTIC}"
    delta_col = f"DeltaBinned{FDR_BINNED_STATISTIC}MinusBase"

    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    required = {"BaseFDR", binned_col, "ComparisonStatus"}
    missing = required - set(header)
    if missing:
        raise ValueError(f"{path} is missing required FDR columns: {sorted(missing)}")

    usecols = ["BaseFDR", binned_col, "ComparisonStatus"]
    if delta_col in header:
        usecols.append(delta_col)

    peptide_col = None
    if "UnmodifiedPeptide" in header:
        peptide_col = "UnmodifiedPeptide"
        usecols.append(peptide_col)
    elif "CorePeptide" in header:
        peptide_col = "CorePeptide"
        usecols.append(peptide_col)

    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    input_rows = len(df)

    if FDR_SHARED_ONLY:
        df = df[df["ComparisonStatus"].astype(str).eq("shared")].copy()

    df["BaseFDR"] = pd.to_numeric(df["BaseFDR"], errors="coerce")
    df[binned_col] = pd.to_numeric(df[binned_col], errors="coerce")
    if delta_col in df.columns:
        df["DeltaFDR"] = pd.to_numeric(df[delta_col], errors="coerce")
    else:
        df["DeltaFDR"] = df[binned_col] - df["BaseFDR"]

    df = df[
        np.isfinite(df["BaseFDR"])
        & np.isfinite(df[binned_col])
        & np.isfinite(df["DeltaFDR"])
    ].copy()

    df["Dataset"] = dataset
    df["PartitionStructure"] = partition_structure_label
    df["BinnedFDR"] = df[binned_col]
    if bin_size is not None:
        df["BinSize"] = bin_size

    metadata = {
        "InputRows": input_rows,
        "AnalysisRows": len(df),
        "SharedOnly": FDR_SHARED_ONLY,
        "BinnedStatistic": FDR_BINNED_STATISTIC,
        "PeptideColumn": peptide_col or "",
    }
    return df, metadata


def summarize_fdr_overall(df: pd.DataFrame) -> dict[str, object]:
    delta = df["DeltaFDR"]
    improved = delta < 0
    unchanged = np.isclose(delta, 0.0, atol=1e-15)
    worsened = delta > 0

    row: dict[str, object] = {
        "Dataset": df["Dataset"].iloc[0],
        "PartitionStructure": df["PartitionStructure"].iloc[0],
        "N": len(df),
        "BaseFDRMean": float(df["BaseFDR"].mean()),
        "BaseFDRMedian": float(df["BaseFDR"].median()),
        "BinnedFDRMean": float(df["BinnedFDR"].mean()),
        "BinnedFDRMedian": float(df["BinnedFDR"].median()),
        "DeltaFDRMean": float(delta.mean()),
        "DeltaFDRMedian": float(delta.median()),
        "DeltaFDRQ05": fdr_quantile(delta, 0.05),
        "DeltaFDRQ25": fdr_quantile(delta, 0.25),
        "DeltaFDRQ75": fdr_quantile(delta, 0.75),
        "DeltaFDRQ95": fdr_quantile(delta, 0.95),
        "ImprovedCount": int(improved.sum()),
        "UnchangedCount": int(unchanged.sum()),
        "WorsenedCount": int(worsened.sum()),
        "ImprovedFraction": float(improved.mean()),
        "WorsenedFraction": float(worsened.mean()),
    }
    if "BinSize" in df.columns:
        row["BinSize"] = df["BinSize"].iloc[0]

    peptide_col = (
        "UnmodifiedPeptide"
        if "UnmodifiedPeptide" in df.columns
        else ("CorePeptide" if "CorePeptide" in df.columns else None)
    )
    row["UniquePeptides"] = int(df[peptide_col].nunique()) if peptide_col else np.nan
    return row


def summarize_fdr_intervals(df: pd.DataFrame) -> list[dict[str, object]]:
    categories = pd.cut(
        df["BaseFDR"],
        bins=FDR_INTERVAL_EDGES,
        include_lowest=True,
        right=True,
        duplicates="raise",
    )
    work = df.assign(BaseFDRInterval=categories)
    rows: list[dict[str, object]] = []

    for interval, sub in work.groupby("BaseFDRInterval", observed=False):
        if pd.isna(interval):
            continue
        left = float(interval.left)
        right = float(interval.right)
        if (
            len(FDR_INTERVAL_EDGES) > 1
            and math.isclose(right, FDR_INTERVAL_EDGES[1])
            and left < FDR_INTERVAL_EDGES[0]
        ):
            left = float(FDR_INTERVAL_EDGES[0])

        delta = sub["DeltaFDR"].dropna()
        row: dict[str, object] = {
            "Dataset": df["Dataset"].iloc[0],
            "PartitionStructure": df["PartitionStructure"].iloc[0],
            "BaseFDRInterval": f"[{left:g}, {right:g}]",
            "IntervalLeft": left,
            "IntervalRight": right,
            "IntervalMidpoint": (left + right) / 2.0,
            "N": len(delta),
            "DeltaFDRMean": float(delta.mean()) if len(delta) else np.nan,
            "DeltaFDRMedian": float(delta.median()) if len(delta) else np.nan,
            "DeltaFDRQ25": fdr_quantile(delta, 0.25),
            "DeltaFDRQ75": fdr_quantile(delta, 0.75),
            "ImprovedFraction": float((delta < 0).mean()) if len(delta) else np.nan,
            "WorsenedFraction": float((delta > 0).mean()) if len(delta) else np.nan,
        }
        if "BinSize" in df.columns:
            row["BinSize"] = df["BinSize"].iloc[0]
        rows.append(row)
    return rows


def summarize_fdr_thresholds(df: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    peptide_col = (
        "UnmodifiedPeptide"
        if "UnmodifiedPeptide" in df.columns
        else ("CorePeptide" if "CorePeptide" in df.columns else None)
    )

    for threshold in FDR_THRESHOLDS:
        base_pass = df["BaseFDR"] <= threshold
        binned_pass = df["BinnedFDR"] <= threshold
        gained = ~base_pass & binned_pass
        lost = base_pass & ~binned_pass

        row: dict[str, object] = {
            "Dataset": df["Dataset"].iloc[0],
            "PartitionStructure": df["PartitionStructure"].iloc[0],
            "Threshold": threshold,
            "N": len(df),
            "BasePassCount": int(base_pass.sum()),
            "PartitionPassCount": int(binned_pass.sum()),
            "NetPassChange": int(binned_pass.sum() - base_pass.sum()),
            "GainedPassCount": int(gained.sum()),
            "LostPassCount": int(lost.sum()),
            "PassBothCount": int((base_pass & binned_pass).sum()),
            "GainedUniquePeptides": (
                int(df.loc[gained, peptide_col].nunique()) if peptide_col else np.nan
            ),
            "LostUniquePeptides": (
                int(df.loc[lost, peptide_col].nunique()) if peptide_col else np.nan
            ),
        }
        if "BinSize" in df.columns:
            row["BinSize"] = df["BinSize"].iloc[0]
        rows.append(row)

    return rows


def collate_fdr_for_run(
    run: NumberedRun,
    out_dir: Path,
    project_root: Path,
    strict: bool,
    include_fdr: bool = False,
) -> bool:
    print(f"[{run.name}] Collating FDR summaries (median, shared rows only)")

    main_root = run.path / f"{run.condition}_Binning"
    try:
        group_structure = partition_structure(
            count_partition_directories(main_root, PARTITION_KIND_INDIVIDUAL),
            PARTITION_KIND_INDIVIDUAL,
        )
        pair_structure = partition_structure(
            count_partition_directories(main_root, PARTITION_KIND_PAIRWISE),
            PARTITION_KIND_PAIRWISE,
        )
    except (OSError, ValueError, FileNotFoundError) as exc:
        warn(f"{run.name}: cannot determine primary FDR partition structure ({exc})")
        return False if strict else True

    main_specs = [
        (
            group_structure,
            main_root
            / "Group_Downstream"
            / "classic_vs_group_fdr_audit"
            / "classic_vs_binned_scan_peptide_fdr_wide.tsv",
        ),
        (
            pair_structure,
            main_root
            / "Pairwise_Downstream"
            / "classic_vs_pair_fdr_audit"
            / "classic_vs_binned_scan_peptide_fdr_wide.tsv",
        ),
    ]

    main_overall: list[dict[str, object]] = []
    main_intervals: list[dict[str, object]] = []
    main_thresholds: list[dict[str, object]] = []
    main_manifest: list[dict[str, object]] = []

    for structure, source in main_specs:
        try:
            df, meta = load_fdr_rows(source, run.condition, structure)
            main_overall.append(summarize_fdr_overall(df))
            main_intervals.extend(summarize_fdr_intervals(df))
            main_thresholds.extend(summarize_fdr_thresholds(df))
            main_manifest.append({
                "Dataset": run.condition,
                "Run": run.name,
                "PartitionStructure": structure,
                "InputPath": project_relative(source, project_root),
                **meta,
                "ValidationStatus": "PASS",
            })
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(f"{run.name}: FDR summary unavailable for {structure} ({exc})")
            if strict:
                return False

    size_overall: list[dict[str, object]] = []
    size_intervals: list[dict[str, object]] = []
    size_thresholds: list[dict[str, object]] = []
    size_manifest: list[dict[str, object]] = []

    for size, condition_root in discover_size_series(run):
        try:
            structure = partition_structure(
                count_partition_directories(
                    condition_root, PARTITION_KIND_INDIVIDUAL
                ),
                PARTITION_KIND_INDIVIDUAL,
            )
            source = (
                condition_root
                / "Group_Downstream"
                / "classic_vs_group_fdr_audit"
                / "classic_vs_binned_scan_peptide_fdr_wide.tsv"
            )
            df, meta = load_fdr_rows(
                source, run.condition, structure, bin_size=size
            )
            size_overall.append(summarize_fdr_overall(df))
            size_intervals.extend(summarize_fdr_intervals(df))
            size_thresholds.extend(summarize_fdr_thresholds(df))
            size_manifest.append({
                "Dataset": run.condition,
                "Run": run.name,
                "BinSize": size,
                "PartitionStructure": structure,
                "InputPath": project_relative(source, project_root),
                **meta,
                "ValidationStatus": "PASS",
            })
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(f"{run.name}: FDR size-series summary unavailable for {size} ({exc})")
            if strict:
                return False

    write_tsv(
        out_dir / "main_partition_fdr_overall_summary.tsv",
        main_overall,
        [
            "Dataset", "PartitionStructure", "N", "UniquePeptides",
            "BaseFDRMean", "BaseFDRMedian", "BinnedFDRMean", "BinnedFDRMedian",
            "DeltaFDRMean", "DeltaFDRMedian", "DeltaFDRQ05", "DeltaFDRQ25",
            "DeltaFDRQ75", "DeltaFDRQ95", "ImprovedCount", "UnchangedCount",
            "WorsenedCount", "ImprovedFraction", "WorsenedFraction",
        ],
    )
    write_tsv(
        out_dir / "main_partition_fdr_interval_summary.tsv",
        main_intervals,
        [
            "Dataset", "PartitionStructure", "BaseFDRInterval", "IntervalLeft",
            "IntervalRight", "IntervalMidpoint", "N", "DeltaFDRMean",
            "DeltaFDRMedian", "DeltaFDRQ25", "DeltaFDRQ75",
            "ImprovedFraction", "WorsenedFraction",
        ],
    )
    write_tsv(
        out_dir / "main_partition_fdr_threshold_summary.tsv",
        main_thresholds,
        [
            "Dataset", "PartitionStructure", "Threshold", "N", "BasePassCount",
            "PartitionPassCount", "NetPassChange", "GainedPassCount",
            "LostPassCount", "PassBothCount", "GainedUniquePeptides",
            "LostUniquePeptides",
        ],
    )
    write_tsv(
        out_dir / "main_partition_fdr_manifest.tsv",
        main_manifest,
        [
            "Dataset", "Run", "PartitionStructure", "InputPath",
            "InputRows", "AnalysisRows", "SharedOnly", "BinnedStatistic",
            "PeptideColumn", "ValidationStatus",
        ],
    )
    write_tsv(
        out_dir / "size_series_fdr_overall_summary.tsv",
        size_overall,
        [
            "Dataset", "BinSize", "PartitionStructure", "N", "UniquePeptides",
            "BaseFDRMean", "BaseFDRMedian", "BinnedFDRMean", "BinnedFDRMedian",
            "DeltaFDRMean", "DeltaFDRMedian", "DeltaFDRQ05", "DeltaFDRQ25",
            "DeltaFDRQ75", "DeltaFDRQ95", "ImprovedCount", "UnchangedCount",
            "WorsenedCount", "ImprovedFraction", "WorsenedFraction",
        ],
    )
    write_tsv(
        out_dir / "size_series_fdr_interval_summary.tsv",
        size_intervals,
        [
            "Dataset", "BinSize", "PartitionStructure", "BaseFDRInterval",
            "IntervalLeft", "IntervalRight", "IntervalMidpoint", "N",
            "DeltaFDRMean", "DeltaFDRMedian", "DeltaFDRQ25", "DeltaFDRQ75",
            "ImprovedFraction", "WorsenedFraction",
        ],
    )
    write_tsv(
        out_dir / "size_series_fdr_threshold_summary.tsv",
        size_thresholds,
        [
            "Dataset", "BinSize", "PartitionStructure", "Threshold", "N",
            "BasePassCount", "PartitionPassCount", "NetPassChange",
            "GainedPassCount", "LostPassCount", "PassBothCount",
            "GainedUniquePeptides", "LostUniquePeptides",
        ],
    )
    write_tsv(
        out_dir / "size_series_fdr_manifest.tsv",
        size_manifest,
        [
            "Dataset", "Run", "BinSize", "PartitionStructure", "InputPath",
            "InputRows", "AnalysisRows", "SharedOnly", "BinnedStatistic",
            "PeptideColumn", "ValidationStatus",
        ],
    )
    return True


def write_results_guide(run: NumberedRun, out_dir: Path) -> None:
    """Write a compact, run-specific guide to the collated result files."""
    rel_run = Path("..") / ".." / run.name
    rel_full = Path("..") / ".." / f"{run.condition}_Full"
    rows = [
        {
            "File": "binning_results_summary.tsv",
            "Meaning": (
                "Combined summary of scan recovery and peptide/protein "
                "gain metrics for the main Individual/Pairwise comparison and all tested "
                "Individual-bin sizes."
            ),
            "Relative_Source_or_Derivation": (
                "Derived from the recovery summaries and peptide/protein crosstabs described "
                "by the other rows in this guide."
            ),
        },
        {
            "File": "main_partition_scan_recovery_summary.tsv",
            "Meaning": (
                "Scan recovery for the primary Individual-bin and Pairwise-bin analyses."
            ),
            "Relative_Source_or_Derivation": (
                f"{rel_run}/{run.condition}_Binning/Group_Downstream/classical_guided_culling/"
                "classical_guided_group_reference_recovery_summary.tsv and "
                f"{rel_run}/{run.condition}_Binning/Pairwise_Downstream/classical_guided_culling/"
                "classical_guided_pair_reference_recovery_summary.tsv"
            ),
        },
        {
            "File": "group_size_scan_recovery_summary.tsv",
            "Meaning": (
                "Scan recovery across dynamically discovered Individual-bin-size-series datasets."
            ),
            "Relative_Source_or_Derivation": (
                f"{rel_run}/{run.condition}_Groups_<LABEL>/Group_Downstream/"
                "classical_guided_culling/classical_guided_group_reference_recovery_summary.tsv"
            ),
        },
        {
            "File": "main_partition_peptide_protein_gains.tsv",
            "Meaning": (
                "Distinct peptide and protein counts for the main Individual and Pairwise analyses, "
                "with gains calculated relative to the condition's full-database baseline."
            ),
            "Relative_Source_or_Derivation": (
                f"Partition counts: {rel_run}/{run.condition}_Binning/"
                "<Group_Downstream|Pairwise_Downstream>/protein_rollup/. "
                f"Baseline counts: {rel_full}/peptide_crosstab_annotated.tsv and "
                f"{rel_full}/protein_crosstab_annotated.tsv."
            ),
        },
        {
            "File": "group_size_peptide_protein_gains.tsv",
            "Meaning": (
                "Distinct peptide and protein counts across the Individual-bin-size series, "
                "with gains calculated relative to the condition's full-database baseline."
            ),
            "Relative_Source_or_Derivation": (
                f"Partition counts: {rel_run}/{run.condition}_Groups_<LABEL>/"
                f"Group_Downstream/protein_rollup/. Baseline counts: "
                f"{rel_full}/peptide_crosstab_annotated.tsv and {rel_full}/protein_crosstab_annotated.tsv."
            ),
        },
        {
            "File": "input_manifest.tsv",
            "Meaning": (
                "Provenance manifest listing the exact source paths used by the collator, "
                "validation status, copied recovery inputs, and row-count notes."
            ),
            "Relative_Source_or_Derivation": "Generated by collate_binning_results.py during collation; stored paths are relative to the project root.",
        },
        {
            "File": "recovery_inputs/*.tsv",
            "Meaning": (
                "Copies of the small original scan-recovery summary files used to build "
                "the collated recovery tables."
            ),
            "Relative_Source_or_Derivation": (
                f"Copied from {rel_run}/<condition>/.../classical_guided_culling/"
                "classical_guided_*_reference_recovery_summary.tsv"
            ),
        },
        {
            "File": "main_partition_fdr_overall_summary.tsv",
            "Meaning": "Overall FDR-change statistics for the primary Individual-bin and Pairwise-bin analyses.",
            "Relative_Source_or_Derivation": "Reduced from <Condition>_Binning/<Group_Downstream|Pairwise_Downstream>/classic_vs_*_fdr_audit/classic_vs_binned_scan_peptide_fdr_wide.tsv when --include-fdr is used.",
        },
        {
            "File": "main_partition_fdr_interval_summary.tsv",
            "Meaning": "Median and distributional ΔFDR statistics across predefined full-database FDR intervals for the primary comparison.",
            "Relative_Source_or_Derivation": "Generated from the primary classic-vs-binned FDR audit tables using the interval boundaries defined in this collator.",
        },
        {
            "File": "main_partition_fdr_threshold_summary.tsv",
            "Meaning": "Pass/fail transition counts at FDR thresholds 0.01 and 0.05 for the primary comparison.",
            "Relative_Source_or_Derivation": "Generated from the primary classic-vs-binned FDR audit tables.",
        },
        {
            "File": "size_series_fdr_overall_summary.tsv",
            "Meaning": "Overall FDR-change statistics across the Individual-bin-size series.",
            "Relative_Source_or_Derivation": "Reduced from <Condition>_Groups_<LABEL>/Group_Downstream/classic_vs_group_fdr_audit/classic_vs_binned_scan_peptide_fdr_wide.tsv when --include-fdr is used.",
        },
        {
            "File": "size_series_fdr_interval_summary.tsv",
            "Meaning": "Median and distributional ΔFDR statistics across predefined full-database FDR intervals for each bin size.",
            "Relative_Source_or_Derivation": "Generated from the size-series classic-vs-binned FDR audit tables using the interval boundaries defined in this collator.",
        },
        {
            "File": "size_series_fdr_threshold_summary.tsv",
            "Meaning": "Pass/fail transition counts at FDR thresholds 0.01 and 0.05 for each bin size.",
            "Relative_Source_or_Derivation": "Generated from the size-series classic-vs-binned FDR audit tables.",
        },
        {
            "File": "*_fdr_manifest.tsv",
            "Meaning": "Provenance for the large FDR audit inputs reduced by the collator.",
            "Relative_Source_or_Derivation": "Generated when --include-fdr is used; records project-relative input paths, input/analysis row counts, statistic, and shared-row filtering.",
        },
        {
            "File": "results_guide.tsv",
            "Meaning": "This run-specific file guide.",
            "Relative_Source_or_Derivation": "Generated by collate_binning_results.py.",
        },
    ]
    write_tsv(
        out_dir / "results_guide.tsv",
        rows,
        ["File", "Meaning", "Relative_Source_or_Derivation"],
    )


def summarize_run(
    run: NumberedRun,
    results_root: Path,
    project_root: Path,
    strict: bool,
    include_fdr: bool = False,
) -> bool:
    print(f"\n[{run.name}] condition={run.condition}, replicate={run.replicate}")
    out_dir = results_root / run.name
    recovery_dir = out_dir / "recovery_inputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest: list[dict[str, object]] = []
    main_recovery: list[dict[str, object]] = []
    size_recovery: list[dict[str, object]] = []
    main_gains: list[dict[str, object]] = []
    size_gains: list[dict[str, object]] = []
    master: list[dict[str, object]] = []

    try:
        baseline, base_paths = baseline_counts(run, project_root)
    except (OSError, ValueError, FileNotFoundError) as exc:
        warn(f"{run.name}: baseline could not be read ({exc}); skipping run")
        return False

    manifest.extend([
        {
            "Run": run.name,
            "Category": "baseline",
            "Analysis": "full_database",
            "SourcePath": base_paths["BaselinePeptidePath"],
            "CopiedPath": "",
            "ValidationStatus": "PASS",
            "Notes": (f"source rows={baseline['BaselinePeptideSourceRows']}; " f"distinct Peptide={baseline['BaselinePeptides']}; " f"distinct PeptideFlanked={baseline['BaselinePeptideFlanked']}"),
        },
        {
            "Run": run.name,
            "Category": "baseline",
            "Analysis": "full_database",
            "SourcePath": base_paths["BaselineProteinPath"],
            "CopiedPath": "",
            "ValidationStatus": "PASS",
            "Notes": (f"source rows={baseline['BaselineProteinSourceRows']}; " f"distinct Protein={baseline['BaselineProteins']}"),
        },
    ])

    main_root = run.path / f"{run.condition}_Binning"
    try:
        main_group_count = count_partition_directories(
            main_root, PARTITION_KIND_INDIVIDUAL
        )
        main_pair_count = count_partition_directories(
            main_root, PARTITION_KIND_PAIRWISE
        )
    except (OSError, ValueError, FileNotFoundError) as exc:
        warn(f"{run.name}: primary partition structure could not be determined ({exc})")
        return False

    main_group_structure = partition_structure(
        main_group_count, PARTITION_KIND_INDIVIDUAL
    )
    main_pair_structure = partition_structure(
        main_pair_count, PARTITION_KIND_PAIRWISE
    )

    main_specs = [
        (
            PARTITION_KIND_INDIVIDUAL,
            main_group_structure,
            main_root / RECOVERY_GROUP_REL,
            "main_individual_recovery_summary.tsv",
        ),
        (
            PARTITION_KIND_PAIRWISE,
            main_pair_structure,
            main_root / RECOVERY_PAIR_REL,
            "main_pairwise_recovery_summary.tsv",
        ),
    ]

    for kind, structure, source, copy_name in main_specs:
        try:
            recovery, rec_manifest = load_recovery(source, kind)
            main_recovery.append({
                "Dataset": run.condition,
                "Run": run.name,
                "PartitionStructure": structure,
                **recovery,
            })
            copied = recovery_dir / copy_name
            safe_copy(source, copied)
            manifest.append({
                "Run": run.name,
                "Category": "recovery_input",
                "Analysis": structure.replace(" ", "_"),
                "SourcePath": project_relative(source, project_root),
                "CopiedPath": project_relative(copied, project_root),
                "ValidationStatus": rec_manifest["ValidationStatus"],
                "Notes": (
                    f"metrics={rec_manifest['MetricRows']}; "
                    f"{rec_manifest['RecoveredMetric']}/{rec_manifest['TotalMetric']}"
                ),
            })
            master.append({
                "Dataset": run.condition,
                "Run": run.name,
                "AnalysisFamily": "Pairwise",
                "PartitionStructure": structure,
                "BinSize": "",
                **recovery,
                **baseline,
            })
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(f"{run.name}: {structure} recovery unavailable ({exc})")
            if strict:
                return False

    size_specs: list[tuple[str, Path, str]] = []
    for size, condition_root in discover_size_series(run):
        try:
            group_count = count_partition_directories(
                condition_root, PARTITION_KIND_INDIVIDUAL
            )
            structure = partition_structure(
                group_count, PARTITION_KIND_INDIVIDUAL
            )
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(f"{run.name}: bin size {size} partition structure unavailable ({exc})")
            if strict:
                return False
            continue

        size_specs.append((size, condition_root, structure))
        source = condition_root / RECOVERY_GROUP_REL
        try:
            recovery, rec_manifest = load_recovery(
                source, PARTITION_KIND_INDIVIDUAL
            )
            size_recovery.append({
                "Dataset": run.condition,
                "Run": run.name,
                "BinSize": size,
                "PartitionStructure": structure,
                **recovery,
            })
            copied = recovery_dir / f"group_{size}_recovery_summary.tsv"
            safe_copy(source, copied)
            manifest.append({
                "Run": run.name,
                "Category": "recovery_input",
                "Analysis": f"size_{size}_{structure.replace(' ', '_')}",
                "SourcePath": project_relative(source, project_root),
                "CopiedPath": project_relative(copied, project_root),
                "ValidationStatus": rec_manifest["ValidationStatus"],
                "Notes": (
                    f"metrics={rec_manifest['MetricRows']}; "
                    f"{rec_manifest['RecoveredMetric']}/{rec_manifest['TotalMetric']}"
                ),
            })
            master.append({
                "Dataset": run.condition,
                "Run": run.name,
                "AnalysisFamily": "size_series",
                "PartitionStructure": structure,
                "BinSize": size,
                **recovery,
                **baseline,
            })
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(f"{run.name}: bin size {size} recovery unavailable ({exc})")
            if strict:
                return False

    for structure, downstream in [
        (main_group_structure, main_root / "Group_Downstream"),
        (main_pair_structure, main_root / "Pairwise_Downstream"),
    ]:
        try:
            row, man = gain_row(
                run, downstream, structure, baseline, base_paths, project_root
            )
            main_gains.append(row)
            manifest.extend(man)
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(f"{run.name}: {structure} peptide/protein gain unavailable ({exc})")
            if strict:
                return False

    for size, condition_root, structure in size_specs:
        try:
            row, man = gain_row(
                run,
                condition_root / "Group_Downstream",
                structure,
                baseline,
                base_paths,
                project_root,
                bin_size=size,
            )
            size_gains.append(row)
            manifest.extend(man)
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(
                f"{run.name}: bin size {size} peptide/protein gain unavailable ({exc})"
            )
            if strict:
                return False

    main_gain_lookup = {
        str(r["PartitionStructure"]): r for r in main_gains
    }
    size_gain_lookup = {
        str(r["BinSize"]): r for r in size_gains
    }
    for row in master:
        gain = None
        if row["AnalysisFamily"] == "Pairwise":
            gain = main_gain_lookup.get(str(row["PartitionStructure"]))
        elif row["AnalysisFamily"] == "size_series":
            gain = size_gain_lookup.get(str(row["BinSize"]))
        if gain:
            for key in [
                "CorrectedPeptides",
                "CorrectedPeptideFlanked",
                "CorrectedProteins",
                "PeptideGain",
                "PeptideFlankedGain",
                "ProteinGain",
            ]:
                row[key] = gain[key]

    main_recovery_fields = [
        "Dataset", "Run", "PartitionStructure", "TotalReferenceScans",
        "RecoveredScans", "UnrecoveredScans", "RecoveredFraction",
        "UnrecoveredFraction", "RecoveryPercent", "LossPercent",
    ]
    size_recovery_fields = [
        "Dataset", "Run", "BinSize", "PartitionStructure",
        "TotalReferenceScans", "RecoveredScans", "UnrecoveredScans",
        "RecoveredFraction", "UnrecoveredFraction",
        "RecoveryPercent", "LossPercent",
    ]
    main_gain_fields = [
        "Dataset", "Run", "PartitionStructure", "DownstreamRoot",
        "CorrectedPeptides", "CorrectedPeptideFlanked", "CorrectedProteins",
        "BaselinePeptides", "BaselinePeptideFlanked", "BaselineProteins",
        "PeptideGain", "PeptideFlankedGain", "ProteinGain",
        "PeptideSourceRows", "ProteinSourceRows", "PeptidePath", "ProteinPath",
        "BaselinePeptidePath", "BaselineProteinPath",
    ]
    size_gain_fields = [
        "Dataset", "Run", "BinSize", "PartitionStructure",
        "DownstreamRoot", "CorrectedPeptides", "CorrectedPeptideFlanked",
        "CorrectedProteins", "BaselinePeptides", "BaselinePeptideFlanked",
        "BaselineProteins", "PeptideGain", "PeptideFlankedGain", "ProteinGain",
        "PeptideSourceRows", "ProteinSourceRows", "PeptidePath", "ProteinPath",
        "BaselinePeptidePath", "BaselineProteinPath",
    ]
    master_fields = [
        "Dataset", "Run", "AnalysisFamily", "PartitionStructure", "BinSize",
        "TotalReferenceScans", "RecoveredScans", "UnrecoveredScans",
        "RecoveredFraction", "UnrecoveredFraction", "RecoveryPercent",
        "LossPercent", "BaselinePeptides", "BaselinePeptideFlanked", "BaselineProteins",
        "CorrectedPeptides", "CorrectedPeptideFlanked", "CorrectedProteins",
        "PeptideGain", "PeptideFlankedGain", "ProteinGain",
    ]
    manifest_fields = [
        "Run", "Category", "Analysis", "SourcePath", "CopiedPath",
        "ValidationStatus", "Notes",
    ]

    write_tsv(
        out_dir / "main_partition_scan_recovery_summary.tsv",
        main_recovery, main_recovery_fields
    )
    write_tsv(
        out_dir / "group_size_scan_recovery_summary.tsv",
        size_recovery, size_recovery_fields
    )
    write_tsv(
        out_dir / "main_partition_peptide_protein_gains.tsv",
        main_gains, main_gain_fields
    )
    write_tsv(
        out_dir / "group_size_peptide_protein_gains.tsv",
        size_gains, size_gain_fields
    )
    write_tsv(
        out_dir / "binning_results_summary.tsv",
        master, master_fields
    )
    write_tsv(
        out_dir / "input_manifest.tsv",
        manifest, manifest_fields
    )
    write_results_guide(run, out_dir)

    if include_fdr:
        fdr_ok = collate_fdr_for_run(run, out_dir, project_root, strict)
        if not fdr_ok and strict:
            return False

    print(
        f"[OK] {run.name}: recovery={len(main_recovery) + len(size_recovery)} rows; "
        f"gains={len(main_gains) + len(size_gains)} rows -> {out_dir}"
    )
    return True



def _matched_analysis_order(name: str) -> tuple[int, object]:
    if name == "Individual":
        return (0, 0)
    if name == "Pairwise":
        return (1, 0)
    size_key = size_label_sort_key(name)
    if size_key[0] == 0:
        return (2, size_key[1])
    return (3, size_key[1])


def collate_matched_formula_outputs(
    project_root: Path,
    results_root: Path,
    runs: list[NumberedRun],
    matched_root: Path,
    strict: bool,
) -> bool:
    """Collate matched-formula downstream crosstabs into one compact table."""
    conditions = sorted({r.condition for r in runs}, key=natural_key)
    rows: list[dict[str, object]] = []
    manifest: list[dict[str, object]] = []

    for condition in conditions:
        condition_runs = sorted([r for r in runs if r.condition == condition], key=lambda r: r.replicate)
        if not condition_runs:
            continue
        condition_out = matched_root / condition
        if not condition_out.is_dir():
            warn(f"{condition}: matched output directory not found: {condition_out}")
            if strict:
                return False
            continue

        full_path = project_root / f"{condition}_Full"
        try:
            baseline, base_paths = baseline_counts_from_full(full_path, project_root)
        except (OSError, ValueError, FileNotFoundError) as exc:
            warn(f"{condition}: matched baseline could not be read ({exc})")
            if strict:
                return False
            continue

        analysis_names = sorted(
            [p.name for p in condition_out.iterdir() if p.is_dir() and not p.name.startswith("_")],
            key=_matched_analysis_order,
        )
        for analysis_name in analysis_names:
            if analysis_name in {"Individual", "Pairwise"}:
                analysis_family = "Pairwise"
                formula_design = f"{analysis_name} Bins"
                bin_size = ""
                partition_root = condition_runs[0].path / f"{condition}_Binning"
                kind = PARTITION_KIND_INDIVIDUAL if analysis_name == "Individual" else PARTITION_KIND_PAIRWISE
            else:
                analysis_family = "Bin Size Series"
                formula_design = f"{analysis_name} Individual Bins"
                bin_size = analysis_name
                partition_root = condition_runs[0].path / f"{condition}_Groups_{analysis_name}"
                kind = PARTITION_KIND_INDIVIDUAL

            try:
                structure = partition_structure(count_partition_directories(partition_root, kind), kind)
            except (OSError, ValueError, FileNotFoundError) as exc:
                warn(f"{condition} matched {analysis_name}: partition structure unavailable ({exc})")
                if strict:
                    return False
                structure = ""

            analysis_root = condition_out / analysis_name
            peptide_path = analysis_root / "protein_rollup" / PEPTIDE_CROSSTAB_NAME
            protein_path = analysis_root / "protein_rollup" / PROTEIN_CROSSTAB_NAME
            try:
                pep = count_peptide_entities(peptide_path)
                pro = count_protein_entities(protein_path)
            except (OSError, ValueError, FileNotFoundError) as exc:
                warn(f"{condition} matched {analysis_name}: gain outputs unavailable ({exc})")
                if strict:
                    return False
                continue

            row = {
                "Dataset": condition,
                "ApplicationMode": "Matched formulas applied to source samples",
                "AnalysisFamily": analysis_family,
                "FormulaDesign": formula_design,
                "PartitionStructure": structure,
                "BinSize": bin_size,
                "MatchedOutputRoot": project_relative(analysis_root, project_root),
                "BaselinePeptides": baseline["BaselinePeptides"],
                "CorrectedPeptides": pep["Peptides"],
                "PeptideGain": pep["Peptides"] - baseline["BaselinePeptides"],
                "BaselinePeptideFlanked": baseline["BaselinePeptideFlanked"],
                "CorrectedPeptideFlanked": pep["PeptideFlanked"],
                "PeptideFlankedGain": pep["PeptideFlanked"] - baseline["BaselinePeptideFlanked"],
                "BaselineProteins": baseline["BaselineProteins"],
                "CorrectedProteins": pro["Proteins"],
                "ProteinGain": pro["Proteins"] - baseline["BaselineProteins"],
                "PeptideSourceRows": pep["SourceRows"],
                "ProteinSourceRows": pro["SourceRows"],
                "PeptidePath": project_relative(peptide_path, project_root),
                "ProteinPath": project_relative(protein_path, project_root),
                "BaselinePeptidePath": base_paths["BaselinePeptidePath"],
                "BaselineProteinPath": base_paths["BaselineProteinPath"],
            }
            rows.append(row)
            manifest.extend([
                {
                    "Dataset": condition, "Analysis": analysis_name, "FileType": "peptide",
                    "SourcePath": row["PeptidePath"], "ValidationStatus": "PASS",
                    "CountingRule": "distinct Peptide; distinct PeptideFlanked retained as QC",
                    "Notes": f"source rows={pep['SourceRows']}; distinct Peptide={pep['Peptides']}; distinct PeptideFlanked={pep['PeptideFlanked']}",
                },
                {
                    "Dataset": condition, "Analysis": analysis_name, "FileType": "protein",
                    "SourcePath": row["ProteinPath"], "ValidationStatus": "PASS",
                    "CountingRule": "distinct Protein",
                    "Notes": f"source rows={pro['SourceRows']}; distinct Protein={pro['Proteins']}",
                },
            ])

    if not rows:
        warn(f"No matched-formula peptide/protein outputs were collated from {matched_root}")
        return not strict

    fields = [
        "Dataset", "ApplicationMode", "AnalysisFamily", "FormulaDesign", "PartitionStructure",
        "BinSize", "MatchedOutputRoot", "BaselinePeptides", "CorrectedPeptides", "PeptideGain",
        "BaselinePeptideFlanked", "CorrectedPeptideFlanked", "PeptideFlankedGain",
        "BaselineProteins", "CorrectedProteins", "ProteinGain", "PeptideSourceRows",
        "ProteinSourceRows", "PeptidePath", "ProteinPath", "BaselinePeptidePath", "BaselineProteinPath",
    ]
    manifest_fields = ["Dataset", "Analysis", "FileType", "SourcePath", "ValidationStatus", "CountingRule", "Notes"]
    write_tsv(results_root / "matched_formula_peptide_protein_gains.tsv", rows, fields)
    write_tsv(results_root / "matched_formula_input_manifest.tsv", manifest, manifest_fields)
    write_tsv(
        results_root / "matched_formula_results_guide.tsv",
        [
            {
                "File": "matched_formula_peptide_protein_gains.tsv",
                "Meaning": "Dataset-level peptide/protein changes after each numbered-series formula is applied only to its matched source sample and one combined downstream crosstab is generated.",
                "CountingRule": "Peptides = distinct Peptide; PeptideFlanked = distinct PeptideFlanked QC; proteins = distinct Protein.",
            },
            {
                "File": "matched_formula_input_manifest.tsv",
                "Meaning": "Provenance and counting-rule audit for matched-formula crosstab inputs.",
                "CountingRule": "Recorded per file.",
            },
        ],
        ["File", "Meaning", "CountingRule"],
    )
    print(f"[OK] matched formula outputs: {len(rows)} peptide/protein gain rows -> {results_root / 'matched_formula_peptide_protein_gains.tsv'}")
    return True

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collate scan-recovery, peptide/protein gain, and optional FDR "
            "summaries from numbered Binning runs into results/<run>/ directories."
        )
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=None,
        help=(
            "Binning data root containing <condition>_<integer> runs and matching "
            "<condition>_Full baselines. Default: <Binning>/data inferred from "
            "this script's installed location."
        ),
    )
    parser.add_argument(
        "-o", "--output-dir", default=None,
        help="Results root (default: <root>/results)",
    )
    parser.add_argument(
        "--run", action="append", default=[],
        help="Only collate this numbered run (repeatable; e.g. --run Fecal_1)",
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="Treat any missing per-analysis input as failure for that numbered run.",
    )
    parser.add_argument(
        "--include-fdr",
        action="store_true",
        help=(
            "Also reduce large FDR audit tables into compact overall, interval, "
            "and threshold summaries. This is slower than recovery/gain collation."
        ),
    )
    parser.add_argument(
        "--matched-output-dir",
        default=None,
        help=(
            "Matched formula output root. Default: <root>/Matched_Formula_Output. "
            "When present, matched peptide/protein gains are collated automatically."
        ),
    )
    parser.add_argument(
        "--skip-matched",
        action="store_true",
        help="Do not collate Matched_Formula_Output even when it exists.",
    )
    parser.add_argument(
        "--require-matched",
        action="store_true",
        help="Fail if matched-formula outputs are absent or incomplete.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {SCRIPT_VERSION}")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    default_root = Path(__file__).resolve().parents[2] / "data"
    root = (
        Path(args.root).expanduser().resolve()
        if args.root
        else default_root.resolve()
    )
    results_root = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else root / "results"
    )
    matched_root = (
        Path(args.matched_output_dir).expanduser().resolve()
        if args.matched_output_dir
        else root / MATCHED_OUTPUT_DIRNAME
    )

    if not root.is_dir():
        print(f"ERROR: root directory does not exist: {root}", file=sys.stderr)
        return 2

    runs = discover_numbered_runs(root)
    if args.run:
        requested = set(args.run)
        known = {r.name for r in runs}
        missing = sorted(requested - known, key=natural_key)
        if missing:
            print(
                "ERROR: requested run(s) were not discovered with matching *_Full baselines: "
                + ", ".join(missing),
                file=sys.stderr,
            )
            return 2
        runs = [r for r in runs if r.name in requested]

    if not runs:
        print(
            "ERROR: no numbered <condition>_<integer> directories with matching "
            "<condition>_Full baselines were discovered.",
            file=sys.stderr,
        )
        return 1

    print("============================================================")
    print("BINNING RESULT COLLATION")
    print("============================================================")
    print(f"Script version : {SCRIPT_VERSION}")
    print(f"Root           : {root}")
    print(f"Results root   : {results_root}")
    print(f"Matched root   : {matched_root}")
    print("Discovered runs: " + ", ".join(r.name for r in runs))

    successes = 0
    for run in runs:
        if summarize_run(run, results_root, root, strict=args.strict, include_fdr=args.include_fdr):
            successes += 1

    matched_ok = True
    if not args.skip_matched:
        if matched_root.is_dir():
            matched_ok = collate_matched_formula_outputs(
                root, results_root, runs, matched_root, strict=(args.strict or args.require_matched)
            )
        elif args.require_matched:
            warn(f"Required matched output directory not found: {matched_root}")
            matched_ok = False
        else:
            print(f"Matched formula collation skipped: directory not present: {matched_root}")

    print(f"\nDone: {successes}/{len(runs)} numbered run(s) collated; matched_ok={matched_ok}.")
    return 0 if successes == len(runs) and matched_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
