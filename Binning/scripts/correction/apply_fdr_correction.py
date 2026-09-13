#!/usr/bin/env python3
"""
apply_fdr_correction.py

Apply a previously derived monotonic FDR correction formula to every
``*_withsyn.tsv`` file in the ``fdr_esti`` directory of a classical/full run.

The original files are never modified. Corrected copies are written to a
separate output directory. In each corrected copy:

- ``OriginalFDR`` preserves the source FDR value.
- ``FDR`` is replaced with the corrected value for downstream pipeline use.
- ``CorrectedFDR`` duplicates the corrected value with an explicit name.
- ``FDRCorrectionDelta`` is CorrectedFDR - OriginalFDR.
- Formula provenance columns identify the applied condition and formula.

Because a withsyn table may contain multiple protein rows for one ScanNum,
all rows are corrected, but threshold gains and losses are counted once per
ScanNum. Scan-level pass status is based on the minimum valid FDR represented
for that scan; this is equivalent to the usual case where protein-expanded
rows carry identical FDR values.

Default input columns
---------------------
ScanNum
Peptide
Protein
FDR

Primary outputs
---------------
<outdir>/*_withsyn.tsv
    Corrected copies retaining the original filenames.

<outdir>/correction_application_summary.tsv
    Per-file and combined threshold-crossing statistics, counted once per scan.

<outdir>/correction_application_metadata.tsv
    Formula, path, row-count, clipping, and consistency diagnostics.

<outdir>/gained_scans.tsv
    One row per newly passing ScanNum and threshold, with peptide annotations.

The output directory defaults to ``corrected_fdr_esti`` beside the
``correctional_formula`` directory containing the formula. Use ``--outdir`` to
override this placement.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

SCRIPT_VERSION = "1.0.0"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Apply a monotonic correction formula to *_withsyn.tsv files from "
            "a classical/full run while preserving the originals."
        )
    )
    p.add_argument(
        "--formula",
        required=True,
        type=Path,
        help="correction_formula.json produced by derive_fdr_correction.py.",
    )
    p.add_argument(
        "--full_dir",
        required=True,
        type=Path,
        help="Classical/full-run directory containing fdr_esti.",
    )
    p.add_argument(
        "--input_dir",
        type=Path,
        default=None,
        help="Optional direct input directory. Default: <full_dir>/fdr_esti.",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=(
            "Corrected output directory. Default: corrected_fdr_esti beside the "
            "correctional_formula directory containing --formula."
        ),
    )
    p.add_argument(
        "--pattern",
        default="*_withsyn.tsv",
        help="Input filename glob. Default: *_withsyn.tsv",
    )
    p.add_argument(
        "--thresholds",
        default="0.01,0.05",
        help="Comma-separated FDR thresholds for gain/loss summaries. Default: 0.01,0.05",
    )
    p.add_argument("--scan_col", default="ScanNum")
    p.add_argument("--peptide_col", default="Peptide")
    p.add_argument("--protein_col", default="Protein")
    p.add_argument("--fdr_col", default="FDR")
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Permit replacement of files already present in the output directory.",
    )
    p.add_argument(
        "--no_gained_scan_table",
        action="store_true",
        help="Do not write gained_scans.tsv.",
    )
    return p.parse_args()


def parse_thresholds(spec: str) -> list[float]:
    values = sorted({float(x.strip()) for x in spec.split(",") if x.strip()})
    if not values:
        raise ValueError("At least one threshold is required.")
    if any((not np.isfinite(x)) or x < 0.0 or x > 1.0 for x in values):
        raise ValueError("Thresholds must be finite and within [0, 1].")
    return values


def infer_default_outdir(formula_path: Path) -> Path:
    formula_path = formula_path.resolve()
    for parent in formula_path.parents:
        if parent.name == "correctional_formula":
            return parent.parent / "corrected_fdr_esti"
    return formula_path.parent / "corrected_fdr_esti"


def load_formula(path: Path) -> tuple[dict, np.ndarray, np.ndarray]:
    with path.open("r", encoding="utf-8") as handle:
        formula = json.load(handle)

    required = {"fdr_anchor", "corrected_fdr_anchor"}
    missing = sorted(required - set(formula))
    if missing:
        raise ValueError(f"Formula JSON is missing required fields: {missing}")

    x = np.asarray(formula["fdr_anchor"], dtype=float)
    y = np.asarray(formula["corrected_fdr_anchor"], dtype=float)

    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError("fdr_anchor and corrected_fdr_anchor must be equal-length 1D arrays.")
    if len(x) < 2:
        raise ValueError("At least two formula anchors are required.")
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(y)):
        raise ValueError("Formula anchors must be finite.")
    if np.any(np.diff(x) <= 0):
        raise ValueError("fdr_anchor must be strictly increasing.")
    if np.any(np.diff(y) < -1e-12):
        raise ValueError("corrected_fdr_anchor must be nondecreasing.")
    if not (np.isclose(x[0], 0.0) and np.isclose(x[-1], 1.0)):
        raise ValueError("Formula anchors must span the complete [0, 1] FDR domain.")
    if np.any((y < -1e-12) | (y > 1.0 + 1e-12)):
        raise ValueError("corrected_fdr_anchor values must lie within [0, 1].")

    return formula, x, np.clip(y, 0.0, 1.0)


def atomic_to_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        df.to_csv(temp_path, sep="\t", index=False)
        os.replace(temp_path, path)
    finally:
        temp_path.unlink(missing_ok=True)


def unique_join(values: pd.Series, max_items: int = 100) -> str:
    unique = sorted(
        {
            str(v)
            for v in values
            if pd.notna(v) and str(v).strip() and str(v).strip().lower() != "nan"
        }
    )
    if len(unique) > max_items:
        return ";".join(unique[:max_items]) + f";...(+{len(unique) - max_items} more)"
    return ";".join(unique)


def prepare_scan_table(
    df: pd.DataFrame,
    scan_col: str,
    peptide_col: str,
) -> pd.DataFrame:
    work = df.loc[df[scan_col].notna() & df["OriginalFDR"].notna()].copy()
    if work.empty:
        return pd.DataFrame(
            columns=[
                "ScanNum",
                "OriginalFDR",
                "CorrectedFDR",
                "Peptide",
                "SourceRows",
                "DistinctOriginalFDRValues",
            ]
        )

    aggregations: dict[str, tuple[str, object]] = {
        "OriginalFDR": ("OriginalFDR", "min"),
        "CorrectedFDR": ("CorrectedFDR", "min"),
        "SourceRows": (scan_col, "size"),
        "DistinctOriginalFDRValues": ("OriginalFDR", "nunique"),
    }
    if peptide_col in work.columns:
        aggregations["Peptide"] = (peptide_col, unique_join)

    scans = work.groupby(scan_col, as_index=False, dropna=False).agg(**aggregations)
    scans = scans.rename(columns={scan_col: "ScanNum"})
    if "Peptide" not in scans.columns:
        scans["Peptide"] = ""
    return scans


def threshold_rows(
    scans: pd.DataFrame,
    source_file: str,
    thresholds: Sequence[float],
) -> tuple[list[dict[str, object]], list[pd.DataFrame]]:
    rows: list[dict[str, object]] = []
    gained_tables: list[pd.DataFrame] = []

    for threshold in thresholds:
        original_pass = scans["OriginalFDR"] <= threshold
        corrected_pass = scans["CorrectedFDR"] <= threshold
        gained = (~original_pass) & corrected_pass
        lost = original_pass & (~corrected_pass)

        gained_sub = scans.loc[gained].copy()
        if not gained_sub.empty:
            gained_sub.insert(0, "SourceFile", source_file)
            gained_sub.insert(1, "Threshold", threshold)
            gained_sub["OriginalPass"] = False
            gained_sub["CorrectedPass"] = True
            gained_tables.append(gained_sub)

        rows.append(
            {
                "SourceFile": source_file,
                "Threshold": threshold,
                "UniqueScans": int(len(scans)),
                "OriginalPassScans": int(original_pass.sum()),
                "CorrectedPassScans": int(corrected_pass.sum()),
                "GainedScans": int(gained.sum()),
                "LostScans": int(lost.sum()),
                "NetPassScanChange": int(corrected_pass.sum() - original_pass.sum()),
                "PassBothScans": int((original_pass & corrected_pass).sum()),
                "FailBothScans": int((~original_pass & ~corrected_pass).sum()),
                "GainedUniquePeptideStrings": int(
                    gained_sub["Peptide"].replace("", pd.NA).nunique(dropna=True)
                ),
            }
        )
    return rows, gained_tables


def process_file(
    input_path: Path,
    output_path: Path,
    formula: dict,
    x: np.ndarray,
    y: np.ndarray,
    args: argparse.Namespace,
    thresholds: Sequence[float],
) -> tuple[list[dict[str, object]], dict[str, object], list[pd.DataFrame], pd.DataFrame]:
    df = pd.read_csv(input_path, sep="\t", low_memory=False)

    required = {args.scan_col, args.fdr_col}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{input_path} is missing required columns: {missing}")

    original_numeric = pd.to_numeric(df[args.fdr_col], errors="coerce")
    invalid_count = int(original_numeric.isna().sum())
    out_of_range = original_numeric.notna() & ~original_numeric.between(0.0, 1.0)
    out_of_range_count = int(out_of_range.sum())
    if out_of_range_count:
        examples = original_numeric.loc[out_of_range].head(5).tolist()
        raise ValueError(
            f"{input_path} contains {out_of_range_count} FDR values outside [0,1]; "
            f"examples: {examples}"
        )

    valid = original_numeric.notna()
    corrected = pd.Series(np.nan, index=df.index, dtype=float)
    corrected.loc[valid] = np.clip(
        np.interp(original_numeric.loc[valid].to_numpy(dtype=float), x, y),
        0.0,
        1.0,
    )

    # Preserve source FDR and make the corrected copy immediately usable by
    # downstream code that expects the canonical FDR column.
    if "OriginalFDR" in df.columns or "CorrectedFDR" in df.columns:
        raise ValueError(
            f"{input_path} already contains OriginalFDR or CorrectedFDR; "
            "refusing to ambiguously re-correct an existing output."
        )

    insert_at = list(df.columns).index(args.fdr_col)
    df.insert(insert_at, "OriginalFDR", original_numeric)
    df[args.fdr_col] = corrected
    df.insert(insert_at + 2, "CorrectedFDR", corrected)
    df.insert(insert_at + 3, "FDRCorrectionDelta", corrected - original_numeric)
    df.insert(
        insert_at + 4,
        "FDRCorrectionCondition",
        str(formula.get("condition_name", "")),
    )
    df.insert(
        insert_at + 5,
        "FDRCorrectionFormulaVersion",
        str(formula.get("script_version", "")),
    )

    atomic_to_csv(df, output_path)

    scans = prepare_scan_table(df, args.scan_col, args.peptide_col)
    summary_rows, gained_tables = threshold_rows(scans, input_path.name, thresholds)

    metadata = {
        "SourceFile": input_path.name,
        "InputPath": str(input_path.resolve()),
        "OutputPath": str(output_path.resolve()),
        "InputRows": int(len(df)),
        "RowsWithValidFDR": int(valid.sum()),
        "RowsWithMissingOrNonNumericFDR": invalid_count,
        "UniqueScansWithValidFDR": int(len(scans)),
        "ScansWithMultipleOriginalFDRValues": int(
            (scans["DistinctOriginalFDRValues"] > 1).sum()
        ),
        "MinimumOriginalFDR": float(original_numeric.min()) if valid.any() else np.nan,
        "MaximumOriginalFDR": float(original_numeric.max()) if valid.any() else np.nan,
        "MinimumCorrectedFDR": float(corrected.min()) if valid.any() else np.nan,
        "MaximumCorrectedFDR": float(corrected.max()) if valid.any() else np.nan,
        "RowsCorrectedToZero": int((corrected == 0.0).sum()),
        "RowsCorrectedToOne": int((corrected == 1.0).sum()),
        "RowsWithLowerCorrectedFDR": int((corrected < original_numeric).sum()),
        "RowsWithEqualCorrectedFDR": int(
            np.isclose(
                corrected.to_numpy(dtype=float),
                original_numeric.to_numpy(dtype=float),
                equal_nan=False,
            ).sum()
        ),
        "RowsWithHigherCorrectedFDR": int((corrected > original_numeric).sum()),
    }
    return summary_rows, metadata, gained_tables, scans


def add_combined_summary(
    scan_tables: list[pd.DataFrame],
    thresholds: Sequence[float],
) -> list[dict[str, object]]:
    if not scan_tables:
        return []

    combined = pd.concat(scan_tables, ignore_index=True)
    # Scan numbers can repeat across separate files. Treat SourceFile + ScanNum
    # as the combined identity for the combined summary.
    combined["ScanNum"] = (
        combined["SourceFile"].astype(str) + "|" + combined["ScanNum"].astype(str)
    )
    combined = combined.drop(columns=["SourceFile"], errors="ignore")
    rows, _ = threshold_rows(combined, "ALL_FILES", thresholds)
    return rows


def main() -> int:
    args = parse_args()
    formula_path = args.formula.expanduser().resolve()
    full_dir = args.full_dir.expanduser().resolve()
    input_dir = (
        args.input_dir.expanduser().resolve()
        if args.input_dir is not None
        else full_dir / "fdr_esti"
    )
    outdir = (
        args.outdir.expanduser().resolve()
        if args.outdir is not None
        else infer_default_outdir(formula_path)
    )

    if not formula_path.is_file():
        raise SystemExit(f"ERROR: formula not found: {formula_path}")
    if not input_dir.is_dir():
        raise SystemExit(f"ERROR: input directory not found: {input_dir}")
    if input_dir == outdir:
        raise SystemExit("ERROR: input and output directories must differ.")

    try:
        thresholds = parse_thresholds(args.thresholds)
        formula, x, y = load_formula(formula_path)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    input_files = sorted(input_dir.glob(args.pattern))
    if not input_files:
        raise SystemExit(
            f"ERROR: no files matching {args.pattern!r} were found in {input_dir}"
        )

    outdir.mkdir(parents=True, exist_ok=True)
    output_paths = [outdir / path.name for path in input_files]
    existing = [path for path in output_paths if path.exists()]
    if existing and not args.overwrite:
        preview = "\n".join(f"  {p}" for p in existing[:10])
        raise SystemExit(
            "ERROR: corrected output files already exist. Use --overwrite to replace them:\n"
            + preview
        )

    summary_rows: list[dict[str, object]] = []
    metadata_rows: list[dict[str, object]] = []
    gained_tables: list[pd.DataFrame] = []
    combined_scan_tables: list[pd.DataFrame] = []

    for input_path, output_path in zip(input_files, output_paths):
        try:
            file_summary, file_metadata, file_gained, scan_table = process_file(
                input_path,
                output_path,
                formula,
                x,
                y,
                args,
                thresholds,
            )
        except (ValueError, OSError, pd.errors.ParserError) as exc:
            raise SystemExit(f"ERROR: {exc}") from exc

        summary_rows.extend(file_summary)
        metadata_rows.append(file_metadata)
        gained_tables.extend(file_gained)
        scan_table = scan_table.copy()
        scan_table.insert(0, "SourceFile", input_path.name)
        combined_scan_tables.append(scan_table)
        print(
            f"[OK] {input_path.name}: {len(scan_table):,} unique scans -> {output_path}"
        )

    combined_summary = add_combined_summary(combined_scan_tables, thresholds)
    summary_rows.extend(combined_summary)

    summary = pd.DataFrame(summary_rows)
    atomic_to_csv(summary, outdir / "correction_application_summary.tsv")

    run_metadata = pd.DataFrame(
        [
            ("script", Path(__file__).name),
            ("script_version", SCRIPT_VERSION),
            ("generated_utc", datetime.now(timezone.utc).isoformat()),
            ("python_version", platform.python_version()),
            ("pandas_version", pd.__version__),
            ("numpy_version", np.__version__),
            ("formula_path", str(formula_path)),
            ("formula_condition_name", formula.get("condition_name", "")),
            ("formula_derivation_script_version", formula.get("script_version", "")),
            ("formula_model", formula.get("model", "")),
            ("formula_monotonic_method", formula.get("monotonic_method", "")),
            ("formula_anchor_count", len(x)),
            ("input_directory", str(input_dir)),
            ("input_pattern", args.pattern),
            ("output_directory", str(outdir)),
            ("input_file_count", len(input_files)),
            ("thresholds", ",".join(f"{x:g}" for x in thresholds)),
            ("scan_column", args.scan_col),
            ("peptide_column", args.peptide_col),
            ("protein_column", args.protein_col),
            ("source_fdr_column", args.fdr_col),
            ("source_files_preserved", 1),
            ("output_fdr_column_replaced_with_corrected_values", 1),
            ("original_fdr_preserved_as", "OriginalFDR"),
        ],
        columns=["Metric", "Value"],
    )
    atomic_to_csv(run_metadata, outdir / "correction_application_metadata.tsv")

    file_metadata = pd.DataFrame(metadata_rows)
    atomic_to_csv(file_metadata, outdir / "correction_file_diagnostics.tsv")

    if not args.no_gained_scan_table:
        if gained_tables:
            gained = pd.concat(gained_tables, ignore_index=True)
            preferred = [
                "SourceFile",
                "Threshold",
                "ScanNum",
                "Peptide",
                "OriginalFDR",
                "CorrectedFDR",
                "SourceRows",
                "DistinctOriginalFDRValues",
                "OriginalPass",
                "CorrectedPass",
            ]
            columns = [c for c in preferred if c in gained.columns] + [
                c for c in gained.columns if c not in preferred
            ]
            gained = gained[columns].sort_values(
                ["Threshold", "SourceFile", "OriginalFDR", "ScanNum"]
            )
        else:
            gained = pd.DataFrame(
                columns=[
                    "SourceFile",
                    "Threshold",
                    "ScanNum",
                    "Peptide",
                    "OriginalFDR",
                    "CorrectedFDR",
                    "SourceRows",
                    "DistinctOriginalFDRValues",
                    "OriginalPass",
                    "CorrectedPass",
                ]
            )
        atomic_to_csv(gained, outdir / "gained_scans.tsv")

    print("\nFDR correction application completed successfully.")
    print(f"Formula condition:  {formula.get('condition_name', '')}")
    print(f"Input files:        {len(input_files):,}")
    print(f"Output directory:   {outdir}")
    for threshold in thresholds:
        combined_row = summary.loc[
            (summary["SourceFile"] == "ALL_FILES")
            & np.isclose(summary["Threshold"], threshold)
        ].iloc[0]
        print(
            f"FDR <= {threshold:g}:       "
            f"{int(combined_row['OriginalPassScans']):,} original -> "
            f"{int(combined_row['CorrectedPassScans']):,} corrected "
            f"(+{int(combined_row['GainedScans']):,} gained; "
            f"-{int(combined_row['LostScans']):,} lost)"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        raise SystemExit(130)
