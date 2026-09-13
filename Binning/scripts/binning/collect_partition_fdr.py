#!/usr/bin/env python3
"""
collect_partition_bin_fdr_values.py

Collect and summarize partition FDR values after partition-level post-processing.

Typical upstream workflow:
    partition SIC preparation
    -> optional/classical-guided partition culling
    -> partition-level FDR recalculation
    -> optional SYN expansion limited to surviving scan/peptide assignments

Reads:
    <dataset_dir>/*_partition_i_j/<fdr_dir_name>/<pattern>

Default read target:
    <dataset_dir>/*_partition_*_*/fdr_esti/*_withsyn.tsv

Writes, by default into <dataset_dir>/Binning_Downstream/:
    partition_bin_fdr_long.tsv
        One row per partition input row, with RunFolder/PartitionGroup metadata.
        The input FDR column is renamed to PartitionFDR to make clear that it is
        the FDR estimated within a single partition context.

    partition_bin_fdr_annotated.tsv
        Partition-level rows with binned FDR summary statistics merged onto each row.
        No generic column named FDR is written, because a scan/peptide/protein
        can have multiple partition FDR observations.

    scan_core_peptide_fdr_summary.tsv
        Primary comparison input. One row per ScanNum + CorePeptide, calculated
        from exactly one normalized observation per represented partition bin.

    scan_peptide_fdr_summary.tsv
        Diagnostic exact-flanked summary: one row per ScanNum + Peptide.

    scan_peptide_protein_fdr_summary.tsv
        One row per ScanNum + Peptide + Protein, summarizing FDR across
        represented bins.

Notes:
    - Peptide strings are preserved exactly except for whitespace stripping.
    - Modification notation is not removed.
    - Binned statistics are explicit: BinnedFDRMean, BinnedFDRMedian,
      BinnedFDRMin, BinnedFDRMax, BinnedFDRStd, BinnedFDRRange, and
      BinnedFDRCount.
    - The script intentionally does not choose a primary/representative FDR.
      Downstream plotting, formula fitting, or filtering should explicitly
      choose which binned statistic to use.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, List

import pandas as pd


PAIR_RE = re.compile(r".*_pair_(\d+)_(\d+)$")
GROUP_RE = re.compile(r".*_group_(\d+)$")


def parse_partition(name: str, partition_type: str):
    if partition_type == "pair":
        m = PAIR_RE.match(name)
        if not m:
            return None
        return {
            "PartitionType": "pair",
            "PartitionID": f"{m.group(1)}_{m.group(2)}",
            "PairGroupA": int(m.group(1)),
            "PairGroupB": int(m.group(2)),
        }
    if partition_type == "group":
        m = GROUP_RE.match(name)
        if not m:
            return None
        return {
            "PartitionType": "group",
            "PartitionID": m.group(1),
            "Group": int(m.group(1)),
        }
    raise ValueError(f"Unsupported partition_type: {partition_type}")


def normalize_core_peptide(value: object, strip_mods: bool = False) -> str:
    """Strip PHRP/MS-GF+ flanking notation; optionally remove PTM symbols."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.count(".") >= 2:
        first = text.find(".")
        last = text.rfind(".")
        if first < last:
            text = text[first + 1:last]
    text = text.strip()
    if strip_mods:
        text = re.sub(r"[^A-Z]", "", text)
    return text

def normalize_bool_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def require_columns(df: pd.DataFrame, required: set[str], path: Path) -> None:
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: {path} missing required columns: {sorted(missing)}")


def ordered_existing_columns(df: pd.DataFrame, wanted: Iterable[str]) -> List[str]:
    return [c for c in wanted if c in df.columns]


def choose_representative_rows(long_df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """
    Choose one representative row per group for a compact downstream table.

    This option is provided only for convenience when users need a compact table.
    It does not choose a primary FDR. Binned FDR statistics are still reported in
    explicit columns, and the representative row only supplies metadata fields.
    Lower SpecEValue is preferred when available, then lower PartitionFDR.
    """
    work = long_df.copy()
    sort_cols = list(group_cols)
    ascending = [True] * len(sort_cols)

    if "SpecEValue" in work.columns:
        work["SpecEValue_num_for_sort"] = pd.to_numeric(work["SpecEValue"], errors="coerce")
        sort_cols.append("SpecEValue_num_for_sort")
        ascending.append(True)

    sort_cols.append("PartitionFDR_num")
    ascending.append(True)

    sort_cols.append("RunFolder")
    ascending.append(True)

    work = work.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    rep = work.groupby(group_cols, dropna=False, as_index=False).head(1).copy()
    rep = rep.drop(columns=["SpecEValue_num_for_sort"], errors="ignore")
    return rep


def summarize_fdr(long_df: pd.DataFrame, group_cols: list[str], prefix: str = "Binned") -> pd.DataFrame:
    summary = (
        long_df
        .groupby(group_cols, dropna=False)
        .agg(
            RowsAcrossBins=("RunFolder", "size"),
            PartitionBinsRepresented=("RunFolder", "nunique"),
            **{
                f"{prefix}FDRMean": ("PartitionFDR_num", "mean"),
                f"{prefix}FDRMedian": ("PartitionFDR_num", "median"),
                f"{prefix}FDRMin": ("PartitionFDR_num", "min"),
                f"{prefix}FDRMax": ("PartitionFDR_num", "max"),
                f"{prefix}FDRStd": ("PartitionFDR_num", "std"),
                f"{prefix}FDRCount": ("PartitionFDR_num", "count"),
                "BinnedRunFolders": ("RunFolder", lambda x: ";".join(sorted(x.unique()))),
            }
        )
        .reset_index()
    )

    summary[f"{prefix}FDRRange"] = summary[f"{prefix}FDRMax"] - summary[f"{prefix}FDRMin"]

    if "MSGFScore" in long_df.columns:
        tmp = long_df.copy()
        tmp["MSGFScore_num"] = pd.to_numeric(tmp["MSGFScore"], errors="coerce")
        msgf = (
            tmp.groupby(group_cols, dropna=False)
            .agg(
                MSGFScoreMin=("MSGFScore_num", "min"),
                MSGFScoreMax=("MSGFScore_num", "max"),
                MSGFScoreRange=("MSGFScore_num", lambda x: x.max() - x.min()),
            )
            .reset_index()
        )
        summary = summary.merge(msgf, on=group_cols, how="left")

    return summary


def normalize_partition_scan_core_rows(long_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reduce post-SYN annotation expansions to one observation per
    RunFolder + ScanNum + CorePeptide.

    PartitionFDR should be identical within each group because alternative SYN rows
    inherit the FDR of the scored representative. SpecEValue is checked when
    available. Conflicting groups are returned for audit and cause a hard stop
    in main(). Representative metadata is selected deterministically by lowest
    SpecEValue, then lowest PartitionFDR, then exact peptide string.
    """
    work = long_df.copy()
    work["CorePeptide"] = work["Peptide"].map(normalize_core_peptide)
    work["UnmodifiedPeptide"] = work["Peptide"].map(
        lambda x: normalize_core_peptide(x, strip_mods=True)
    )
    work = work[work["CorePeptide"] != ""].copy()

    partition_group_cols = ["RunFolder", "ScanNum", "CorePeptide"]
    work["PartitionFDR_num"] = pd.to_numeric(work["PartitionFDR_num"], errors="coerce")
    if "SpecEValue" in work.columns:
        work["SpecEValue_num_for_check"] = pd.to_numeric(
            work["SpecEValue"], errors="coerce"
        )
    else:
        work["SpecEValue_num_for_check"] = pd.NA

    conflict_summary = (
        work.groupby(partition_group_cols, dropna=False)
        .agg(
            ExactPeptideForms=("Peptide", "nunique"),
            ProteinMappings=("Protein", "nunique"),
            PartitionFDRUnique=("PartitionFDR_num", lambda x: x.dropna().nunique()),
            PartitionFDRMin=("PartitionFDR_num", "min"),
            PartitionFDRMax=("PartitionFDR_num", "max"),
            SpecEValueUnique=(
                "SpecEValue_num_for_check",
                lambda x: pd.to_numeric(x, errors="coerce").dropna().nunique(),
            ),
            SpecEValueMin=("SpecEValue_num_for_check", "min"),
            SpecEValueMax=("SpecEValue_num_for_check", "max"),
        )
        .reset_index()
    )
    conflicts = conflict_summary[
        (conflict_summary["PartitionFDRUnique"] > 1)
        | (conflict_summary["SpecEValueUnique"] > 1)
    ].copy()

    work["SpecEValue_num_for_sort"] = pd.to_numeric(
        work["SpecEValue_num_for_check"], errors="coerce"
    )
    sort_cols = partition_group_cols + [
        "SpecEValue_num_for_sort", "PartitionFDR_num", "Peptide", "Protein"
    ]
    work = work.sort_values(sort_cols, kind="mergesort", na_position="last")
    representatives = (
        work.groupby(partition_group_cols, dropna=False, as_index=False)
        .head(1)
        .copy()
    )

    metadata = (
        work.groupby(partition_group_cols, dropna=False)
        .agg(
            PartitionExactPeptideForms=("Peptide", lambda x: ";".join(sorted(set(map(str, x))))),
            PartitionExactPeptideFormCount=("Peptide", "nunique"),
            PartitionProteinMappingCount=("Protein", "nunique"),
        )
        .reset_index()
    )
    representatives = representatives.merge(
        metadata, on=partition_group_cols, how="left"
    )
    representatives = representatives.drop(
        columns=["SpecEValue_num_for_check", "SpecEValue_num_for_sort"],
        errors="ignore",
    )
    return representatives, conflicts


def build_core_summary(partition_core_df: pd.DataFrame) -> pd.DataFrame:
    """Summarize one normalized partition observation per ScanNum + CorePeptide."""
    group_cols = ["ScanNum", "CorePeptide"]
    summary = summarize_fdr(partition_core_df, group_cols, prefix="Binned")

    extra = (
        partition_core_df.groupby(group_cols, dropna=False)
        .agg(
            UnmodifiedPeptide=(
                "UnmodifiedPeptide",
                lambda x: sorted({str(v) for v in x if str(v)})[0]
                if any(str(v) for v in x) else "",
            ),
            BinnedPeptideForms=(
                "PartitionExactPeptideForms",
                lambda x: ";".join(
                    sorted({
                        form
                        for value in x
                        for form in str(value).split(";")
                        if form
                    })
                ),
            ),
            BinnedPeptideFormCount=(
                "PartitionExactPeptideForms",
                lambda x: len({
                    form
                    for value in x
                    for form in str(value).split(";")
                    if form
                }),
            ),
            BinnedProteinCount=("Protein", "nunique"),
            BinnedProteins=(
                "Protein", lambda x: ";".join(sorted(set(map(str, x))))
            ),
            MaxExactFormsWithinPartition=("PartitionExactPeptideFormCount", "max"),
            PartitionsWithMultipleExactForms=(
                "PartitionExactPeptideFormCount", lambda x: int((x > 1).sum())
            ),
        )
        .reset_index()
    )
    return summary.merge(extra, on=group_cols, how="left")


def collect_partition_files(dataset_dir: Path, fdr_dir_name: str, pattern: str, partition_type: str) -> pd.DataFrame:
    dfs: list[pd.DataFrame] = []

    folder_glob = "*_pair_*" if partition_type == "pair" else "*_group_*"
    partition_folders = sorted(
        p for p in dataset_dir.glob(folder_glob)
        if p.is_dir() and parse_partition(p.name, partition_type) is not None
    )

    if not partition_folders:
        raise SystemExit(f"ERROR: no partition folders found in {dataset_dir}")

    for run_folder in partition_folders:
        partition = parse_partition(run_folder.name, partition_type)
        assert partition is not None

        fdr_dir = run_folder / fdr_dir_name
        if not fdr_dir.is_dir():
            print(f"WARNING: missing {fdr_dir}, skipping")
            continue

        files = sorted(p for p in fdr_dir.glob(pattern) if not p.name.startswith("._"))
        if not files:
            print(f"WARNING: no files matching {pattern} in {fdr_dir}, skipping")
            continue

        for path in files:
            df = pd.read_csv(path, sep="\t")
            require_columns(df, {"ScanNum", "Peptide", "Protein", "FDR"}, path)

            # Rename the input FDR immediately. In this workflow, FDR is a
            # partition-specific estimate, not a final global/collapsed value.
            df["PartitionFDR"] = pd.to_numeric(df["FDR"], errors="coerce")
            df["PartitionFDR_num"] = df["PartitionFDR"]
            df = df.drop(columns=["FDR"], errors="ignore")

            # Preserve a stable, familiar column order while keeping useful extras if present.
            preferred_cols = [
                "ScanNum", "Peptide", "Protein", "absPPM",
                "MSGFScore", "SpecEValue", "EValue", "QValue",
                "MSMSScore", "PepQValue", "ElutionTime",
                "ParentIonIntensity", "PeakArea", "StatMomentsArea",
                "IsDecoy", "PartitionFDR",
            ]
            keep_cols = ordered_existing_columns(df, preferred_cols)
            extra_cols = [c for c in df.columns if c not in keep_cols]
            df = df[keep_cols + extra_cols].copy()

            df["RunFolder"] = run_folder.name
            for key, value in partition.items():
                df[key] = value
            df["SourceFile"] = path.name

            df["Peptide"] = df["Peptide"].astype(str).str.strip()
            df["Protein"] = df["Protein"].astype(str).str.strip()
            df["CorePeptide"] = df["Peptide"].map(normalize_core_peptide)
            df["UnmodifiedPeptide"] = df["Peptide"].map(
                lambda x: normalize_core_peptide(x, strip_mods=True)
            )

            if "IsDecoy" in df.columns:
                df["IsDecoy"] = normalize_bool_series(df["IsDecoy"])

            dfs.append(df)

    if not dfs:
        raise SystemExit("No usable partition FDR files found.")

    return pd.concat(dfs, ignore_index=True, sort=False)


def build_annotated_table(
    long_df: pd.DataFrame,
    summary: pd.DataFrame,
    group_cols: list[str],
    collapse_to_representatives: bool,
) -> pd.DataFrame:
    if collapse_to_representatives:
        base = choose_representative_rows(long_df, group_cols=group_cols)
    else:
        base = long_df.copy()

    annotated = base.merge(summary, on=group_cols, how="left")

    # Put familiar columns first. There is intentionally no generic FDR column.
    front = [
        "ScanNum", "Peptide", "Protein", "absPPM",
        "MSGFScore", "SpecEValue", "EValue", "QValue",
        "MSMSScore", "PepQValue", "ElutionTime",
        "ParentIonIntensity", "PeakArea", "StatMomentsArea",
        "IsDecoy", "RunFolder", "PartitionType", "PartitionID", "PairGroupA", "PairGroupB", "Group", "SourceFile", "PartitionFDR",
        "BinnedFDRMean", "BinnedFDRMedian", "BinnedFDRMin", "BinnedFDRMax",
        "BinnedFDRStd", "BinnedFDRRange", "BinnedFDRCount",
        "PartitionBinsRepresented", "BinnedRunFolders", "RowsAcrossBins",
        "MSGFScoreMin", "MSGFScoreMax", "MSGFScoreRange",
    ]
    ordered = ordered_existing_columns(annotated, front)
    rest = [c for c in annotated.columns if c not in ordered and c != "PartitionFDR_num"]
    return annotated[ordered + rest]


def write_table(df: pd.DataFrame, path: Path) -> None:
    # Keep internal numeric helper columns out of user-facing output.
    out = df.drop(columns=["PartitionFDR_num"], errors="ignore")
    out.to_csv(path, sep="\t", index=False)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Collect partition FDR observations and write explicit binned FDR summaries. "
            "No primary FDR statistic is chosen; downstream steps must explicitly select "
            "BinnedFDRMean, BinnedFDRMedian, or another summary column."
        )
    )
    ap.add_argument("--dataset_dir", required=True, type=Path)
    ap.add_argument("--partition_type", required=True, choices=["pair", "group"])
    ap.add_argument("--fdr_dir_name", default="fdr_esti")
    ap.add_argument("--pattern", default="*_withsyn.tsv")
    ap.add_argument("--downstream_dir", default=None)

    ap.add_argument("--long_out", default=None)
    ap.add_argument("--annotated_out", default=None)
    ap.add_argument(
        "--scan_peptide_summary_out",
        default="scan_peptide_fdr_summary.tsv",
        help="Legacy exact-flanked ScanNum+Peptide diagnostic summary.",
    )
    ap.add_argument(
        "--scan_core_summary_out",
        default="scan_core_peptide_fdr_summary.tsv",
        help="Primary one-row-per-ScanNum+CorePeptide summary for comparison.",
    )
    ap.add_argument(
        "--partition_core_observations_out",
        default=None,
        help="One normalized observation per partition+ScanNum+CorePeptide.",
    )
    ap.add_argument(
        "--partition_core_conflicts_out",
        default=None,
        help="Audit table for within-partition FDR/SpecEValue conflicts.",
    )
    ap.add_argument("--triad_summary_out", default="scan_peptide_protein_fdr_summary.tsv")

    ap.add_argument(
        "--collapse_annotated_to_representatives",
        action="store_true",
        help=(
            "Write one representative metadata row per ScanNum+Peptide+Protein in the annotated table. "
            "Binned FDR summaries remain explicit; no primary FDR value is chosen."
        ),
    )
    args = ap.parse_args()

    dataset_dir = args.dataset_dir.expanduser().resolve()
    prefix = args.partition_type
    downstream_name = args.downstream_dir or ("Pairwise_Downstream" if prefix == "pair" else "Group_Downstream")
    long_out = args.long_out or f"{prefix}_fdr_long.tsv"
    annotated_out = args.annotated_out or f"{prefix}_fdr_annotated.tsv"
    observations_out = args.partition_core_observations_out or f"{prefix}_scan_core_fdr_observations.tsv"
    conflicts_out = args.partition_core_conflicts_out or f"{prefix}_scan_core_conflicts.tsv"
    downstream_dir = dataset_dir / downstream_name
    downstream_dir.mkdir(parents=True, exist_ok=True)

    print("Collecting partition-level FDR rows...")
    long_df = collect_partition_files(
        dataset_dir=dataset_dir,
        fdr_dir_name=args.fdr_dir_name,
        pattern=args.pattern,
        partition_type=args.partition_type,
    )

    long_path = downstream_dir / long_out
    write_table(long_df, long_path)

    print("Normalizing to one observation per partition + ScanNum + CorePeptide...")
    partition_core_df, partition_core_conflicts = normalize_partition_scan_core_rows(long_df)
    conflict_path = downstream_dir / conflicts_out
    write_table(partition_core_conflicts, conflict_path)
    if not partition_core_conflicts.empty:
        preview = partition_core_conflicts.head(20).to_string(index=False)
        raise SystemExit(
            "ERROR: conflicting PartitionFDR or SpecEValue values were found within "
            "partition+ScanNum+CorePeptide groups. No core summary was written. "
            f"See {conflict_path}. Example rows:\n{preview}"
        )

    partition_core_path = downstream_dir / observations_out
    write_table(partition_core_df, partition_core_path)

    print("Summarizing by ScanNum + CorePeptide...")
    scan_core_summary = build_core_summary(partition_core_df)
    core_path = downstream_dir / args.scan_core_summary_out
    write_table(scan_core_summary, core_path)

    print("Summarizing by ScanNum + exact Peptide (diagnostic)...")
    scan_peptide_summary = summarize_fdr(long_df, ["ScanNum", "Peptide"], prefix="Binned")
    sp_path = downstream_dir / args.scan_peptide_summary_out
    write_table(scan_peptide_summary, sp_path)

    print("Summarizing by ScanNum + Peptide + Protein...")
    triad_summary = summarize_fdr(long_df, ["ScanNum", "Peptide", "Protein"], prefix="Binned")
    triad_path = downstream_dir / args.triad_summary_out
    write_table(triad_summary, triad_path)

    print("Building annotated partition-level table...")
    annotated = build_annotated_table(
        long_df=long_df,
        summary=triad_summary,
        group_cols=["ScanNum", "Peptide", "Protein"],
        collapse_to_representatives=args.collapse_annotated_to_representatives,
    )
    annotated_path = downstream_dir / annotated_out
    write_table(annotated, annotated_path)

    print()
    print(f"Wrote long table:                 {long_path}")
    print(f"Wrote annotated table:            {annotated_path}")
    print(f"Wrote normalized partition/core rows:  {partition_core_path}")
    print(f"Wrote partition/core conflict audit:   {conflict_path}")
    print(f"Wrote scan/core summary:          {core_path}")
    print(f"Wrote exact scan/peptide summary: {sp_path}")
    print(f"Wrote scan/peptide/protein table: {triad_path}")
    print()
    print(f"Input rows:             {len(long_df):,}")
    print(f"Normalized partition/core observations: {len(partition_core_df):,}")
    print(f"Unique ScanNum/CorePeptide:         {scan_core_summary.shape[0]:,}")
    print(f"Unique ScanNum/exact Peptide:       {scan_peptide_summary.shape[0]:,}")
    print(f"Unique triads:                      {triad_summary.shape[0]:,}")
    print(f"Unique scans:           {long_df['ScanNum'].nunique():,}")
    print()
    print("PartitionBinsRepresented distribution, triad-level:")
    print(triad_summary["PartitionBinsRepresented"].value_counts().sort_index().to_string())
    print()
    print("Binned FDR summary, triad-level:")
    cols = [
        "BinnedFDRMean", "BinnedFDRMedian", "BinnedFDRMin",
        "BinnedFDRMax", "BinnedFDRStd", "BinnedFDRRange",
    ]
    print(triad_summary[cols].describe().to_string())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
