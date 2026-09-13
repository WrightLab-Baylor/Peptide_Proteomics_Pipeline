#!/usr/bin/env python3
"""
Internal Individual-bin engine for classical-guided SIC culling.

This module is used by ``scripts/binning/classical_guided_cull.py`` when the
Individual-bin workflow is selected. It is not intended to be the primary
user-facing entry point.

For each Individual-bin folder, the engine retains SIC rows whose
``(ScanNum, flank-stripped peptide)`` assignment matches the classical/full-
database ``*_withsyn.tsv`` reference. Peptide modification notation is retained
during matching. When enabled, a matching PHRP SYN row can provide a fallback
for a reference scan that was not recovered from SICs in that bin.

Key behavior
------------
- SIC-derived rows take precedence over SYN fallback rows.
- SYN fallback contributes at most one representative row per scan per bin.
- Only reference-matching scan/peptide assignments are retained.
- Fallback rows must contain a numeric ``SpecEValue`` before FDR estimation.
- Internal helper/provenance columns are omitted from final culled outputs by
  default.
- Original input files are never modified.

Primary outputs
---------------
Each Individual-bin folder receives a classical-guided SIC table under the
configured output directory. Dataset-level audit and recovery tables are also
written, including culling summaries, scan audits, unrecovered-reference scans,
and overall reference-recovery summaries.

Use ``classical_guided_cull.py --help`` for the supported public interface.
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd


GROUP_RE = re.compile(r".*_group_(\d+)$")


def parse_group_number(run_name: str) -> Optional[int]:
    match = GROUP_RE.match(run_name)
    if not match:
        return None
    return int(match.group(1))


def normalize_core_peptide(value: object) -> str:
    """
    Remove MS-GF+/PHRP flanking residues but retain modification notation.

    Examples:
      -.DCFNGYCYGCC.R   -> DCFNGYCYGCC
      F.YANYEEYFGY.-    -> YANYEEYFGY
      A.HPAPPM*R.-      -> HPAPPM*R

    This is the identity used for culling. PTMs are intentionally retained.
    """
    if pd.isna(value):
        return ""
    pep = str(value).strip()
    if not pep:
        return ""

    if pep.count(".") >= 2:
        first = pep.find(".")
        last = pep.rfind(".")
        if first < last:
            pep = pep[first + 1:last]

    return pep.strip()


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def is_decoy_protein(protein: object, decoy_prefix: str) -> bool:
    if pd.isna(protein):
        return False
    # Some PHRP fields can contain multiple proteins; any decoy-looking protein
    # should mark the row as decoy for row-level accounting.
    text = str(protein)
    parts = re.split(r"[;,\s]+", text)
    return any(part.startswith(decoy_prefix) for part in parts if part)


def require_columns(df: pd.DataFrame, required: Set[str], path: Path) -> None:
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"{path} is missing required columns: {sorted(missing)}")


def find_group_folders(dataset_dir: Path) -> List[Path]:
    group_folders = sorted(
        p for p in dataset_dir.glob("*_group_*")
        if p.is_dir() and parse_group_number(p.name) is not None
    )
    if not group_folders:
        raise RuntimeError(
            f"No group folders matching '*_group_i' found in {dataset_dir}"
        )
    return group_folders


def pick_first_existing_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def load_reference_withsyn(
    reference_withsyn: Path,
    scan_col: str,
    peptide_col: str,
    decoy_col: Optional[str],
    protein_col: str,
    decoy_prefix: str,
) -> Tuple[Set[Tuple[object, str]], pd.DataFrame, pd.DataFrame, Dict[object, Set[str]]]:
    if not reference_withsyn.exists():
        raise FileNotFoundError(f"Reference withsyn file not found: {reference_withsyn}")

    ref = pd.read_csv(reference_withsyn, sep="\t")
    required = {scan_col, peptide_col}
    if protein_col:
        required.add(protein_col)
    require_columns(ref, required, reference_withsyn)

    ref = ref.copy()
    ref["CorePeptide"] = ref[peptide_col].map(normalize_core_peptide)
    ref = ref[ref["CorePeptide"] != ""].copy()

    if decoy_col and decoy_col in ref.columns:
        ref["ReferenceIsDecoy"] = ref[decoy_col].map(truthy)
    elif protein_col and protein_col in ref.columns:
        ref["ReferenceIsDecoy"] = ref[protein_col].map(lambda x: is_decoy_protein(x, decoy_prefix))
    else:
        ref["ReferenceIsDecoy"] = False

    reference_groups = (
        ref[[scan_col, peptide_col, protein_col, "CorePeptide", "ReferenceIsDecoy"]]
        .drop_duplicates()
        .rename(columns={scan_col: "ScanNum", peptide_col: "ReferencePeptide", protein_col: "ReferenceProtein"})
        .sort_values(["ScanNum", "CorePeptide", "ReferencePeptide", "ReferenceProtein"], kind="mergesort")
        .reset_index(drop=True)
    )

    # Collapse to the scan/core identity for allowed matching, while preserving
    # decoy status if any row for that identity is decoy.
    reference_groups = (
        reference_groups.groupby(["ScanNum", "CorePeptide"], as_index=False)
        .agg(
            ReferencePeptide=("ReferencePeptide", lambda x: ";".join(sorted(map(str, set(x))))),
            ReferenceProtein=("ReferenceProtein", lambda x: ";".join(sorted(map(str, set(x))))),
            ReferenceIsDecoy=("ReferenceIsDecoy", "any"),
        )
        .sort_values(["ScanNum", "CorePeptide"], kind="mergesort")
        .reset_index(drop=True)
    )

    allowed_groups = set(zip(reference_groups["ScanNum"], reference_groups["CorePeptide"]))

    reference_scan_to_peptides: Dict[object, Set[str]] = defaultdict(set)
    for scan, pep in zip(reference_groups["ScanNum"], reference_groups["CorePeptide"]):
        reference_scan_to_peptides[scan].add(pep)

    scan_summary = (
        reference_groups.groupby("ScanNum", as_index=False)
        .agg(
            ReferenceCorePeptideCount=("CorePeptide", "nunique"),
            ReferenceHasDecoy=("ReferenceIsDecoy", "any"),
        )
        .sort_values("ScanNum", kind="mergesort")
        .reset_index(drop=True)
    )

    return allowed_groups, reference_groups, scan_summary, reference_scan_to_peptides


def load_sic_file(path: Path, scan_col: str, peptide_col: str, protein_col: str) -> pd.DataFrame:
    """
    Load a PlusSICStats file and normalize common column-name variants.

    Older/newer prepare_SICs outputs may use Scan instead of ScanNum, while
    peptide and protein columns can also vary slightly. Normalize the detected
    columns onto the names requested by the downstream culler so the remaining
    logic uses one stable schema.
    """
    df = pd.read_csv(path, sep="\t")

    actual_scan_col = scan_col if scan_col in df.columns else pick_first_existing_column(
        df,
        ["ScanNum", "Scan", "Scan_Number", "ScanNumber", "Scan_Number_Start"],
    )
    actual_peptide_col = peptide_col if peptide_col in df.columns else pick_first_existing_column(
        df,
        ["Peptide", "PeptideSequence", "Sequence", "Peptide_Sequence"],
    )
    actual_protein_col = protein_col if protein_col in df.columns else pick_first_existing_column(
        df,
        ["Protein", "ProteinName", "Protein_Name", "ProteinID", "Protein_ID"],
    )

    missing = []
    if actual_scan_col is None:
        missing.append(scan_col)
    if actual_peptide_col is None:
        missing.append(peptide_col)
    if actual_protein_col is None:
        missing.append(protein_col)

    if missing:
        raise RuntimeError(
            f"{path} is missing required SIC fields. Could not identify: "
            f"{missing}. Available columns: {list(df.columns)}"
        )

    df = df.copy()

    if actual_scan_col != scan_col:
        df[scan_col] = df[actual_scan_col]
    if actual_peptide_col != peptide_col:
        df[peptide_col] = df[actual_peptide_col]
    if actual_protein_col != protein_col:
        df[protein_col] = df[actual_protein_col]

    df["CorePeptide"] = df[peptide_col].map(normalize_core_peptide)
    df["SourceLayer"] = "SIC"
    df["SourcePath"] = str(path)
    return df


def find_syn_file(group_folder: Path, syn_dir_name: str, syn_glob: str) -> Optional[Path]:
    syn_dir = group_folder / syn_dir_name
    if not syn_dir.is_dir():
        return None
    files = sorted(p for p in syn_dir.glob(syn_glob) if not p.name.startswith("._"))
    return files[0] if files else None


def load_syn_file(
    path: Path,
    scan_col: str,
    peptide_col: str,
    protein_col: str,
    decoy_prefix: str,
) -> pd.DataFrame:
    """
    Load a PHRP SYN file and standardize the minimum fields needed for fallback.

    The group-level SYN file may not have exactly the same columns as PlusSICStats.
    We keep all SYN columns in memory, then later project onto the SIC output
    schema so FDR estimation receives a consistent table.
    """
    df = pd.read_csv(path, sep="\t")

    actual_scan_col = scan_col if scan_col in df.columns else pick_first_existing_column(
        df, ["ScanNum", "Scan", "Scan_Number", "ScanNumber"]
    )
    actual_peptide_col = peptide_col if peptide_col in df.columns else pick_first_existing_column(
        df, ["Peptide", "PeptideSequence", "Sequence"]
    )
    actual_protein_col = protein_col if protein_col in df.columns else pick_first_existing_column(
        df, ["Protein", "ProteinName", "Protein_Name", "ProteinID", "Protein_ID"]
    )

    missing = []
    if actual_scan_col is None:
        missing.append(scan_col)
    if actual_peptide_col is None:
        missing.append(peptide_col)
    if missing:
        raise RuntimeError(
            f"{path} is missing required SYN columns. Could not identify: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    df = df.copy()
    if actual_scan_col != scan_col:
        df[scan_col] = df[actual_scan_col]
    if actual_peptide_col != peptide_col:
        df[peptide_col] = df[actual_peptide_col]
    if actual_protein_col is None:
        df[protein_col] = ""
    elif actual_protein_col != protein_col:
        df[protein_col] = df[actual_protein_col]

    # PHRP SYN files retain the MS-GF+ spectral score under the original
    # MSGFDB_SpecEValue name, while PlusSICStats and the downstream FDR
    # estimator expect SpecEValue. Preserve the scored PSM value explicitly
    # so SYN fallback rows participate in exactly the same ranking field as
    # genuine SIC-derived rows.
    if "SpecEValue" not in df.columns and "MSGFDB_SpecEValue" in df.columns:
        df["SpecEValue"] = df["MSGFDB_SpecEValue"]

    df["CorePeptide"] = df[peptide_col].map(normalize_core_peptide)
    df["SourceLayer"] = "SYN"
    df["SourcePath"] = str(path)
    df["IsDecoy"] = df[protein_col].map(lambda x: is_decoy_protein(x, decoy_prefix))
    return df


def project_syn_rows_to_sic_schema(
    syn_rows: pd.DataFrame,
    sic_columns: Sequence[str],
    scan_col: str,
    peptide_col: str,
    protein_col: str,
) -> pd.DataFrame:
    """Return SYN fallback rows with the same user-facing columns as SIC output."""
    if syn_rows.empty:
        return pd.DataFrame(columns=list(sic_columns))

    out = pd.DataFrame(index=syn_rows.index)
    for col in sic_columns:
        if col in syn_rows.columns:
            out[col] = syn_rows[col]
        else:
            out[col] = pd.NA

    # These should always exist after load_syn_file, but force them just in case.
    out[scan_col] = syn_rows[scan_col]
    out[peptide_col] = syn_rows[peptide_col]
    out[protein_col] = syn_rows[protein_col]

    if "SourceLayer" in sic_columns:
        out["SourceLayer"] = "SYN"
    if "SourcePath" in sic_columns:
        out["SourcePath"] = syn_rows["SourcePath"]

    return out


def choose_one_syn_fallback_per_scan(
    syn_matches: pd.DataFrame,
    already_kept_scans: Set[object],
    allowed_groups: Set[Tuple[object, str]],
    scan_col: str,
    peptide_col: str,
) -> pd.DataFrame:
    """
    From matching SYN candidates, choose at most one fallback row per ScanNum.

    We only fill scans that are absent from retained SIC rows in the same group.
    If multiple SYN rows match the same reference scan/peptide, choose a stable
    representative without score-based biological interpretation. If common
    score columns exist, sort by them only as deterministic tie-breakers.
    """
    if syn_matches.empty:
        return syn_matches.copy()

    work = syn_matches[~syn_matches[scan_col].isin(already_kept_scans)].copy()
    if work.empty:
        return work

    # Safety: ensure only reference matching rows are eligible.
    mask = [
        (scan, pep) in allowed_groups
        for scan, pep in zip(work[scan_col], work["CorePeptide"])
    ]
    work = work.loc[mask].copy()
    if work.empty:
        return work

    # Deterministic ordering. Use score columns if present, but do not use them
    # to choose between SIC and SYN; SIC has already won above.
    sort_cols: List[str] = [scan_col, "CorePeptide"]
    ascending: List[bool] = [True, True]

    for col, asc in [
        ("SpecEValue", True),
        ("MSGFDB_SpecEValue", True),
        ("EValue", True),
        ("QValue", True),
        ("PepQValue", True),
        ("MSGFScore", False),
        ("MSMSScore", False),
    ]:
        if col in work.columns:
            numeric_col = f"__sort_{col}"
            work[numeric_col] = pd.to_numeric(work[col], errors="coerce")
            sort_cols.append(numeric_col)
            ascending.append(asc)

    sort_cols.extend([peptide_col, "SourcePath"])
    ascending.extend([True, True])

    work = work.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    chosen = work.drop_duplicates(subset=[scan_col], keep="first").copy()
    chosen = chosen.drop(columns=[c for c in chosen.columns if c.startswith("__sort_")], errors="ignore")
    return chosen


def scan_category_counts(
    file_df: pd.DataFrame,
    kept_sic_df: pd.DataFrame,
    added_syn_df: pd.DataFrame,
    reference_scan_to_peptides: Dict[object, Set[str]],
    scan_col: str,
) -> Dict[str, int]:
    group_scan_to_peptides: Dict[object, Set[str]] = defaultdict(set)
    for scan, pep in zip(file_df[scan_col], file_df["CorePeptide"]):
        if pep:
            group_scan_to_peptides[scan].add(pep)

    kept_sic_scans = set(kept_sic_df[scan_col].unique()) if not kept_sic_df.empty else set()
    added_syn_scans = set(added_syn_df[scan_col].unique()) if not added_syn_df.empty else set()
    recovered_scans = kept_sic_scans | added_syn_scans

    group_scans = set(group_scan_to_peptides)
    reference_scans = set(reference_scan_to_peptides)

    scans_in_reference_and_group = group_scans & reference_scans
    alternative_only = sum(1 for scan in scans_in_reference_and_group if scan not in recovered_scans)

    return {
        "GroupUniqueScans": len(group_scans),
        "SICRecoveredScans": len(kept_sic_scans),
        "SYNFallbackScans": len(added_syn_scans),
        "RecoveredScans": len(recovered_scans),
        "AlternativeOnlyScans": alternative_only,
        "GroupScansNotInReference": len(group_scans - reference_scans),
    }


def build_recovery_outputs(
    reference_groups: pd.DataFrame,
    recovered_scan_to_groups: Dict[object, Set[str]],
    sic_recovered_scan_to_groups: Dict[object, Set[str]],
    syn_recovered_scan_to_groups: Dict[object, Set[str]],
    syn_observed_scan_to_groups: Dict[object, Set[str]],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create global recovery tables from all group-level kept rows."""
    ref_scan_to_peptides: Dict[object, Set[str]] = defaultdict(set)
    ref_scan_to_decoy: Dict[object, bool] = defaultdict(bool)

    for scan, pep, is_decoy in zip(
        reference_groups["ScanNum"],
        reference_groups["CorePeptide"],
        reference_groups["ReferenceIsDecoy"],
    ):
        ref_scan_to_peptides[scan].add(pep)
        ref_scan_to_decoy[scan] = bool(ref_scan_to_decoy[scan] or is_decoy)

    rows = []
    for scan in sorted(ref_scan_to_peptides):
        all_groups = sorted(recovered_scan_to_groups.get(scan, set()))
        sic_groups = sorted(sic_recovered_scan_to_groups.get(scan, set()))
        syn_added_groups = sorted(syn_recovered_scan_to_groups.get(scan, set()))
        syn_observed_groups = sorted(syn_observed_scan_to_groups.get(scan, set()))
        rows.append({
            "ScanNum": scan,
            "ReferenceCorePeptides": ";".join(sorted(ref_scan_to_peptides[scan])),
            "ReferenceHasDecoy": ref_scan_to_decoy[scan],
            "RecoveredInAnyGroup": len(all_groups) > 0,
            "RecoveredGroupCount": len(all_groups),
            "RecoveredGroupFolders": ";".join(all_groups),
            "SICRecoveredGroupCount": len(sic_groups),
            "SICRecoveredGroupFolders": ";".join(sic_groups),
            "SYNFallbackAddedGroupCount": len(syn_added_groups),
            "SYNFallbackAddedGroupFolders": ";".join(syn_added_groups),
            "SYNObservedMatchingGroupCount": len(syn_observed_groups),
            "SYNObservedMatchingGroupFolders": ";".join(syn_observed_groups),
            "RecoveredBySIC": len(sic_groups) > 0,
            "RecoveredBySYNFallback": len(syn_added_groups) > 0,
            "ObservedInSYNAtLeastOnce": len(syn_observed_groups) > 0,
        })

    recovery_by_scan = pd.DataFrame(rows)
    unrecovered = recovery_by_scan[~recovery_by_scan["RecoveredInAnyGroup"]].copy()

    total_reference_scans = len(recovery_by_scan)
    recovered_scans = int(recovery_by_scan["RecoveredInAnyGroup"].sum()) if total_reference_scans else 0
    sic_recovered_scans = int(recovery_by_scan["RecoveredBySIC"].sum()) if total_reference_scans else 0
    syn_fallback_recovered_scans = int(recovery_by_scan["RecoveredBySYNFallback"].sum()) if total_reference_scans else 0
    syn_observed_scans = int(recovery_by_scan["ObservedInSYNAtLeastOnce"].sum()) if total_reference_scans else 0
    unrecovered_scans = total_reference_scans - recovered_scans
    recovery_fraction = recovered_scans / total_reference_scans if total_reference_scans else float("nan")

    def distribution(col: str) -> str:
        if not total_reference_scans:
            return ""
        d = recovery_by_scan[col].value_counts().sort_index().to_dict()
        return ";".join(f"{k}:{v}" for k, v in d.items())

    summary = pd.DataFrame([
        {"Metric": "TotalReferenceScans", "Value": total_reference_scans},
        {"Metric": "RecoveredInAtLeastOneGroup", "Value": recovered_scans},
        {"Metric": "RecoveredBySICInAtLeastOneGroup", "Value": sic_recovered_scans},
        {"Metric": "RecoveredBySYNFallbackInAtLeastOneGroup", "Value": syn_fallback_recovered_scans},
        {"Metric": "ObservedMatchingSYNInAtLeastOneGroup", "Value": syn_observed_scans},
        {"Metric": "UnrecoveredInAllGroups", "Value": unrecovered_scans},
        {"Metric": "RecoveredFraction", "Value": recovery_fraction},
        {"Metric": "RecoveredGroupCountDistribution", "Value": distribution("RecoveredGroupCount")},
        {"Metric": "SICRecoveredGroupCountDistribution", "Value": distribution("SICRecoveredGroupCount")},
        {"Metric": "SYNFallbackAddedGroupCountDistribution", "Value": distribution("SYNFallbackAddedGroupCount")},
        {"Metric": "SYNObservedMatchingGroupCountDistribution", "Value": distribution("SYNObservedMatchingGroupCount")},
    ])

    return recovery_by_scan, unrecovered, summary


def write_reference_outputs(
    dataset_dir: Path,
    reference_groups: pd.DataFrame,
    scan_summary: pd.DataFrame,
    dry_run: bool,
) -> Tuple[Path, Path]:
    ref_groups_out = dataset_dir / "classical_guided_group_reference_scan_peptides.tsv"
    ref_summary_out = dataset_dir / "classical_guided_group_reference_scan_summary.tsv"

    if not dry_run:
        reference_groups.to_csv(ref_groups_out, sep="\t", index=False)
        scan_summary.to_csv(ref_summary_out, sep="\t", index=False)

    return ref_groups_out, ref_summary_out


def cull_group_files(
    dataset_dir: Path,
    allowed_groups: Set[Tuple[object, str]],
    reference_groups: pd.DataFrame,
    reference_scan_to_peptides: Dict[object, Set[str]],
    sic_dir_name: str,
    sic_glob: str,
    syn_dir_name: str,
    syn_glob: str,
    use_syn_fallback: bool,
    out_dir_name: str,
    out_suffix: str,
    decoy_prefix: str,
    scan_col: str,
    peptide_col: str,
    protein_col: str,
    overwrite: bool,
    dry_run: bool,
    write_source_layer: bool,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    group_folders = find_group_folders(dataset_dir)

    recovered_scan_to_groups: Dict[object, Set[str]] = defaultdict(set)
    sic_recovered_scan_to_groups: Dict[object, Set[str]] = defaultdict(set)
    syn_recovered_scan_to_groups: Dict[object, Set[str]] = defaultdict(set)
    syn_observed_scan_to_groups: Dict[object, Set[str]] = defaultdict(set)

    summary_rows: List[Dict[str, object]] = []
    audit_rows: List[Dict[str, object]] = []

    for group_folder in group_folders:
        group_number = parse_group_number(group_folder.name)
        assert group_number is not None

        sic_dir = group_folder / sic_dir_name
        if not sic_dir.is_dir():
            summary_rows.append({
                "RunFolder": group_folder.name,
                "Group": group_number,
                "SourceFile": "",
                "SynFile": "",
                "Status": "missing_sic_dir",
                "InputRows": 0,
                "SICRowsKept": 0,
                "SYNFallbackRowsAdded": 0,
                "OutputRows": 0,
                "TargetRows": 0,
                "DecoyRows": 0,
                "RowFDR": "",
                "GroupUniqueScans": 0,
                "SICRecoveredScans": 0,
                "SYNFallbackScans": 0,
                "RecoveredScans": 0,
                "AlternativeOnlyScans": 0,
                "GroupScansNotInReference": 0,
                "SYNMatchingRowsAvailable": 0,
                "SYNMatchingScansAvailable": 0,
                "OutputFile": "",
            })
            print(f"WARNING: Missing {sic_dir}, skipping", file=sys.stderr)
            continue

        sic_files = sorted(p for p in sic_dir.glob(sic_glob) if not p.name.startswith("._"))
        if not sic_files:
            summary_rows.append({
                "RunFolder": group_folder.name,
                "Group": group_number,
                "SourceFile": "",
                "SynFile": "",
                "Status": "missing_sic_files",
                "InputRows": 0,
                "SICRowsKept": 0,
                "SYNFallbackRowsAdded": 0,
                "OutputRows": 0,
                "TargetRows": 0,
                "DecoyRows": 0,
                "RowFDR": "",
                "GroupUniqueScans": 0,
                "SICRecoveredScans": 0,
                "SYNFallbackScans": 0,
                "RecoveredScans": 0,
                "AlternativeOnlyScans": 0,
                "GroupScansNotInReference": 0,
                "SYNMatchingRowsAvailable": 0,
                "SYNMatchingScansAvailable": 0,
                "OutputFile": "",
            })
            print(f"WARNING: No files matching {sic_glob} in {sic_dir}, skipping", file=sys.stderr)
            continue

        syn_file = find_syn_file(group_folder, syn_dir_name=syn_dir_name, syn_glob=syn_glob) if use_syn_fallback else None
        syn_df: Optional[pd.DataFrame] = None
        if use_syn_fallback and syn_file is not None:
            try:
                syn_df = load_syn_file(
                    syn_file,
                    scan_col=scan_col,
                    peptide_col=peptide_col,
                    protein_col=protein_col,
                    decoy_prefix=decoy_prefix,
                )
                syn_match_mask = [
                    (scan, pep) in allowed_groups
                    for scan, pep in zip(syn_df[scan_col], syn_df["CorePeptide"])
                ]
                syn_matching_all = syn_df.loc[syn_match_mask].copy()
                for scan in syn_matching_all[scan_col].unique():
                    syn_observed_scan_to_groups[scan].add(group_folder.name)
            except Exception as exc:
                print(f"WARNING: Failed to load SYN fallback file {syn_file}: {exc}", file=sys.stderr)
                syn_df = None
        else:
            syn_matching_all = pd.DataFrame()

        out_dir = group_folder / out_dir_name
        if not dry_run:
            out_dir.mkdir(parents=True, exist_ok=True)

        for sic_file in sic_files:
            file_df = load_sic_file(sic_file, scan_col=scan_col, peptide_col=peptide_col, protein_col=protein_col)

            sic_keep_mask = [
                (scan, pep) in allowed_groups
                for scan, pep in zip(file_df[scan_col], file_df["CorePeptide"])
            ]
            kept_sic = file_df.loc[sic_keep_mask].copy()
            kept_sic["IsDecoy"] = kept_sic[protein_col].map(lambda x: is_decoy_protein(x, decoy_prefix))

            kept_sic_scans = set(kept_sic[scan_col].unique()) if not kept_sic.empty else set()

            added_syn_projected = pd.DataFrame(columns=file_df.columns)
            chosen_syn = pd.DataFrame()
            if use_syn_fallback and syn_df is not None and not syn_df.empty:
                syn_match_mask = [
                    (scan, pep) in allowed_groups
                    for scan, pep in zip(syn_df[scan_col], syn_df["CorePeptide"])
                ]
                syn_matches = syn_df.loc[syn_match_mask].copy()
                chosen_syn = choose_one_syn_fallback_per_scan(
                    syn_matches=syn_matches,
                    already_kept_scans=kept_sic_scans,
                    allowed_groups=allowed_groups,
                    scan_col=scan_col,
                    peptide_col=peptide_col,
                )
                if not chosen_syn.empty:
                    # Match output schema to SIC output. Include internal helper columns for
                    # accounting now; they are dropped from final output below.
                    output_schema = list(file_df.columns)
                    if "IsDecoy" not in output_schema:
                        output_schema.append("IsDecoy")
                    if "SourceLayer" not in output_schema:
                        output_schema.append("SourceLayer")
                    if "SourcePath" not in output_schema:
                        output_schema.append("SourcePath")
                    added_syn_projected = project_syn_rows_to_sic_schema(
                        chosen_syn,
                        sic_columns=output_schema,
                        scan_col=scan_col,
                        peptide_col=peptide_col,
                        protein_col=protein_col,
                    )
                    added_syn_projected["CorePeptide"] = chosen_syn["CorePeptide"].values
                    added_syn_projected["IsDecoy"] = chosen_syn["IsDecoy"].values
                    added_syn_projected["SourceLayer"] = "SYN"
                    added_syn_projected["SourcePath"] = str(syn_file)

                    # A fallback row must never enter score-ranked FDR
                    # estimation without the score used by that estimator.
                    # Fail loudly instead of silently sorting an NA score to
                    # the bottom of the target-decoy ranking.
                    if "SpecEValue" not in added_syn_projected.columns:
                        raise RuntimeError(
                            f"Projected SYN fallback rows from {syn_file} lack "
                            "the required SpecEValue column."
                        )
                    fallback_scores = pd.to_numeric(
                        added_syn_projected["SpecEValue"], errors="coerce"
                    )
                    missing_score_mask = fallback_scores.isna()
                    if missing_score_mask.any():
                        bad_cols = [
                            c for c in [scan_col, peptide_col, protein_col, "SpecEValue"]
                            if c in added_syn_projected.columns
                        ]
                        bad_preview = (
                            added_syn_projected.loc[missing_score_mask, bad_cols]
                            .head(10)
                            .to_string(index=False)
                        )
                        raise RuntimeError(
                            f"{int(missing_score_mask.sum())} projected SYN fallback "
                            f"row(s) from {syn_file} lack a numeric SpecEValue after "
                            "PHRP-to-SIC column normalization. Example rows:\n"
                            f"{bad_preview}"
                        )

            if "SourceLayer" not in kept_sic.columns:
                kept_sic["SourceLayer"] = "SIC"
            if "SourcePath" not in kept_sic.columns:
                kept_sic["SourcePath"] = str(sic_file)

            # Avoid concatenating an empty/all-NA fallback frame. This preserves
            # the current result while preventing pandas FutureWarning noise and
            # future dtype changes when one side is empty.
            if added_syn_projected.empty:
                combined = kept_sic.copy()
            elif kept_sic.empty:
                combined = added_syn_projected.copy()
            else:
                # pandas currently emits a FutureWarning when one input contains
                # all-NA columns. Suppress only that exact warning for this concat;
                # do not alter either frame, since dtype coercion could affect the
                # controlled comparison.
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore",
                        message=(
                            "The behavior of DataFrame concatenation with empty or "
                            "all-NA entries is deprecated.*"
                        ),
                        category=FutureWarning,
                    )
                    combined = pd.concat(
                        [kept_sic, added_syn_projected],
                        ignore_index=True,
                        sort=False,
                    )

            # Safety check: one representative per ScanNum going to FDR estimation.
            # If SIC itself has duplicate rows for a scan, retain one stable row to
            # protect the FDR input from row multiplication. This should normally be
            # no-op for PlusSICStats.
            if not combined.empty and combined[scan_col].duplicated().any():
                source_rank = combined["SourceLayer"].map({"SIC": 0, "SYN": 1}).fillna(2)
                combined = combined.assign(__source_rank=source_rank)
                combined = combined.sort_values([scan_col, "__source_rank", "CorePeptide"], kind="mergesort")
                combined = combined.drop_duplicates(subset=[scan_col], keep="first")
                combined = combined.drop(columns=["__source_rank"], errors="ignore")

            # Update recovery maps after duplicate protection.
            for scan, source_layer in zip(combined[scan_col], combined.get("SourceLayer", pd.Series([""] * len(combined)))):
                recovered_scan_to_groups[scan].add(group_folder.name)
                if source_layer == "SYN":
                    syn_recovered_scan_to_groups[scan].add(group_folder.name)
                else:
                    sic_recovered_scan_to_groups[scan].add(group_folder.name)

            scan_counts = scan_category_counts(
                file_df=file_df,
                kept_sic_df=kept_sic,
                added_syn_df=added_syn_projected,
                reference_scan_to_peptides=reference_scan_to_peptides,
                scan_col=scan_col,
            )

            if sic_file.name.endswith("_PlusSICStats.tsv"):
                out_name = sic_file.name.replace("_PlusSICStats.tsv", out_suffix)
            else:
                out_name = sic_file.stem + out_suffix
            out_file = out_dir / out_name

            if out_file.exists() and not overwrite and not dry_run:
                status = "skipped_destination_exists"
            else:
                status = "dry_run" if dry_run else "written"
                if not dry_run:
                    final_cols_to_drop = ["CorePeptide"]
                    if not write_source_layer:
                        final_cols_to_drop.extend(["SourceLayer", "SourcePath"])
                    final_combined = combined.drop(columns=final_cols_to_drop, errors="ignore")
                    final_combined.to_csv(out_file, sep="\t", index=False)

            n_output = len(combined)
            if "IsDecoy" in combined.columns:
                n_decoy = int(combined["IsDecoy"].map(truthy).sum()) if n_output else 0
            else:
                n_decoy = int(combined[protein_col].map(lambda x: is_decoy_protein(x, decoy_prefix)).sum()) if n_output else 0
            n_target = n_output - n_decoy
            row_fdr = (n_decoy / n_output) if n_output else ""

            syn_matching_rows_available = int(len(syn_matching_all)) if use_syn_fallback and syn_file is not None else 0
            syn_matching_scans_available = int(syn_matching_all[scan_col].nunique()) if syn_matching_rows_available else 0

            summary_rows.append({
                "RunFolder": group_folder.name,
                "Group": group_number,
                "SourceFile": str(sic_file),
                "SynFile": str(syn_file) if syn_file is not None else "",
                "Status": status,
                "InputRows": len(file_df),
                "SICRowsKept": len(kept_sic),
                "SYNFallbackRowsAdded": len(added_syn_projected),
                "OutputRows": n_output,
                "TargetRows": n_target,
                "DecoyRows": n_decoy,
                "RowFDR": row_fdr,
                "GroupUniqueScans": scan_counts["GroupUniqueScans"],
                "SICRecoveredScans": scan_counts["SICRecoveredScans"],
                "SYNFallbackScans": scan_counts["SYNFallbackScans"],
                "RecoveredScans": scan_counts["RecoveredScans"],
                "AlternativeOnlyScans": scan_counts["AlternativeOnlyScans"],
                "GroupScansNotInReference": scan_counts["GroupScansNotInReference"],
                "SYNMatchingRowsAvailable": syn_matching_rows_available,
                "SYNMatchingScansAvailable": syn_matching_scans_available,
                "OutputFile": str(out_file),
            })

            group_scan_to_peptides: Dict[object, Set[str]] = defaultdict(set)
            for scan, pep in zip(file_df[scan_col], file_df["CorePeptide"]):
                if pep:
                    group_scan_to_peptides[scan].add(pep)

            kept_sic_scans_after = set(kept_sic[scan_col].unique()) if not kept_sic.empty else set()
            added_syn_scans_after = set(added_syn_projected[scan_col].unique()) if not added_syn_projected.empty else set()
            recovered_scans_after = kept_sic_scans_after | added_syn_scans_after

            for scan in sorted(group_scan_to_peptides):
                ref_peps = reference_scan_to_peptides.get(scan, set())
                group_peps = group_scan_to_peptides[scan]
                if scan in kept_sic_scans_after:
                    category = "RecoveredBySIC"
                elif scan in added_syn_scans_after:
                    category = "RecoveredBySYNFallback"
                elif ref_peps:
                    category = "AlternativeOnlyInSIC"
                else:
                    category = "GroupScanNotInReference"
                audit_rows.append({
                    "RunFolder": group_folder.name,
                    "SourceFile": sic_file.name,
                    "SynFile": str(syn_file) if syn_file is not None else "",
                    "ScanNum": scan,
                    "Category": category,
                    "GroupCorePeptides": ";".join(sorted(group_peps)),
                    "ReferenceCorePeptides": ";".join(sorted(ref_peps)),
                })

            # Add audit rows for SYN fallback scans that were absent from SIC scan list.
            for scan in sorted(added_syn_scans_after - set(group_scan_to_peptides)):
                ref_peps = reference_scan_to_peptides.get(scan, set())
                syn_peps = set(chosen_syn.loc[chosen_syn[scan_col] == scan, "CorePeptide"]) if not chosen_syn.empty else set()
                audit_rows.append({
                    "RunFolder": group_folder.name,
                    "SourceFile": sic_file.name,
                    "SynFile": str(syn_file) if syn_file is not None else "",
                    "ScanNum": scan,
                    "Category": "RecoveredBySYNFallback_AbsentFromSICScanList",
                    "GroupCorePeptides": "",
                    "ReferenceCorePeptides": ";".join(sorted(ref_peps)),
                    "SYNCorePeptides": ";".join(sorted(syn_peps)),
                })

    recovery_by_scan, unrecovered, recovery_summary = build_recovery_outputs(
        reference_groups=reference_groups,
        recovered_scan_to_groups=recovered_scan_to_groups,
        sic_recovered_scan_to_groups=sic_recovered_scan_to_groups,
        syn_recovered_scan_to_groups=syn_recovered_scan_to_groups,
        syn_observed_scan_to_groups=syn_observed_scan_to_groups,
    )

    return (
        pd.DataFrame(summary_rows),
        pd.DataFrame(audit_rows),
        recovery_by_scan,
        unrecovered,
        recovery_summary,
    )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Internal Individual-bin culling engine. Use classical_guided_cull.py as the public entry point."
    )
    p.add_argument("--dataset_dir", required=True, type=Path,
                   help="Dataset directory containing Individual-bin folders named *_group_i.")
    p.add_argument("--reference_withsyn", required=True, type=Path,
                   help="Classical/full-database *_withsyn.tsv reference file.")

    p.add_argument("--sic_dir_name", default="SICs",
                   help="Input SIC directory inside each group folder. Default: SICs")
    p.add_argument("--sic_glob", default="*_PlusSICStats.tsv",
                   help="Glob for group-level SIC files. Default: *_PlusSICStats.tsv")
    p.add_argument("--syn_dir_name", default="results/PHRPOut",
                   help="SYN directory inside each group folder. Default: results/PHRPOut")
    p.add_argument("--syn_glob", default="*_syn.txt",
                   help="Glob for group-level SYN files. Default: *_syn.txt")
    p.add_argument("--no_syn_fallback", action="store_true",
                   help="Disable SYN fallback and behave like SIC-only culling.")

    p.add_argument("--out_dir_name", default="SICs_classical_guided",
                   help="Output directory inside each group folder. Default: SICs_classical_guided")
    p.add_argument("--out_suffix", default="_classical_guided.tsv",
                   help="Output suffix. Default: _classical_guided.tsv")

    p.add_argument("--scan_col", default="ScanNum")
    p.add_argument("--peptide_col", default="Peptide")
    p.add_argument("--protein_col", default="Protein")
    p.add_argument("--decoy_col", default="IsDecoy",
                   help="Reference decoy column, if present. Default: IsDecoy. If absent, protein prefix is used.")
    p.add_argument("--decoy_prefix", default="XXX_")

    p.add_argument("--write_source_layer", action="store_true",
                   help="Write SourceLayer/SourcePath columns to output. Default: omit helper provenance columns from final output.")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing output files.")
    p.add_argument("--dry_run", action="store_true",
                   help="Report what would be done without writing group-level output files.")

    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    dataset_dir = args.dataset_dir.expanduser().resolve()
    reference_withsyn = args.reference_withsyn.expanduser().resolve()
    use_syn_fallback = not args.no_syn_fallback

    if not dataset_dir.is_dir():
        print(f"ERROR: dataset_dir not found or not a directory: {dataset_dir}", file=sys.stderr)
        return 2

    print("Classical-guided Individual-bin culling engine")
    print("=" * 60)
    print(f"Dataset directory:   {dataset_dir}")
    print(f"Reference withsyn:   {reference_withsyn}")
    print(f"Input SIC dir name:  {args.sic_dir_name}")
    print(f"SYN fallback:        {use_syn_fallback}")
    print(f"SYN dir/glob:        {args.syn_dir_name}/{args.syn_glob}")
    print(f"Output dir name:     {args.out_dir_name}")
    print(f"Dry run:             {args.dry_run}")
    print(f"Overwrite:           {args.overwrite}")
    print()

    print("Loading classical reference...")
    allowed_groups, reference_groups, scan_summary, reference_scan_to_peptides = load_reference_withsyn(
        reference_withsyn=reference_withsyn,
        scan_col=args.scan_col,
        peptide_col=args.peptide_col,
        decoy_col=args.decoy_col,
        protein_col=args.protein_col,
        decoy_prefix=args.decoy_prefix,
    )

    ref_groups_out, ref_summary_out = write_reference_outputs(
        dataset_dir=dataset_dir,
        reference_groups=reference_groups,
        scan_summary=scan_summary,
        dry_run=args.dry_run,
    )

    print(f"Reference rows loaded as unique ScanNum/CorePeptide groups: {len(reference_groups):,}")
    print(f"Reference unique scans: {scan_summary['ScanNum'].nunique():,}")
    multi = int((scan_summary["ReferenceCorePeptideCount"] > 1).sum())
    print(f"Reference scans with >1 CorePeptide: {multi:,}")
    print(f"Reference peptide table: {ref_groups_out}")
    print(f"Reference scan summary: {ref_summary_out}")
    print()

    print("Culling group SIC files, adding SYN fallback rows where needed, and auditing recovery...")
    summary, audit, recovery_by_scan, unrecovered, recovery_summary = cull_group_files(
        dataset_dir=dataset_dir,
        allowed_groups=allowed_groups,
        reference_groups=reference_groups,
        reference_scan_to_peptides=reference_scan_to_peptides,
        sic_dir_name=args.sic_dir_name,
        sic_glob=args.sic_glob,
        syn_dir_name=args.syn_dir_name,
        syn_glob=args.syn_glob,
        use_syn_fallback=use_syn_fallback,
        out_dir_name=args.out_dir_name,
        out_suffix=args.out_suffix,
        decoy_prefix=args.decoy_prefix,
        scan_col=args.scan_col,
        peptide_col=args.peptide_col,
        protein_col=args.protein_col,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        write_source_layer=args.write_source_layer,
    )

    summary_out = dataset_dir / "classical_guided_group_culling_summary.tsv"
    audit_out = dataset_dir / "classical_guided_group_scan_audit.tsv"
    recovery_by_scan_out = dataset_dir / "classical_guided_group_reference_recovery_by_scan.tsv"
    unrecovered_out = dataset_dir / "classical_guided_group_reference_unrecovered_scans.tsv"
    recovery_summary_out = dataset_dir / "classical_guided_group_reference_recovery_summary.tsv"

    if not args.dry_run:
        summary.to_csv(summary_out, sep="\t", index=False)
        audit.to_csv(audit_out, sep="\t", index=False)
        recovery_by_scan.to_csv(recovery_by_scan_out, sep="\t", index=False)
        unrecovered.to_csv(unrecovered_out, sep="\t", index=False)
        recovery_summary.to_csv(recovery_summary_out, sep="\t", index=False)

    print(f"Group summary: {summary_out}")
    print(f"Group scan audit: {audit_out}")
    print(f"Reference recovery by scan: {recovery_by_scan_out}")
    print(f"Unrecovered reference scans: {unrecovered_out}")
    print(f"Reference recovery summary: {recovery_summary_out}")
    print()

    total_input = int(summary["InputRows"].sum()) if not summary.empty else 0
    total_sic_kept = int(summary["SICRowsKept"].sum()) if not summary.empty else 0
    total_syn_added = int(summary["SYNFallbackRowsAdded"].sum()) if not summary.empty else 0
    total_output = int(summary["OutputRows"].sum()) if not summary.empty else 0
    total_target = int(summary["TargetRows"].sum()) if not summary.empty else 0
    total_decoy = int(summary["DecoyRows"].sum()) if not summary.empty else 0
    total_fdr = total_decoy / total_output if total_output else float("nan")

    status_counts = summary["Status"].value_counts().to_dict() if not summary.empty else {}

    rec_metrics = dict(zip(recovery_summary["Metric"], recovery_summary["Value"])) if not recovery_summary.empty else {}
    total_ref_scans = int(rec_metrics.get("TotalReferenceScans", 0))
    recovered_any = int(rec_metrics.get("RecoveredInAtLeastOneGroup", 0))
    recovered_by_sic = int(rec_metrics.get("RecoveredBySICInAtLeastOneGroup", 0))
    recovered_by_syn = int(rec_metrics.get("RecoveredBySYNFallbackInAtLeastOneGroup", 0))
    syn_observed = int(rec_metrics.get("ObservedMatchingSYNInAtLeastOneGroup", 0))
    unrecovered_count = int(rec_metrics.get("UnrecoveredInAllGroups", 0))
    recovered_fraction = float(rec_metrics.get("RecoveredFraction", float("nan")))
    recovery_distribution = rec_metrics.get("RecoveredGroupCountDistribution", "")
    sic_distribution = rec_metrics.get("SICRecoveredGroupCountDistribution", "")
    syn_added_distribution = rec_metrics.get("SYNFallbackAddedGroupCountDistribution", "")
    syn_observed_distribution = rec_metrics.get("SYNObservedMatchingGroupCountDistribution", "")

    print("Run summary")
    print("-" * 60)
    print(f"Group files represented: {len(summary):,}")
    print(f"Statuses: {status_counts}")
    print(f"Rows before culling:     {total_input:,}")
    print(f"SIC rows kept:           {total_sic_kept:,}")
    print(f"SYN fallback rows added: {total_syn_added:,}")
    print(f"Rows after culling:      {total_output:,}")
    print(f"Target rows kept:        {total_target:,}")
    print(f"Decoy rows kept:         {total_decoy:,}")
    print(f"Overall row FDR:         {total_fdr:.6g}" if total_output else "Overall row FDR:         NA")
    print()
    print("Reference recovery")
    print("-" * 60)
    print(f"Reference scans:                    {total_ref_scans:,}")
    print(f"Recovered in >=1 group:              {recovered_any:,}")
    print(f"Recovered by SIC in >=1 group:       {recovered_by_sic:,}")
    print(f"Recovered by SYN fallback >=1 group: {recovered_by_syn:,}")
    print(f"Observed matching SYN >=1 group:     {syn_observed:,}")
    print(f"Unrecovered in all groups:           {unrecovered_count:,}")
    print(f"Reference scan recovery fraction:   {recovered_fraction:.6g}" if total_ref_scans else "Reference scan recovery fraction:   NA")
    print(f"Recovered group count distribution:  {recovery_distribution}")
    print(f"SIC group count distribution:        {sic_distribution}")
    print(f"SYN added group count distribution:  {syn_added_distribution}")
    print(f"SYN observed group count distribution:{syn_observed_distribution}")
    print("=" * 60)
    print("Done.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
