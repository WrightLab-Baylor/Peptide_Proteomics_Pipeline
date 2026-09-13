#!/usr/bin/env python3
"""Unified classical-guided SIC culling for Individual-bin and Pairwise-bin analyses.

This is the public culling interface used by the Binning pipeline. It dispatches
to the required implementation engines under ``scripts/binning/_internal`` and
writes dataset-level reports with explicit group/pair filename prefixes so the
two analysis types cannot overwrite one another.

Public terminology uses Individual bins and Pairwise bins. The command-line
values ``group`` and ``pair`` and the corresponding on-disk directory names are
retained for compatibility with the validated workflow.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Optional

import pandas as pd


def load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import engine: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_reports(
    report_dir: Path,
    prefix: str,
    reference_rows: pd.DataFrame,
    scan_summary: pd.DataFrame,
    summary: pd.DataFrame,
    audit: pd.DataFrame,
    recovery_by_scan: pd.DataFrame,
    unrecovered: pd.DataFrame,
    recovery_summary: pd.DataFrame,
    dry_run: bool,
) -> dict[str, Path]:
    paths = {
        "reference": report_dir / f"classical_guided_{prefix}_reference_scan_peptides.tsv",
        "reference_summary": report_dir / f"classical_guided_{prefix}_reference_scan_summary.tsv",
        "culling_summary": report_dir / f"classical_guided_{prefix}_culling_summary.tsv",
        "scan_audit": report_dir / f"classical_guided_{prefix}_scan_audit.tsv",
        "recovery_by_scan": report_dir / f"classical_guided_{prefix}_reference_recovery_by_scan.tsv",
        "unrecovered": report_dir / f"classical_guided_{prefix}_reference_unrecovered_scans.tsv",
        "recovery_summary": report_dir / f"classical_guided_{prefix}_reference_recovery_summary.tsv",
    }
    if not dry_run:
        report_dir.mkdir(parents=True, exist_ok=True)
        reference_rows.to_csv(paths["reference"], sep="\t", index=False)
        scan_summary.to_csv(paths["reference_summary"], sep="\t", index=False)
        summary.to_csv(paths["culling_summary"], sep="\t", index=False)
        audit.to_csv(paths["scan_audit"], sep="\t", index=False)
        recovery_by_scan.to_csv(paths["recovery_by_scan"], sep="\t", index=False)
        unrecovered.to_csv(paths["unrecovered"], sep="\t", index=False)
        recovery_summary.to_csv(paths["recovery_summary"], sep="\t", index=False)
    return paths


def run_engine(args: argparse.Namespace) -> int:
    dataset_dir = args.dataset_dir.expanduser().resolve()
    reference_withsyn = args.reference_withsyn.expanduser().resolve()
    tools_dir = Path(__file__).resolve().parent
    prefix = args.output_prefix or args.partition_type

    if args.partition_type == "pair":
        engine_path = (args.pair_engine or tools_dir / "_internal/classical_guided_cull_pairwise.py").resolve()
        module = load_module(engine_path, "pair_cull_engine")
        allowed, reference_rows, scan_summary, ref_scan_to_peptides = module.load_reference_withsyn(
            reference_withsyn, args.scan_col, args.peptide_col, args.decoy_col,
            args.protein_col, args.decoy_prefix,
        )
        result = module.cull_pair_files(
            dataset_dir=dataset_dir,
            allowed_pairs=allowed,
            reference_pairs=reference_rows,
            reference_scan_to_peptides=ref_scan_to_peptides,
            sic_dir_name=args.sic_dir_name,
            sic_glob=args.sic_glob,
            syn_dir_name=args.syn_dir_name,
            syn_glob=args.syn_glob,
            use_syn_fallback=not args.no_syn_fallback,
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
        singular, plural = "Pairwise-bin", "Pairwise bins"
    else:
        engine_path = (args.group_engine or tools_dir / "_internal/classical_guided_cull_individual.py").resolve()
        module = load_module(engine_path, "group_cull_engine")
        allowed, reference_rows, scan_summary, ref_scan_to_peptides = module.load_reference_withsyn(
            reference_withsyn, args.scan_col, args.peptide_col, args.decoy_col,
            args.protein_col, args.decoy_prefix,
        )
        result = module.cull_group_files(
            dataset_dir=dataset_dir,
            allowed_groups=allowed,
            reference_groups=reference_rows,
            reference_scan_to_peptides=ref_scan_to_peptides,
            sic_dir_name=args.sic_dir_name,
            sic_glob=args.sic_glob,
            syn_dir_name=args.syn_dir_name,
            syn_glob=args.syn_glob,
            use_syn_fallback=not args.no_syn_fallback,
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
        singular, plural = "Individual-bin", "Individual bins"

    summary, audit, recovery_by_scan, unrecovered, recovery_summary = result
    report_dir = (args.report_dir or dataset_dir / f"{args.partition_type.capitalize()}_Downstream" / "classical_guided_culling").expanduser().resolve()
    paths = write_reports(
        report_dir, prefix, reference_rows, scan_summary, summary, audit,
        recovery_by_scan, unrecovered, recovery_summary, args.dry_run,
    )

    print(f"Classical-guided {singular} culling")
    print("=" * 64)
    print(f"Dataset:          {dataset_dir}")
    print(f"Reference:        {reference_withsyn}")
    print(f"Partition type:   {args.partition_type}")
    print(f"Output prefix:    {prefix}")
    print(f"Engine:           {engine_path}")
    print(f"Prepared SICs:    {args.sic_dir_name}/{args.sic_glob}")
    print(f"Culled directory: {args.out_dir_name}")
    print(f"Report directory: {report_dir}")
    print(f"Dry/audit run:    {args.dry_run}")
    print()
    print(f"Reference scan-peptide keys: {len(reference_rows):,}")
    print(f"Partition files represented: {len(summary):,}")
    if not summary.empty and "Status" in summary:
        print(f"Statuses: {summary['Status'].value_counts().to_dict()}")
    metrics = dict(zip(recovery_summary.get("Metric", []), recovery_summary.get("Value", [])))
    print(f"Recovered fraction: {metrics.get('RecoveredFraction', 'NA')}")
    print(f"Reports {'would be written' if args.dry_run else 'written'}:")
    for path in paths.values():
        print(f"  {path}")
    print("=" * 64)
    return 0


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Cull classical-guided SIC outputs for Individual-bin (group) or Pairwise-bin (pair) partitions.")
    p.add_argument("--partition-type", choices=["pair", "group"], required=True)
    p.add_argument("--dataset-dir", "--dataset_dir", dest="dataset_dir", type=Path, required=True)
    p.add_argument("--reference-withsyn", "--reference_withsyn", dest="reference_withsyn", type=Path, required=True)
    p.add_argument("--output-prefix", default=None, help="Report filename prefix; default is partition type.")
    p.add_argument("--report-dir", type=Path, default=None, help="Directory for culling reports. Default: <dataset>/<Pairwise|Group>_Downstream/classical_guided_culling (legacy on-disk names retained for compatibility).")
    p.add_argument("--pair-engine", type=Path)
    p.add_argument("--group-engine", type=Path)
    p.add_argument("--sic-dir-name", "--sic_dir_name", dest="sic_dir_name", default="SICs")
    p.add_argument("--sic-glob", "--sic_glob", dest="sic_glob", default="*_PlusSICStats.tsv")
    p.add_argument("--syn-dir-name", "--syn_dir_name", dest="syn_dir_name", default="results/PHRPOut")
    p.add_argument("--syn-glob", "--syn_glob", dest="syn_glob", default="*_syn.txt")
    p.add_argument("--out-dir-name", "--out_dir_name", dest="out_dir_name", default=None)
    p.add_argument("--out-suffix", "--out_suffix", dest="out_suffix", default="_classical_guided.tsv")
    p.add_argument("--scan-col", default="ScanNum")
    p.add_argument("--peptide-col", default="Peptide")
    p.add_argument("--protein-col", default="Protein")
    p.add_argument("--decoy-col", default="IsDecoy")
    p.add_argument("--decoy-prefix", default="XXX_")
    p.add_argument("--no-syn-fallback", action="store_true")
    p.add_argument("--write-source-layer", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="Compute an in-memory audit but write no outputs.")
    args = p.parse_args(argv)
    if args.out_dir_name is None:
        args.out_dir_name = f"SICs_classical_guided_{args.partition_type}"
    return args


if __name__ == "__main__":
    try:
        raise SystemExit(run_engine(parse_args()))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
