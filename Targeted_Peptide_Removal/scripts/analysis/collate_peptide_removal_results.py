#!/usr/bin/env python3
"""
collate_peptide_removal_results.py

Collate targeted peptide-removal database size and peptide/protein gain/loss
statistics relative to each condition's unmodified baseline run.

Expected data layout
--------------------
<Targeted_Peptide_Removal>/data/
├── Condition/
│   ├── database/
│   ├── peptide_crosstab_annotated.tsv
│   └── protein_crosstab_annotated.tsv
├── Condition_1.0/
├── Condition_0.9/
├── ...
└── results/

A baseline is discovered generically as a top-level condition directory with the
required final crosstabs and at least one sibling named ``Condition_*``. Each
condition is written to ``results/Condition.tsv``. Incomplete removal runs warn
and are skipped; an incomplete baseline prevents that condition from being
summarized.

Database size is the number of FASTA records in ``database/``. Peptide counts
are distinct non-empty ``Peptide`` values and protein counts are distinct
non-empty ``Protein`` values.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Iterable

PEPTIDE_CROSSTAB = "peptide_crosstab_annotated.tsv"
PROTEIN_CROSSTAB = "protein_crosstab_annotated.tsv"
FASTA_SUFFIXES = {".fa", ".faa", ".fasta"}


def warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)


def natural_key(text: str):
    """Sort strings naturally, so 0.9 precedes 1.0 and 2 precedes 10."""
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


def find_peptide_fasta(run_dir: Path) -> Path | None:
    """
    Find the run's peptide FASTA in database/.

    If exactly one FASTA-like file exists, use it. If several exist, prefer a
    single filename containing 'peptide'. Ambiguous cases are rejected.
    """
    database_dir = run_dir / "database"
    if not database_dir.is_dir():
        warn(f"{run_dir.name}: missing database/ directory")
        return None

    candidates = sorted(
        (
            p
            for p in database_dir.iterdir()
            if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES
        ),
        key=lambda p: natural_key(p.name),
    )

    if len(candidates) == 1:
        return candidates[0]

    peptide_named = [p for p in candidates if "peptide" in p.name.lower()]
    if len(peptide_named) == 1:
        return peptide_named[0]

    if not candidates:
        warn(f"{run_dir.name}: no FASTA found in {database_dir}")
    else:
        warn(
            f"{run_dir.name}: ambiguous database FASTAs "
            f"({', '.join(p.name for p in candidates)})"
        )
    return None


def count_fasta_entries(path: Path) -> int:
    """Count FASTA records exactly as grep -c '^>' would."""
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def read_distinct_column(path: Path, column: str) -> set[str]:
    """Return distinct non-empty identifiers from a named TSV column."""
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("file is empty or lacks a header")
        if column not in reader.fieldnames:
            raise ValueError(
                f"required column {column!r} not found; available columns: {reader.fieldnames}"
            )
        entities: set[str] = set()
        for row in reader:
            entity = (row.get(column) or "").strip()
            if entity:
                entities.add(entity)
    return entities


def validate_run(run_dir: Path) -> tuple[Path, Path, Path] | None:
    """Locate all inputs required to summarize one run."""
    fasta = find_peptide_fasta(run_dir)
    peptide = run_dir / PEPTIDE_CROSSTAB
    protein = run_dir / PROTEIN_CROSSTAB

    missing = []
    if fasta is None:
        missing.append("peptide FASTA")
    if not peptide.is_file() or peptide.stat().st_size == 0:
        missing.append(PEPTIDE_CROSSTAB)
    if not protein.is_file() or protein.stat().st_size == 0:
        missing.append(PROTEIN_CROSSTAB)

    if missing:
        warn(f"{run_dir.name}: incomplete; skipping ({', '.join(missing)} missing)")
        return None

    return fasta, peptide, protein


def discover_conditions(root: Path) -> list[tuple[Path, list[Path]]]:
    """
    Discover baseline conditions and their sibling removal runs.

    A baseline is any top-level directory with both final annotated crosstabs
    and at least one sibling whose name begins with '<baseline>_'.
    """
    top_dirs = [p for p in root.iterdir() if p.is_dir()]
    discovered = []

    for baseline in top_dirs:
        if not (baseline / PEPTIDE_CROSSTAB).is_file():
            continue
        if not (baseline / PROTEIN_CROSSTAB).is_file():
            continue

        prefix = baseline.name + "_"
        variants = [
            p for p in top_dirs
            if p.name.startswith(prefix) and p != baseline
        ]
        if variants:
            discovered.append(
                (baseline, sorted(variants, key=lambda p: natural_key(p.name)))
            )

    return sorted(discovered, key=lambda item: natural_key(item[0].name))


def summarize_condition(
    baseline_dir: Path,
    variant_dirs: Iterable[Path],
    output_dir: Path,
) -> bool:
    """Write one Condition.tsv. Returns True if the condition was processed."""
    condition = baseline_dir.name
    baseline_inputs = validate_run(baseline_dir)
    if baseline_inputs is None:
        warn(
            f"{condition}: baseline is incomplete; cannot calculate "
            "condition-relative gains/losses"
        )
        return False

    baseline_fasta, baseline_peptide_file, baseline_protein_file = baseline_inputs

    try:
        baseline_db_count = count_fasta_entries(baseline_fasta)
        baseline_peptides = read_distinct_column(baseline_peptide_file, "Peptide")
        baseline_proteins = read_distinct_column(baseline_protein_file, "Protein")
    except (OSError, ValueError) as exc:
        warn(f"{condition}: failed reading baseline inputs: {exc}")
        return False

    rows = []

    # Include baseline explicitly so later graphing never needs to reopen it.
    rows.append(
        {
            "Condition": condition,
            "Run": condition,
            "Removal_Level": "baseline",
            "Is_Baseline": "True",
            "Database_Peptides": baseline_db_count,
            "Database_Peptides_Removed": 0,
            "Database_Fraction_Removed": "0.000000",
            "Peptides_Passing": len(baseline_peptides),
            "Peptides_Gained": 0,
            "Peptides_Lost": 0,
            "Net_Peptide_Change": 0,
            "Proteins_Passing": len(baseline_proteins),
            "Proteins_Gained": 0,
            "Proteins_Lost": 0,
            "Net_Protein_Change": 0,
        }
    )

    processed_variants = 0

    for run_dir in variant_dirs:
        inputs = validate_run(run_dir)
        if inputs is None:
            continue

        fasta, peptide_file, protein_file = inputs
        try:
            db_count = count_fasta_entries(fasta)
            peptides = read_distinct_column(peptide_file, "Peptide")
            proteins = read_distinct_column(protein_file, "Protein")
        except (OSError, ValueError) as exc:
            warn(f"{run_dir.name}: failed reading inputs; skipping ({exc})")
            continue

        peptide_gained = len(peptides - baseline_peptides)
        peptide_lost = len(baseline_peptides - peptides)
        protein_gained = len(proteins - baseline_proteins)
        protein_lost = len(baseline_proteins - proteins)

        db_removed = baseline_db_count - db_count
        if baseline_db_count:
            fraction_removed = db_removed / baseline_db_count
        else:
            fraction_removed = 0.0

        removal_level = run_dir.name[len(condition) + 1 :]

        rows.append(
            {
                "Condition": condition,
                "Run": run_dir.name,
                "Removal_Level": removal_level,
                "Is_Baseline": "False",
                "Database_Peptides": db_count,
                "Database_Peptides_Removed": db_removed,
                "Database_Fraction_Removed": f"{fraction_removed:.6f}",
                "Peptides_Passing": len(peptides),
                "Peptides_Gained": peptide_gained,
                "Peptides_Lost": peptide_lost,
                "Net_Peptide_Change": len(peptides) - len(baseline_peptides),
                "Proteins_Passing": len(proteins),
                "Proteins_Gained": protein_gained,
                "Proteins_Lost": protein_lost,
                "Net_Protein_Change": len(proteins) - len(baseline_proteins),
            }
        )
        processed_variants += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{condition}.tsv"

    fieldnames = [
        "Condition",
        "Run",
        "Removal_Level",
        "Is_Baseline",
        "Database_Peptides",
        "Database_Peptides_Removed",
        "Database_Fraction_Removed",
        "Peptides_Passing",
        "Peptides_Gained",
        "Peptides_Lost",
        "Net_Peptide_Change",
        "Proteins_Passing",
        "Proteins_Gained",
        "Proteins_Lost",
        "Net_Protein_Change",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[OK] {condition}: wrote {output_path} "
        f"(baseline + {processed_variants} complete removal run(s))"
    )
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collate peptide-removal database size and peptide/protein "
            "gain/loss statistics relative to baseline condition runs."
        )
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=None,
        help=(
            "Targeted peptide-removal data root containing baseline and Condition_* "
            "folders. Default: <Targeted_Peptide_Removal>/data inferred from this script."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory (default: <root>/results)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    default_root = Path(__file__).resolve().parents[2] / "data"
    root = (
        Path(args.root).expanduser().resolve()
        if args.root
        else default_root.resolve()
    )
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else root / "results"
    )

    if not root.is_dir():
        print(f"ERROR: root directory does not exist: {root}", file=sys.stderr)
        return 2

    conditions = discover_conditions(root)
    if not conditions:
        print(
            "ERROR: no baseline Condition directories with matching Condition_* "
            "runs were discovered.",
            file=sys.stderr,
        )
        return 1

    print(f"Root: {root}")
    print(f"Output: {output_dir}")
    print(
        "Discovered conditions: "
        + ", ".join(baseline.name for baseline, _ in conditions)
    )

    successes = 0
    for baseline, variants in conditions:
        if summarize_condition(baseline, variants, output_dir):
            successes += 1

    print(
        f"Done: {successes}/{len(conditions)} condition(s) summarized. "
        "Warnings above identify incomplete runs that were skipped."
    )
    return 0 if successes else 1


if __name__ == "__main__":
    raise SystemExit(main())
