#!/usr/bin/env python3
"""
conflict_culler_fast.py

Generate a peptide FASTA with all peptides participating in qualifying conflict
edges removed. Conflict edges are read from ``conflicts.tsv`` produced by
``peptide_fasta_evaluator.py``.

The strategy is intentionally conservative: if either endpoint participates in
a conflict meeting the selected Jaccard and LCS thresholds, both endpoints are
removed. This is not a minimum vertex-cover solver.

Validated defaults are Jaccard >= 0.5 and LCS fraction >= 0.5. Thresholds below
0.5 are intentionally unsupported. Reciprocal-conflict filtering is not part of
the current workflow. Input FASTA headers are preserved for surviving peptides.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Iterator, Tuple


def similarity_threshold(value: str) -> float:
    try:
        threshold = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Similarity threshold must be numeric.") from exc
    if not 0.5 <= threshold <= 1.0:
        raise argparse.ArgumentTypeError(
            "Similarity threshold must be between 0.5 and 1.0 inclusive."
        )
    return threshold


# ----------------------------
# FASTA streaming I/O
# ----------------------------

def iter_fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    """
    Stream FASTA records as (header, sequence) preserving header text exactly (minus leading '>').

    - Concatenates multi-line sequences
    - Uppercases sequences
    - Skips empty records
    """
    header = None
    seq_parts = []

    def flush():
        nonlocal header, seq_parts
        if header is None:
            seq_parts = []
            return None
        seq = "".join(seq_parts).strip().upper()
        seq_parts = []
        if not seq:
            return None
        return header, seq

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                rec = flush()
                if rec is not None:
                    yield rec
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    rec = flush()
    if rec is not None:
        yield rec


def write_fasta_records(records: Iterable[Tuple[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


# ----------------------------
# Conflicts parsing
# ----------------------------

def load_removal_set(
    conflicts_tsv: Path,
    min_jaccard: float,
    min_lcs: float,
) -> Tuple[Dict[str, int], int]:
    """
    Returns:
      - counts: peptide -> number of qualifying edges it participates in
      - qualifying_edges: number of qualifying edges observed (row-level, not deduped undirected)
    """
    counts: Counter[str] = Counter()
    qualifying_edges = 0

    with conflicts_tsv.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"t1", "t2", "jaccard_kmers", "lcs_frac"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"conflicts.tsv missing required columns: {sorted(missing)}")

        for row in r:
            t1 = row["t1"].strip()
            t2 = row["t2"].strip()
            if not t1 or not t2:
                continue

            jac = float(row["jaccard_kmers"])
            lcs = float(row["lcs_frac"])
            if jac < min_jaccard:
                continue
            if min_lcs > 0 and lcs < min_lcs:
                continue

            qualifying_edges += 1
            counts[t1] += 1
            counts[t2] += 1

    return dict(counts), qualifying_edges


def write_kv_tsv(rows: Iterable[Tuple[str, object]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in rows:
            w.writerow([k, v])


def write_remaining_edges(
    conflicts_tsv: Path,
    removed: set[str],
    min_jaccard: float,
    min_lcs: float,
    out_path: Path,
) -> int:
    """
    Re-scan conflicts.tsv and write any qualifying edges whose endpoints both survived.
    Returns number of remaining qualifying edges.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    remaining = 0

    with conflicts_tsv.open("r", encoding="utf-8", newline="") as f_in, out_path.open("w", encoding="utf-8", newline="") as f_out:
        r = csv.DictReader(f_in, delimiter="\t")
        w = csv.writer(f_out, delimiter="\t")
        w.writerow(["t1", "t2", "jaccard_kmers", "lcs_frac"])

        for row in r:
            t1 = row.get("t1", "").strip()
            t2 = row.get("t2", "").strip()
            if not t1 or not t2:
                continue

            jac = float(row["jaccard_kmers"])
            lcs = float(row["lcs_frac"])
            if jac < min_jaccard:
                continue
            if min_lcs > 0 and lcs < min_lcs:
                continue

            if (t1 in removed) or (t2 in removed):
                continue

            remaining += 1
            w.writerow([t1, t2, jac, lcs])

    return remaining


# ----------------------------
# Main
# ----------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Fast conservative culling: remove any peptide involved in a qualifying conflict edge."
    )
    ap.add_argument("--peptide_fasta", required=True, type=Path)
    ap.add_argument("--conflicts_tsv", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)

    ap.add_argument("--min_jaccard", type=similarity_threshold, default=0.5, help="Minimum Jaccard similarity (default: 0.5; minimum: 0.5)")
    ap.add_argument("--min_lcs", type=similarity_threshold, default=0.5, help="Minimum LCS fraction (default: 0.5; minimum: 0.5)")

    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # 1) Determine which peptides must be removed
    remove_counts, qualifying_edges = load_removal_set(
        conflicts_tsv=args.conflicts_tsv,
        min_jaccard=args.min_jaccard,
        min_lcs=args.min_lcs,
    )
    removed = set(remove_counts.keys())

    # 2) Stream FASTA -> write kept records (preserving headers)
    out_fa = args.out_dir / "filtered_target.fa"

    total_records = 0
    kept_records = 0
    removed_records = 0

    def kept_iter():
        nonlocal total_records, kept_records, removed_records
        for hdr, pep in iter_fasta_records(args.peptide_fasta):
            total_records += 1
            if pep in removed:
                removed_records += 1
                continue
            kept_records += 1
            yield (hdr, pep)

    write_fasta_records(kept_iter(), out_fa)

    # 3) Write removed log (with counts)
    removed_tsv = args.out_dir / "removed_peptides.tsv"
    with removed_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "qualifying_conflict_count"])
        for p in sorted(removed):
            w.writerow([p, remove_counts.get(p, 0)])

    # 4) Remaining edges (should be zero unless conflicts has peptides not in FASTA, or FASTA contains peptides not referenced)
    remaining_edges_tsv = args.out_dir / "remaining_conflicts_after_cull.tsv"
    edges_remaining = write_remaining_edges(
        conflicts_tsv=args.conflicts_tsv,
        removed=removed,
        min_jaccard=args.min_jaccard,
        min_lcs=args.min_lcs,
        out_path=remaining_edges_tsv,
    )

    # 5) Report
    report_rows = [
        ("total_records_in_fasta", total_records),
        ("qualifying_edges_seen", qualifying_edges),
        ("removed_peptides_unique", len(removed)),
        ("removed_records_in_fasta", removed_records),
        ("kept_records_in_fasta", kept_records),
        ("edges_remaining_after_cull", edges_remaining),
        ("min_jaccard", args.min_jaccard),
        ("min_lcs", args.min_lcs),
        ("strategy", "remove_both_endpoints_of_any_qualifying_edge"),
    ]
    write_kv_tsv(report_rows, args.out_dir / "cull_report.tsv")

    print(f"Wrote: {out_fa}")
    print(f"Wrote: {removed_tsv}")
    print(f"Wrote: {args.out_dir/'cull_report.tsv'}")
    print(f"Wrote: {remaining_edges_tsv}")
    print(f"Removed {len(removed)} unique peptides; kept {kept_records}/{total_records} FASTA records.")
    if edges_remaining != 0:
        print("NOTE: Some qualifying edges remain. This usually means conflicts.tsv references peptides not present in the FASTA,")
        print("      or you are filtering conflicts with thresholds different from the FASTA's universe.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
