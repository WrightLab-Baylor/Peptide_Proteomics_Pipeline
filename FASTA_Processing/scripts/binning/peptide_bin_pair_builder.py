#!/usr/bin/env python3
"""
peptide_bin_pair_builder.py

Randomly shuffles a peptide FASTA once, partitions that master shuffle into
roughly equal-sized groups, and writes group FASTAs, pairwise group-combination
FASTAs, or both in a pipeline-compatible folder structure.

Example output layout for -o binning_experiment and 9 groups:

binning_experiment/
  summary.txt
  group_assignments.tsv
  group_files/                      # default behavior for original groups
    group_1.fasta
    group_2.fasta
    ...
  binning_experiment_pair_1_2/
    database/
      pair_1_2.fasta
    data/
  binning_experiment_pair_1_3/
    database/
      pair_1_3.fasta
    data/
  ...

Optional alternate layout for original groups (--group_subfolders):

binning_experiment/
  binning_experiment_group_1/
    database/
      group_1.fasta
    data/
  ...

Notes:
- Designed for zero-missed-cleavage peptide FASTA inputs where each peptide can
  be shuffled independently.
- Exactly one master shuffle is generated per run. Groups are contiguous slices
  of that shuffle, and pair bins are concatenations of those same groups.
- FASTA headers are ignored on input; output FASTA headers are set equal to the
  peptide sequence for consistency with the existing peptide FASTA workflow.
- Group sizes differ by at most 1.
- Pairwise bins are all unique unordered pairs: (1,2), (1,3), ..., (n-1,n).
- --output-mode both preserves the original behavior and is the default.
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from itertools import combinations
from pathlib import Path
from typing import Iterable, List, Optional, Sequence


# -----------------------------
# FASTA I/O
# -----------------------------
def read_peptide_fasta(path: Path) -> List[str]:
    """Read peptide sequences from a FASTA file, ignoring headers."""
    if not path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {path}")

    peptides: List[str] = []
    seq_chunks: List[str] = []
    in_record = False

    def flush() -> None:
        nonlocal seq_chunks
        if not seq_chunks:
            return
        seq = "".join(seq_chunks).strip().upper()
        seq_chunks = []
        if seq:
            peptides.append(seq)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_record:
                    flush()
                in_record = True
            else:
                in_record = True
                seq_chunks.append(line)

    if in_record:
        flush()

    return peptides


def write_peptide_fasta(peptides: Iterable[str], out_fasta: Path, line_width: int = 80) -> None:
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with out_fasta.open("w", encoding="utf-8") as handle:
        for pep in peptides:
            handle.write(f">{pep}\n")
            for i in range(0, len(pep), line_width):
                handle.write(pep[i:i + line_width] + "\n")


# -----------------------------
# Grouping helpers
# -----------------------------
def partition_evenly(items: Sequence[str], n_groups: int) -> List[List[str]]:
    """
    Partition items into n_groups such that sizes differ by at most 1.
    Earlier groups receive the +1 remainder items.
    """
    total = len(items)
    base = total // n_groups
    remainder = total % n_groups

    groups: List[List[str]] = []
    start = 0
    for idx in range(n_groups):
        size = base + (1 if idx < remainder else 0)
        end = start + size
        groups.append(list(items[start:end]))
        start = end
    return groups


def summarize_group_sizes(groups: Sequence[Sequence[str]]) -> List[int]:
    return [len(g) for g in groups]


# -----------------------------
# Output writers
# -----------------------------
def write_group_assignments(groups: Sequence[Sequence[str]], out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as handle:
        handle.write("Peptide\tGroup\n")
        for idx, group in enumerate(groups, start=1):
            for pep in group:
                handle.write(f"{pep}\t{idx}\n")


def make_empty_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_group_files_flat(groups: Sequence[Sequence[str]], root_dir: Path) -> None:
    group_root = root_dir / "group_files"
    group_root.mkdir(parents=True, exist_ok=True)

    for idx, group in enumerate(groups, start=1):
        out_fasta = group_root / f"group_{idx}.fasta"
        write_peptide_fasta(group, out_fasta)


def write_group_files_subfolders(groups: Sequence[Sequence[str]], root_dir: Path, experiment_name: str) -> None:
    for idx, group in enumerate(groups, start=1):
        run_dir = root_dir / f"{experiment_name}_group_{idx}"
        db_dir = run_dir / "database"
        data_dir = run_dir / "data"
        db_dir.mkdir(parents=True, exist_ok=True)
        make_empty_dir(data_dir)
        out_fasta = db_dir / f"group_{idx}.fasta"
        write_peptide_fasta(group, out_fasta)


def write_pair_bins(groups: Sequence[Sequence[str]], root_dir: Path, experiment_name: str) -> List[tuple[int, int, int]]:
    pair_stats: List[tuple[int, int, int]] = []

    for i, j in combinations(range(1, len(groups) + 1), 2):
        pair_dir = root_dir / f"{experiment_name}_pair_{i}_{j}"
        db_dir = pair_dir / "database"
        data_dir = pair_dir / "data"
        db_dir.mkdir(parents=True, exist_ok=True)
        make_empty_dir(data_dir)

        pair_peptides = list(groups[i - 1]) + list(groups[j - 1])
        out_fasta = db_dir / f"pair_{i}_{j}.fasta"
        write_peptide_fasta(pair_peptides, out_fasta)
        pair_stats.append((i, j, len(pair_peptides)))

    return pair_stats


def write_summary(
    out_txt: Path,
    input_fasta: Path,
    total_peptides: int,
    n_groups: int,
    seed: Optional[int],
    group_sizes: Sequence[int],
    pair_stats: Sequence[tuple[int, int, int]],
    group_subfolders: bool,
    output_mode: str,
) -> None:
    pair_sizes = [size for _, _, size in pair_stats]
    pairs_requested = output_mode in {"pairs", "both"}
    groups_requested = output_mode in {"groups", "both"}
    expected_pairs = math.comb(n_groups, 2) if pairs_requested else 0

    with out_txt.open("w", encoding="utf-8") as handle:
        handle.write("Peptide binning summary\n")
        handle.write("=======================\n\n")
        handle.write(f"Input FASTA:\t{input_fasta}\n")
        handle.write(f"Total peptides read:\t{total_peptides}\n")
        handle.write(f"Requested groups:\t{n_groups}\n")
        handle.write(f"Random seed:\t{seed if seed is not None else 'system_random'}\n")
        handle.write(f"Original group output mode:\t{'subfolders' if group_subfolders else 'group_files'}\n")
        handle.write("\n")

        handle.write("Group sizes\n")
        handle.write("-----------\n")
        for idx, size in enumerate(group_sizes, start=1):
            handle.write(f"group_{idx}\t{size}\n")
        handle.write(f"min_group_size\t{min(group_sizes) if group_sizes else 0}\n")
        handle.write(f"max_group_size\t{max(group_sizes) if group_sizes else 0}\n")
        handle.write("\n")

        handle.write("Pair sizes\n")
        handle.write("----------\n")
        handle.write(f"total_pairs\t{len(pair_stats)}\n")
        handle.write(f"expected_pairs_n_choose_2\t{expected_pairs}\n")
        handle.write(f"min_pair_size\t{min(pair_sizes) if pair_sizes else 0}\n")
        handle.write(f"max_pair_size\t{max(pair_sizes) if pair_sizes else 0}\n")
        handle.write("\n")
        for i, j, size in pair_stats:
            handle.write(f"pair_{i}_{j}\t{size}\n")


# -----------------------------
# CLI
# -----------------------------
def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Shuffle peptide FASTA entries once, divide them into roughly equal groups, "
            "and write group FASTAs, pairwise group-combination FASTAs, or both."
        )
    )
    parser.add_argument(
        "-i",
        "--input_fasta",
        required=True,
        type=Path,
        help="Input peptide FASTA file",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        required=True,
        type=Path,
        help=(
            "Root output directory for the experiment. The final path component is used "
            "as the experiment name in subfolder naming. Example: binning_experiment"
        ),
    )
    parser.add_argument(
        "-g",
        "--num_groups",
        required=True,
        type=int,
        help="Number of random groups to create",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducible shuffling",
    )
    parser.add_argument(
        "--group_subfolders",
        action="store_true",
        help=(
            "When groups are requested, write them using the same run-subfolder layout "
            "as pair bins. Default behavior writes them to <out_dir>/group_files/."
        ),
    )
    parser.add_argument(
        "--output-mode",
        choices=("groups", "pairs", "both"),
        default="both",
        help=(
            "Outputs to write after creating one master shuffle and its groups: "
            "'groups' writes only group FASTAs; 'pairs' writes only pair-bin FASTAs; "
            "'both' writes both and preserves the original behavior (default: both)."
        ),
    )
    return parser.parse_args(argv)


# -----------------------------
# Main
# -----------------------------
def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    minimum_groups = 2 if args.output_mode in {"pairs", "both"} else 1
    if args.num_groups < minimum_groups:
        if args.output_mode in {"pairs", "both"}:
            message = "--num_groups must be at least 2 when pair bins are requested."
        else:
            message = "--num_groups must be at least 1."
        print(f"ERROR: {message}", file=sys.stderr)
        return 2

    peptides = read_peptide_fasta(args.input_fasta)
    if not peptides:
        print("ERROR: No peptide sequences were read from the input FASTA.", file=sys.stderr)
        return 2

    if args.num_groups > len(peptides):
        print(
            f"ERROR: Requested {args.num_groups} groups but only {len(peptides)} peptides were read. "
            "Cannot create non-empty groups.",
            file=sys.stderr,
        )
        return 2

    rng = random.Random(args.seed)
    shuffled = list(peptides)
    rng.shuffle(shuffled)

    groups = partition_evenly(shuffled, args.num_groups)
    root_dir = args.out_dir
    experiment_name = root_dir.name
    root_dir.mkdir(parents=True, exist_ok=True)

    # Group FASTAs are written from the already-created master-shuffle groups.
    groups_requested = args.output_mode in {"groups", "both"}
    pairs_requested = args.output_mode in {"pairs", "both"}

    if groups_requested:
        if args.group_subfolders:
            write_group_files_subfolders(groups, root_dir, experiment_name)
        else:
            write_group_files_flat(groups, root_dir)

    # Always retain the assignment table so pair-only runs remain traceable to
    # the same master shuffle and group boundaries used to construct pair bins.
    write_group_assignments(groups, root_dir / "group_assignments.tsv")

    # Pair bins are direct concatenations of the same in-memory groups.
    pair_stats = (
        write_pair_bins(groups, root_dir, experiment_name)
        if pairs_requested
        else []
    )

    # Summary file
    group_sizes = summarize_group_sizes(groups)
    write_summary(
        out_txt=root_dir / "summary.txt",
        input_fasta=args.input_fasta,
        total_peptides=len(peptides),
        n_groups=args.num_groups,
        seed=args.seed,
        group_sizes=group_sizes,
        pair_stats=pair_stats,
        group_subfolders=args.group_subfolders,
        output_mode=args.output_mode,
    )

    print("Done.")
    print(f"Output root:        {root_dir}")
    print(f"Experiment name:    {experiment_name}")
    print(f"Input peptides:     {len(peptides)}")
    print(f"Output mode:        {args.output_mode}")
    print("Shuffle strategy:   one master shuffle, then contiguous groups")
    print(f"Groups formed:      {len(groups)}")
    print(f"Group FASTAs written: {len(groups) if groups_requested else 0}")
    print(f"Pair bins written:  {len(pair_stats)}")
    print(f"Summary file:       {root_dir / 'summary.txt'}")
    print(f"Assignments TSV:    {root_dir / 'group_assignments.tsv'}")
    if groups_requested:
        if args.group_subfolders:
            print("Original groups:    written as per-group run subfolders")
        else:
            print(f"Original groups:    {root_dir / 'group_files'}")
    else:
        print("Original groups:    FASTA files not requested")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
