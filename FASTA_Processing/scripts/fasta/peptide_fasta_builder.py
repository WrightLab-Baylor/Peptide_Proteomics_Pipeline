#!/usr/bin/env python3
"""
peptide_fasta_builder.py

Digest a protein FASTA into a deduplicated canonical peptide FASTA and write the
peptide-to-protein digest map used by downstream Binning and Targeted Peptide
Removal workflows.

Default digestion:
- Trypsin-like cleavage after K/R unless the next residue is P.
- Peptide length 6-50 aa.
- Peptides containing non-canonical residues are excluded.
- Target/decoy-overlap removal is optional and OFF by default.
- Peptide FASTA headers equal the peptide sequence.

Outputs:
- <input_stem>_peptide.fasta
- peptide_protein_digest_map.tsv
- digest_stats.txt
- removed_peptides.tsv (unless disabled)
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Generator, Iterable, List, Optional, Set, Tuple


@dataclass(frozen=True)
class FastaRecord:
    header: str
    seq: str


CANONICAL_AAS: Set[str] = set("ACDEFGHIKLMNPQRSTVWY")


# -----------------------------
# FASTA I/O
# -----------------------------
def read_fasta(path: Path) -> Generator[FastaRecord, None, None]:
    header: Optional[str] = None
    seq_chunks: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header=header, seq="".join(seq_chunks))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if header is not None:
        yield FastaRecord(header=header, seq="".join(seq_chunks))


def write_peptide_fasta(peptides: Iterable[str], out_fasta: Path, line_width: int = 80) -> None:
    with out_fasta.open("w", encoding="utf-8") as f:
        for pep in peptides:
            f.write(f">{pep}\n")
            for i in range(0, len(pep), line_width):
                f.write(pep[i:i + line_width] + "\n")


def write_mapping_tsv(pep_to_proteins: Dict[str, Set[str]], out_tsv: Path) -> None:
    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("Peptide\tProteinCount\tProteins\n")
        for pep in sorted(pep_to_proteins.keys()):
            prots = sorted(pep_to_proteins[pep])
            f.write(f"{pep}\t{len(prots)}\t{';'.join(prots)}\n")


def write_removed_peptides_tsv(
    removed_pep_to_proteins: Dict[str, Set[str]],
    out_tsv: Path,
    reasons: Dict[str, str],
    partners: Dict[str, str],
) -> None:
    """
    Removed peptide log with reasons.
    Columns:
      Peptide  ProteinCount  Proteins  Reason  Partner
    """
    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("Peptide\tProteinCount\tProteins\tReason\tPartner\n")
        for pep in sorted(removed_pep_to_proteins.keys()):
            prots = sorted(removed_pep_to_proteins[pep])
            reason = reasons.get(pep, "")
            partner = partners.get(pep, "")
            f.write(f"{pep}\t{len(prots)}\t{';'.join(prots)}\t{reason}\t{partner}\n")


# -----------------------------
# Sequence helpers
# -----------------------------
def is_canonical_peptide(seq: str) -> bool:
    return all(c in CANONICAL_AAS for c in seq)


def protein_id_from_header(header: str) -> str:
    return header.split()[0]


def reverse_peptide(seq: str) -> str:
    return seq[::-1]


# -----------------------------
# Trypsin-like digest used by the production workflow
# -----------------------------
def trypsin_cut_sites_with_noncanonical_block(seq: str) -> List[int]:
    """
    Cut sites include 0 and len(seq).
    Trypsin rule: cleave after K/R unless next residue is P.
    Extended rule: treat any non-canonical residue like Proline in the NEXT position
    (i.e., blocks cleavage).
    """
    cut_sites = [0]
    L = len(seq)

    for i, aa in enumerate(seq):
        if aa in ("K", "R"):
            next_aa = seq[i + 1] if i + 1 < L else None
            if next_aa is None:
                cut_sites.append(i + 1)
                continue

            if next_aa == "P" or (next_aa not in CANONICAL_AAS):
                continue

            cut_sites.append(i + 1)

    if cut_sites[-1] != L:
        cut_sites.append(L)

    return cut_sites


def trypsin_digest_v11(
    protein_seq: str,
    min_len: int,
    max_len: int,
) -> List[str]:
    seq = protein_seq.strip().upper()
    if not seq:
        return []

    cut_sites = trypsin_cut_sites_with_noncanonical_block(seq)

    peps: List[str] = []
    for start, end in zip(cut_sites[:-1], cut_sites[1:]):
        pep = seq[start:end]
        if min_len <= len(pep) <= max_len:
            peps.append(pep)
    return peps


# -----------------------------
# Build peptide universe + removed-peptide bookkeeping
# -----------------------------
def build_peptide_universe(
    fasta_path: Path,
    min_len: int,
    max_len: int,
) -> Tuple[
    Set[str],                    # unique canonical peptides (output universe)
    Dict[str, Set[str]],         # canonical peptide -> proteins
    Dict[str, Set[str]],         # removed peptides -> proteins
    Dict[str, str],              # removed peptide -> reason
    Dict[str, str],              # removed peptide -> partner (for overlaps)
    Dict[str, int],              # stats
]:
    unique_peptides: Set[str] = set()
    pep_to_proteins: Dict[str, Set[str]] = {}

    removed_pep_to_proteins: Dict[str, Set[str]] = {}
    removed_reasons: Dict[str, str] = {}
    removed_partners: Dict[str, str] = {}

    n_proteins = 0
    n_proteins_with_noncanonical = 0
    n_total_peptides_pre_dedup = 0

    n_removed_noncanonical_total = 0
    removed_noncanonical_unique: Set[str] = set()

    for rec in read_fasta(fasta_path):
        n_proteins += 1
        prot_id = protein_id_from_header(rec.header)
        seq_u = rec.seq.strip().upper()
        if not seq_u:
            continue

        if any(c not in CANONICAL_AAS for c in seq_u):
            n_proteins_with_noncanonical += 1

        peptides = trypsin_digest_v11(seq_u, min_len=min_len, max_len=max_len)
        n_total_peptides_pre_dedup += len(peptides)

        for pep in peptides:
            if is_canonical_peptide(pep):
                unique_peptides.add(pep)
                pep_to_proteins.setdefault(pep, set()).add(prot_id)
            else:
                n_removed_noncanonical_total += 1
                removed_noncanonical_unique.add(pep)
                removed_pep_to_proteins.setdefault(pep, set()).add(prot_id)
                removed_reasons.setdefault(pep, "noncanonical")
                removed_partners.setdefault(pep, "")

    stats = {
        "proteins_read": n_proteins,
        "proteins_with_noncanonical_residues": n_proteins_with_noncanonical,
        "peptides_generated_pre_dedup": n_total_peptides_pre_dedup,
        "unique_canonical_peptides_output_pre_tda": len(unique_peptides),
        "noncanonical_peptides_removed_total": n_removed_noncanonical_total,
        "noncanonical_peptides_removed_unique": len(removed_noncanonical_unique),
    }
    return unique_peptides, pep_to_proteins, removed_pep_to_proteins, removed_reasons, removed_partners, stats


def apply_tda_overlap_removal(
    unique_peptides: Set[str],
    pep_to_proteins: Dict[str, Set[str]],
    removed_pep_to_proteins: Dict[str, Set[str]],
    removed_reasons: Dict[str, str],
    removed_partners: Dict[str, str],
) -> Dict[str, int]:
    """
    Remove peptides that would overlap with reverse-decoy peptides.
    Removal policy:
      - If p and r=reverse(p) both exist, remove BOTH p and r.
      - Palindromes (p == r) are removed.
    """
    # Identify pairs once
    pairs: List[Tuple[str, str]] = []
    palindromes = 0

    for p in unique_peptides:
        r = reverse_peptide(p)
        if r in unique_peptides and p <= r:
            pairs.append((p, r))
            if p == r:
                palindromes += 1

    # Build removal set (both members)
    to_remove: Set[str] = set()
    for a, b in pairs:
        to_remove.add(a)
        to_remove.add(b)

    # Log + remove
    for p in to_remove:
        partner = reverse_peptide(p)
        removed_partners[p] = partner
        removed_reasons[p] = "tda_overlap"
        if p in pep_to_proteins:
            removed_pep_to_proteins.setdefault(p, set()).update(pep_to_proteins[p])

    for p in to_remove:
        unique_peptides.discard(p)
        pep_to_proteins.pop(p, None)

    return {
        "tda_overlap_pairs": len(pairs),
        "tda_overlap_removed_unique": len(to_remove),
        "tda_overlap_palindromes": palindromes,
    }


# -----------------------------
# CLI
# -----------------------------
def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Digest a protein FASTA with trypsin rules, deduplicate canonical peptides, write peptide FASTA (header==sequence), mapping TSV, and removed peptides log. Optional removal of target/decoy overlaps for reverse decoys."
    )
    p.add_argument("-i", "--input_fasta", required=True, type=Path, help="Input protein FASTA (.fa/.fasta/.faa)")
    p.add_argument("-o", "--out_dir", required=True, type=Path, help="Output directory")
    p.add_argument("--min_len", type=int, default=6, help="Minimum peptide length (default: 6)")
    p.add_argument("--max_len", type=int, default=50, help="Maximum peptide length (default: 50)")

    p.add_argument(
        "--remove_tda_overlaps",
        action="store_true",
        help="Remove peptides that would overlap with reverse-decoy peptides (target/decoy overlap), including palindromes.",
    )

    # New general opt-out flag (keep old one as alias so existing calls don't break)
    p.add_argument(
        "--no_log_removed_peptides",
        action="store_true",
        help="Disable writing removed_peptides.tsv (default: enabled).",
    )
    p.add_argument(
        "--no_log_removed_noncanonical",
        action="store_true",
        help="Alias for --no_log_removed_peptides (kept for backward compatibility).",
    )

    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if not args.input_fasta.exists():
        print(f"ERROR: input FASTA not found: {args.input_fasta}", file=sys.stderr)
        return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)

    unique_peptides, pep_to_proteins, removed_pep_to_proteins, removed_reasons, removed_partners, stats = build_peptide_universe(
        fasta_path=args.input_fasta,
        min_len=args.min_len,
        max_len=args.max_len,
    )

    # Optional TDA overlap removal
    if args.remove_tda_overlaps:
        tda_stats = apply_tda_overlap_removal(
            unique_peptides=unique_peptides,
            pep_to_proteins=pep_to_proteins,
            removed_pep_to_proteins=removed_pep_to_proteins,
            removed_reasons=removed_reasons,
            removed_partners=removed_partners,
        )
        stats.update(tda_stats)
    else:
        stats.update({
            "tda_overlap_pairs": 0,
            "tda_overlap_removed_unique": 0,
            "tda_overlap_palindromes": 0,
        })
    # ---- Derived TDA overlap rate (fraction + percent, capped precision) ----
    pre = stats.get("unique_canonical_peptides_output_pre_tda", 0)
    removed = stats.get("tda_overlap_removed_unique", 0)

    if pre > 0:
        frac = removed / pre
        stats["tda_overlap_removed_fraction"] = round(frac, 10)
        stats["tda_overlap_removed_percent"] = round(100.0 * frac, 10)
    else:
        stats["tda_overlap_removed_fraction"] = 0.0
        stats["tda_overlap_removed_percent"] = 0.0

    # ---- Record toggle state ----
    stats["tda_overlap_removal_enabled"] = int(args.remove_tda_overlaps)

    sorted_peptides = sorted(unique_peptides)

    input_stem = args.input_fasta.stem
    out_fasta = args.out_dir / f"{input_stem}_peptide.fasta"
    out_map = args.out_dir / "peptide_protein_digest_map.tsv"
    out_stats = args.out_dir / "digest_stats.txt"
    out_removed = args.out_dir / "removed_peptides.tsv"

    write_peptide_fasta(sorted_peptides, out_fasta)
    write_mapping_tsv(pep_to_proteins, out_map)

    log_enabled = not (args.no_log_removed_peptides or args.no_log_removed_noncanonical)
    if log_enabled:
        write_removed_peptides_tsv(
            removed_pep_to_proteins=removed_pep_to_proteins,
            out_tsv=out_removed,
            reasons=removed_reasons,
            partners=removed_partners,
        )

    # Final stats fields
    stats["unique_canonical_peptides_output_post_tda"] = len(unique_peptides)
    stats["removed_peptides_log_written"] = int(log_enabled)
    stats["removed_noncanonical_log_written"] = int(log_enabled)  # kept for compatibility with older parsing

    with out_stats.open("w", encoding="utf-8") as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")

    print("Done.")
    print(f"Peptide FASTA:      {out_fasta}")
    print(f"Digest Map (TSV):   {out_map}")
    print(f"Digest Stats:       {out_stats}")
    if log_enabled:
        print(f"Removed peptides:   {out_removed}")
    if args.remove_tda_overlaps:
        print(f"TDA overlaps removed (unique): {stats['tda_overlap_removed_unique']}")
        print(f"Palindromes removed:          {stats['tda_overlap_palindromes']}")
    print(f"Unique peptides (final): {stats['unique_canonical_peptides_output_post_tda']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
