#!/usr/bin/env python3
"""
peptide_fasta_evaluator.py

Memory-scaled peptide conflict evaluator for targeted peptide-removal workflows.

The evaluator compares each target peptide's reverse-decoy sequence against the
other target peptides using sparse k-mer candidate generation followed by
Jaccard and LCS similarity scoring. It can also report target/decoy overlap and
palindrome QC metrics. Those QC observations are reported, not automatically
removed by this script.

Validated default conflict parameters:
- k-mer size: 3
- minimum shared k-mers: 3
- Jaccard threshold: 0.5
- LCS-fraction threshold: 0.5
- maximum candidates per query peptide: 500

Similarity thresholds below 0.5 are intentionally rejected because candidate
explosion makes this method impractical at lower similarity. Reciprocal-conflict
filtering is not part of the current workflow.

For large production datasets, use the ``conflicts`` subcommand and then apply
threshold-specific removal with ``conflict_culler_fast.py``.
"""

from __future__ import annotations

import argparse
import csv
import math
import multiprocessing as mp
import statistics
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

try:
    from tqdm import tqdm  # type: ignore
except Exception:
    tqdm = None


# ----------------------------
# Mass constants (unchanged)
# ----------------------------

AA_MASS_MONO = {
    "A": 71.037113805,
    "R": 156.101111050,
    "N": 114.042927470,
    "D": 115.026943065,
    "C": 103.009184505,
    "E": 129.042593135,
    "Q": 128.058577540,
    "G": 57.021463735,
    "H": 137.058911875,
    "I": 113.084064015,
    "L": 113.084064015,
    "K": 128.094963050,
    "M": 131.040484645,
    "F": 147.068413945,
    "P": 97.052763875,
    "S": 87.032028435,
    "T": 101.047678505,
    "W": 186.079312980,
    "Y": 163.063328575,
    "V": 99.068413945,
}
H2O_MASS = 18.010564684


# ----------------------------
# Data structures
# ----------------------------

@dataclass(frozen=True)
class Conflict:
    t1: str
    d1: str
    t2: str
    d2: str
    jaccard: float
    lcs_frac: float
    shared_kmers: int


# ----------------------------
# CLI
# ----------------------------
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    if v.lower() in ("no", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected (True/False)")


def similarity_threshold(value: str) -> float:
    """Parse a supported similarity threshold in the validated range [0.5, 1.0]."""
    try:
        threshold = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Similarity threshold must be numeric.") from exc
    if not 0.5 <= threshold <= 1.0:
        raise argparse.ArgumentTypeError(
            "Similarity threshold must be between 0.5 and 1.0 inclusive."
        )
    return threshold


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Memory-scaled peptide target/decoy similarity QC, conflict detection, and optional pruning."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--target_pep_fasta", required=True, type=Path, help="Target peptide FASTA (peptides as sequences)")
    common.add_argument("--out_dir", required=True, type=Path, help="Output directory")
    common.add_argument("--min_len", type=int, default=1, help="Min peptide length")
    common.add_argument("--max_len", type=int, default=10_000, help="Max peptide length")
    common.add_argument("--ignore_il_equivalence", action="store_true", help="Map I->L in all comparisons")

    # QC
    p_qc = sub.add_parser("qc", parents=[common], help="Compute target–decoy similarity QC metrics")
    p_qc.add_argument("--mass_tol_ppm", type=float, default=0.0, help="Enable mass-collision metric with ppm tolerance")
    p_qc.add_argument("--mass_tol_da", type=float, default=0.0, help="Enable mass-collision metric with Da tolerance")

    # Conflicts
    p_conf = sub.add_parser("conflicts", parents=[common], help="Find cross-similarity conflicts reverse(t1) ~ t2")
    p_conf.add_argument("--k", type=int, default=3, help="k-mer size for candidate generation (default 3)")
    p_conf.add_argument("--min_shared_kmers", type=int, default=3, help="Require >= this many shared k-mers to score")
    p_conf.add_argument("--jaccard_threshold", type=similarity_threshold, default=0.5, help="Report conflicts with Jaccard >= threshold (default: 0.5; minimum: 0.5)")
    p_conf.add_argument("--lcs_threshold", type=similarity_threshold, default=0.5, help="Require LCS fraction >= threshold (default: 0.5; minimum: 0.5)")
    p_conf.add_argument("--max_candidates", type=int, default=500, help="Cap candidates per query peptide (prevents blowups)")
    p_conf.add_argument("--max_conflicts_per_peptide", type=int, default=25, help="Cap reported conflicts per t1")
    p_conf.add_argument("--self_match", action="store_true", help="Allow matching t2==t1 (usually false)")
    p_conf.add_argument("--workers", type=int, default=1, help="Worker processes (Linux fork recommended). Default: 1")

    # mk4 knobs
    p_conf.add_argument(
        "--precompute-decoy-kmers",
        dest="precompute_decoy_kmers",
        type=str2bool,
        default=True,
        help="Precompute reverse(t) k-mer lists (True/False). Default: True.",
    )


    # Prune (note: for huge runs, prefer using conflict_culler on conflicts.tsv)
    p_prune = sub.add_parser("prune", parents=[common], help="Prune target peptides to reduce conflicts; outputs filtered FASTA")
    p_prune.add_argument("--k", type=int, default=3)
    p_prune.add_argument("--min_shared_kmers", type=int, default=3)
    p_prune.add_argument("--jaccard_threshold", type=similarity_threshold, default=0.5)
    p_prune.add_argument("--lcs_threshold", type=similarity_threshold, default=0.5)
    p_prune.add_argument("--max_candidates", type=int, default=500)
    p_prune.add_argument("--max_remove_fraction", type=float, default=0.005, help="Max fraction of targets to remove (e.g., 0.005=0.5%%)")
    p_prune.add_argument("--prefer_remove_short", action="store_true", help="Tie-break: remove shorter peptides first (recommended)")
    p_prune.add_argument("--write_decoy_fasta", action="store_true", help="Also write explicit decoy FASTA")
    p_prune.add_argument("--workers", type=int, default=1, help="Worker processes (Linux fork recommended). Default: 1")
    p_prune.add_argument(
        "--precompute-decoy-kmers",
        dest="precompute_decoy_kmers",
        type=str2bool,
        default=True,
        help="Precompute reverse(t) k-mer lists (True/False). Default: True.",
    )

    return p.parse_args(argv)


# ----------------------------
# FASTA I/O + normalization
# ----------------------------

def normalize_peptide(pep: str, ignore_il: bool) -> str:
    pep = pep.strip().upper()
    if ignore_il:
        pep = pep.replace("I", "L")
    return pep


def is_valid_peptide(pep: str) -> bool:
    return pep.isalpha() and pep.isupper()


def read_peptide_fasta(path: Path, min_len: int, max_len: int, ignore_il: bool) -> List[str]:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    peps: List[str] = []
    seq_parts: List[str] = []
    in_record = False

    def flush():
        nonlocal seq_parts
        if not seq_parts:
            return
        seq = "".join(seq_parts)
        seq_parts = []
        pep = normalize_peptide(seq, ignore_il)
        if not pep:
            return
        if not is_valid_peptide(pep):
            return
        if not (min_len <= len(pep) <= max_len):
            return
        peps.append(pep)

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_record:
                    flush()
                in_record = True
            else:
                in_record = True
                seq_parts.append(line)

    if in_record:
        flush()

    return peps


def write_peptide_fasta(peptides: Iterable[str], out_path: Path, prefix: str) -> None:
    with out_path.open("w", encoding="utf-8") as f:
        for i, pep in enumerate(peptides, start=1):
            f.write(f">{prefix}{i}\n{pep}\n")


# ----------------------------
# QC metrics (unchanged from mk3)
# ----------------------------

def peptide_mono_mass(pep: str) -> float:
    mass = H2O_MASS
    for aa in pep:
        m = AA_MASS_MONO.get(aa)
        if m is None:
            return float("nan")
        mass += m
    return mass


def ks_statistic(a: List[float], b: List[float]) -> float:
    a = [x for x in a if math.isfinite(x)]
    b = [x for x in b if math.isfinite(x)]
    a.sort()
    b.sort()
    n = len(a)
    m = len(b)
    if n == 0 or m == 0:
        return float("nan")
    values = sorted(set(a) | set(b))
    i = j = 0
    d = 0.0
    for v in values:
        while i < n and a[i] <= v:
            i += 1
        while j < m and b[j] <= v:
            j += 1
        d = max(d, abs(i / n - j / m))
    return d


def normalize_counter(c: Dict[str, int]) -> Dict[str, float]:
    tot = sum(c.values())
    if tot == 0:
        return {}
    return {k: v / tot for k, v in c.items()}


def js_divergence(p: Dict[str, float], q: Dict[str, float], eps: float = 1e-12) -> float:
    keys = set(p) | set(q)
    if not keys:
        return float("nan")

    def kl(a: Dict[str, float], b: Dict[str, float]) -> float:
        s = 0.0
        for k in keys:
            ak = a.get(k, 0.0)
            if ak <= 0:
                continue
            bk = b.get(k, 0.0) + eps
            s += ak * math.log2(ak / bk)
        return s

    m = {k: 0.5 * (p.get(k, 0.0) + q.get(k, 0.0)) for k in keys}
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def aa_freq(peps: Iterable[str]) -> Dict[str, float]:
    c: Dict[str, int] = {}
    for pep in peps:
        for ch in pep:
            c[ch] = c.get(ch, 0) + 1
    return normalize_counter(c)


def dipep_freq(peps: Iterable[str]) -> Dict[str, float]:
    c: Dict[str, int] = {}
    for pep in peps:
        for i in range(len(pep) - 1):
            di = pep[i:i + 2]
            c[di] = c.get(di, 0) + 1
    return normalize_counter(c)


def compute_mass_collision_fraction(target_masses: List[float], decoy_masses: List[float], tol_da: float) -> float:
    t = sorted(target_masses)
    d = sorted(decoy_masses)
    if not t or not d:
        return float("nan")
    j = 0
    hit = 0
    for mt in t:
        while j < len(d) and d[j] < mt - tol_da:
            j += 1
        if j < len(d) and abs(d[j] - mt) <= tol_da:
            hit += 1
        elif j > 0 and abs(d[j - 1] - mt) <= tol_da:
            hit += 1
    return hit / len(t)


def _write_kv_tsv(rows: List[Tuple[str, object]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in rows:
            w.writerow([k, v])


def run_qc(target_peps: Set[str], out_dir: Path, mass_tol_ppm: float, mass_tol_da: float) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    T = set(target_peps)
    D = {p[::-1] for p in T}

    overlap = T & D
    overlap_T = len(overlap) / len(T) if T else float("nan")
    overlap_D = len(overlap) / len(D) if D else float("nan")

    len_T = [len(p) for p in T]
    len_D = [len(p) for p in D]
    ks_len = ks_statistic([float(x) for x in len_T], [float(x) for x in len_D])

    mass_T = [peptide_mono_mass(p) for p in T]
    mass_T = [m for m in mass_T if math.isfinite(m)]
    mass_D = [peptide_mono_mass(p) for p in D]
    mass_D = [m for m in mass_D if math.isfinite(m)]
    ks_mass = ks_statistic(mass_T, mass_D)

    jsd_aa = max(0.0, js_divergence(aa_freq(T), aa_freq(D)))
    jsd_dipep = max(0.0, js_divergence(dipep_freq(T), dipep_freq(D)))

    pal_T = sum(1 for p in T if p == p[::-1])
    pal_D = sum(1 for p in D if p == p[::-1])
    pal_frac_T = pal_T / len(T) if T else float("nan")
    pal_frac_D = pal_D / len(D) if D else float("nan")

    mass_collision = float("nan")
    if (mass_tol_ppm and mass_tol_ppm > 0) or (mass_tol_da and mass_tol_da > 0):
        if mass_tol_da and mass_tol_da > 0:
            tol_da = mass_tol_da
        else:
            ref = statistics.median(mass_T) if mass_T else 1000.0
            tol_da = (mass_tol_ppm * 1e-6) * ref
        mass_collision = compute_mass_collision_fraction(mass_T, mass_D, tol_da)

    report_rows = [
        ("target_unique", len(T)),
        ("decoy_unique", len(D)),
        ("overlap_count", len(overlap)),
        ("overlap_fraction_of_target", overlap_T),
        ("overlap_fraction_of_decoy", overlap_D),
        ("palindrome_fraction_target", pal_frac_T),
        ("palindrome_fraction_decoy", pal_frac_D),
        ("ks_length", ks_len),
        ("ks_mass", ks_mass),
        ("jsd_aa", jsd_aa),
        ("jsd_dipeptide", jsd_dipep),
        ("mass_collision_fraction_target", mass_collision),
        ("assumption_decoys", "peptide-level reverse(peptide)"),
        ("assumption_ptms", "not modeled (upstream library QC)"),
        ("param_mass_tol_ppm", mass_tol_ppm),
        ("param_mass_tol_da", mass_tol_da),
    ]

    out_report = out_dir / "similarity_report.tsv"
    _write_kv_tsv(report_rows, out_report)

    out_summary = out_dir / "similarity_summary.tsv"
    headers = [
        "target_unique", "decoy_unique", "overlap_count",
        "overlap_frac_T", "overlap_frac_D",
        "ks_len", "ks_mass", "jsd_aa", "jsd_dipep",
        "pal_frac_T", "pal_frac_D", "mass_collision_frac_T"
    ]
    values = [
        len(T), len(D), len(overlap),
        overlap_T, overlap_D,
        ks_len, ks_mass, jsd_aa, jsd_dipep,
        pal_frac_T, pal_frac_D, mass_collision
    ]

    write_header = not out_summary.exists()
    with out_summary.open("a", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        if write_header:
            w.writerow(headers)
        w.writerow(values)

    print(f"Wrote: {out_report}")
    print(f"Wrote: {out_summary}")


# ----------------------------
# k-mer encoding (base-20)
# ----------------------------

AA20 = "ACDEFGHIKLMNPQRSTVWY"
AA20_IDX: Dict[str, int] = {aa: i for i, aa in enumerate(AA20)}
AA20_IDX["I"] = AA20_IDX["I"]
AA20_IDX["L"] = AA20_IDX["L"]

def _kmer_index(seq: str, start: int, k: int) -> int:
    v = 0
    for i in range(start, start + k):
        idx = AA20_IDX.get(seq[i])
        if idx is None:
            return -1
        v = v * 20 + idx
    return v

def encode_kmers_unique_sorted(seq: str, k: int) -> List[int]:
    """
    Return unique, sorted k-mer IDs for seq.
    For typical peptide lengths, sorting this small list is cheap and keeps Jaccard fast.
    """
    L = len(seq)
    if k <= 0 or L < k:
        return []
    tmp: List[int] = []
    for i in range(L - k + 1):
        ki = _kmer_index(seq, i, k)
        if ki >= 0:
            tmp.append(ki)
    if not tmp:
        return []
    tmp.sort()
    # unique in-place
    out = [tmp[0]]
    for x in tmp[1:]:
        if x != out[-1]:
            out.append(x)
    return out


# ----------------------------
# Sparse k-mer store: flat buffer + offsets
# ----------------------------

def _choose_kmer_array_type(k: int) -> str:
    # For k=3, universe 8000 fits in uint16. For k=4, 160000 needs uint32.
    if k <= 3:
        return "H"
    return "I"

def build_kmer_store(id_to_pep: List[str], k: int) -> Tuple[array, array]:
    """
    Build (offsets, data) such that peptide pid has kmers:
      data[offsets[pid] : offsets[pid+1]]
    offsets is array('I') length n+1.
    data is array('H') for k<=3 else array('I') for k=4.
    """
    code = _choose_kmer_array_type(k)
    data = array(code)
    offsets = array("I", [0])

    for pep in id_to_pep:
        km = encode_kmers_unique_sorted(pep, k)
        data.extend(km)
        offsets.append(len(data))

    return offsets, data

def kmers_view(offsets: array, data: array, pid: int) -> memoryview:
    start = offsets[pid]
    end = offsets[pid + 1]
    return memoryview(data)[start:end]

def intersect_union_sizes(a: memoryview, b: memoryview) -> Tuple[int, int]:
    """
    Two-pointer merge on sorted unique integer lists.
    Returns (intersection_size, union_size).
    """
    i = j = 0
    inter = 0
    na = len(a)
    nb = len(b)
    while i < na and j < nb:
        av = a[i]
        bv = b[j]
        if av == bv:
            inter += 1
            i += 1
            j += 1
        elif av < bv:
            i += 1
        else:
            j += 1
    union = na + nb - inter
    return inter, union

def jaccard_sparse(a: memoryview, b: memoryview) -> float:
    if len(a) == 0 and len(b) == 0:
        return float("nan")
    inter, union = intersect_union_sizes(a, b)
    return (inter / union) if union else float("nan")


# ----------------------------
# LCS (unchanged; expensive; only for filtered candidates)
# ----------------------------

def lcs_length(a: str, b: str) -> int:
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return 0
    prev = [0] * (m + 1)
    best = 0
    for i in range(1, n + 1):
        cur = [0] * (m + 1)
        ai = a[i - 1]
        for j in range(1, m + 1):
            if ai == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best = cur[j]
        prev = cur
    return best

def lcs_fraction(a: str, b: str) -> float:
    denom = min(len(a), len(b))
    if denom == 0:
        return float("nan")
    return lcs_length(a, b) / denom


# ----------------------------
# Inverted index: kmer_id -> postings
# ----------------------------

def build_inverted_index(offsets: array, data: array, n: int) -> Dict[int, array]:
    idx: Dict[int, array] = {}
    for pid in range(n):
        start = offsets[pid]
        end = offsets[pid + 1]
        for km in data[start:end]:
            arr = idx.get(km)
            if arr is None:
                arr = array("I")
                idx[km] = arr
            arr.append(pid)
    return idx


# ----------------------------
# Globals for forked workers
# ----------------------------

G_ID_TO_PEP: List[str] = []
G_T_OFF: array
G_T_DATA: array
G_D_OFF: Optional[array] = None
G_D_DATA: Optional[array] = None
G_IDX: Dict[int, array] = {}

G_K = 3
G_MIN_SHARED = 3
G_JACC_T = 0.5
G_LCS_T = 0.5
G_MAX_CAND = 500
G_ALLOW_SELF = False
G_MAX_CONFLICTS_PER_T1 = 25
G_PRECOMP_DEC = True


def _init_worker(params: Tuple[int, int, float, float, int, bool, int, bool]) -> None:
    global G_K, G_MIN_SHARED, G_JACC_T, G_LCS_T, G_MAX_CAND, G_ALLOW_SELF, G_MAX_CONFLICTS_PER_T1, G_PRECOMP_DEC
    (G_K, G_MIN_SHARED, G_JACC_T, G_LCS_T, G_MAX_CAND, G_ALLOW_SELF, G_MAX_CONFLICTS_PER_T1, G_PRECOMP_DEC) = params


def _chunkify_ids(n_items: int, chunk_size: int) -> List[List[int]]:
    return [list(range(i, min(i + chunk_size, n_items))) for i in range(0, n_items, chunk_size)]


def _process_chunk_sparse(t1_ids: List[int]) -> Tuple[List[Tuple[int, int, float, float, int]], int]:
    out: List[Tuple[int, int, float, float, int]] = []

    for t1_id in t1_ids:
        t1 = G_ID_TO_PEP[t1_id]

        # Query kmers = reverse(t1) kmers
        if G_PRECOMP_DEC and (G_D_OFF is not None) and (G_D_DATA is not None):
            q_km = kmers_view(G_D_OFF, G_D_DATA, t1_id)
        else:
            # compute on the fly (slower, less RAM)
            q_list = encode_kmers_unique_sorted(t1[::-1], G_K)
            q_km = memoryview(array(_choose_kmer_array_type(G_K), q_list))

        # Candidate counting via postings
        counts: Dict[int, int] = {}
        getc = counts.get
        for km in q_km:
            post = G_IDX.get(int(km))
            if post is None:
                continue
            for t2_id in post:
                if (not G_ALLOW_SELF) and (t2_id == t1_id):
                    continue
                counts[t2_id] = getc(t2_id, 0) + 1

        if not counts:
            continue

        # Filter by min shared kmers
        # Avoid sorting enormous dicts: first keep only those passing min_shared.
        cand = [(pid, c) for pid, c in counts.items() if c >= G_MIN_SHARED]
        if not cand:
            continue

        # Deterministic ordering on the candidates we actually consider:
        # (-shared, len(t2), lex(t2))
        cand.sort(key=lambda pc: (-pc[1], len(G_ID_TO_PEP[pc[0]]), G_ID_TO_PEP[pc[0]]))
        if len(cand) > G_MAX_CAND:
            cand = cand[:G_MAX_CAND]

        best: List[Tuple[float, float, int, int]] = []  # (jac, lcsf, shared_kmers, t2_id)

        for t2_id, shared in cand:
            t2_km = kmers_view(G_T_OFF, G_T_DATA, t2_id)

            # Jaccard on sparse sets
            jac = jaccard_sparse(q_km, t2_km)
            if not math.isfinite(jac) or jac < G_JACC_T:
                continue

            lcsf = 0.0
            if G_LCS_T and G_LCS_T > 0:
                lcsf = lcs_fraction(t1[::-1], G_ID_TO_PEP[t2_id])
                if lcsf < G_LCS_T:
                    continue

            best.append((jac, lcsf, shared, t2_id))

        if not best:
            continue

        best.sort(key=lambda x: (-x[0], -x[1], -x[2], len(G_ID_TO_PEP[x[3]]), G_ID_TO_PEP[x[3]]))
        if len(best) > G_MAX_CONFLICTS_PER_T1:
            best = best[:G_MAX_CONFLICTS_PER_T1]

        for jac, lcsf, sk, t2_id in best:
            out.append((t1_id, t2_id, jac, lcsf, sk))

    return out, len(t1_ids)


# ----------------------------
# Conflict search entrypoints
# ----------------------------

def build_id_space(target_peps: Set[str]) -> List[str]:
    # mk3 behavior: lexical stable ordering
    return sorted(target_peps)

def find_conflicts_parallel(
    target_peps: Set[str],
    k: int,
    min_shared_kmers: int,
    jaccard_threshold: float,
    lcs_threshold: float,
    max_candidates: int,
    max_conflicts_per_peptide: int,
    allow_self: bool,
    workers: int,
    precompute_decoy_kmers: bool,
) -> List[Conflict]:
    global G_ID_TO_PEP, G_T_OFF, G_T_DATA, G_D_OFF, G_D_DATA, G_IDX

    if k <= 0:
        raise ValueError("k must be >= 1")
    if k > 4:
        raise ValueError("mk4 supports k <= 4 (20^k grows too large beyond that).")

    G_ID_TO_PEP = build_id_space(target_peps)
    n = len(G_ID_TO_PEP)

    # Build sparse k-mer stores
    G_T_OFF, G_T_DATA = build_kmer_store(G_ID_TO_PEP, k)
    if precompute_decoy_kmers:
        revs = [p[::-1] for p in G_ID_TO_PEP]
        G_D_OFF, G_D_DATA = build_kmer_store(revs, k)
    else:
        G_D_OFF, G_D_DATA = None, None

    # Inverted index over target kmers
    G_IDX = build_inverted_index(G_T_OFF, G_T_DATA, n)

    if workers <= 0:
        workers = 1

    chunk_size = max(2000, n // max(1, workers * 8))
    chunks = _chunkify_ids(n, chunk_size)

    ctx = mp.get_context("fork")
    params = (
        k, min_shared_kmers, jaccard_threshold, lcs_threshold, max_candidates,
        allow_self, max_conflicts_per_peptide, precompute_decoy_kmers
    )

    raw: List[Tuple[int, int, float, float, int]] = []

    with ctx.Pool(processes=workers, initializer=_init_worker, initargs=(params,)) as pool:
        it = pool.imap_unordered(_process_chunk_sparse, chunks, chunksize=1)
        pbar = tqdm(total=n, desc="Scanning peptides", unit="pep") if tqdm else None
        try:
            for out_chunk, n_done in it:
                raw.extend(out_chunk)
                if pbar:
                    pbar.update(n_done)
        finally:
            if pbar:
                pbar.close()

    conflicts: List[Conflict] = []
    for t1_id, t2_id, jac, lcsf, sk in raw:
        t1 = G_ID_TO_PEP[t1_id]
        t2 = G_ID_TO_PEP[t2_id]
        conflicts.append(
            Conflict(
                t1=t1,
                d1=t1[::-1],
                t2=t2,
                d2=t2[::-1],
                jaccard=jac,
                lcs_frac=lcsf,
                shared_kmers=sk,
            )
        )

    return conflicts


def find_conflicts(
    target_peps: Set[str],
    k: int,
    min_shared_kmers: int,
    jaccard_threshold: float,
    lcs_threshold: float,
    max_candidates: int,
    max_conflicts_per_peptide: int,
    allow_self: bool,
    workers: int,
    precompute_decoy_kmers: bool,
) -> List[Conflict]:
    if workers <= 1:
        return find_conflicts_parallel(
            target_peps=target_peps,
            k=k,
            min_shared_kmers=min_shared_kmers,
            jaccard_threshold=jaccard_threshold,
            lcs_threshold=lcs_threshold,
            max_candidates=max_candidates,
            max_conflicts_per_peptide=max_conflicts_per_peptide,
            allow_self=allow_self,
            workers=1,
            precompute_decoy_kmers=precompute_decoy_kmers,
        )
    return find_conflicts_parallel(
        target_peps=target_peps,
        k=k,
        min_shared_kmers=min_shared_kmers,
        jaccard_threshold=jaccard_threshold,
        lcs_threshold=lcs_threshold,
        max_candidates=max_candidates,
        max_conflicts_per_peptide=max_conflicts_per_peptide,
        allow_self=allow_self,
        workers=workers,
        precompute_decoy_kmers=precompute_decoy_kmers,
    )


def write_conflicts(conflicts: List[Conflict], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    conflicts_sorted = sorted(
        conflicts,
        key=lambda c: (c.t1, c.t2, -c.jaccard, -c.shared_kmers, -c.lcs_frac),
    )

    out_tsv = out_dir / "conflicts.tsv"
    with out_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["t1", "d1", "t2", "d2", "jaccard_kmers", "lcs_frac", "shared_kmers"])
        for c in conflicts_sorted:
            w.writerow([c.t1, c.d1, c.t2, c.d2, f"{c.jaccard:.6f}", f"{c.lcs_frac:.6f}", c.shared_kmers])

    per_t1: Dict[str, int] = {}
    for c in conflicts_sorted:
        per_t1[c.t1] = per_t1.get(c.t1, 0) + 1

    stats_rows = [
        ("conflict_edges_written", len(conflicts_sorted)),
        ("unique_t1_with_conflicts", len(per_t1)),
        ("max_conflicts_for_single_t1", max(per_t1.values()) if per_t1 else 0),
    ]
    out_stats = out_dir / "conflict_stats.tsv"
    _write_kv_tsv(stats_rows, out_stats)

    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {out_stats}")


# ----------------------------
# Pruning (unchanged high-level; note scalability)
# ----------------------------

def build_conflict_graph(conflicts: List[Conflict]) -> Dict[str, Set[str]]:
    g: Dict[str, Set[str]] = {}
    for c in conflicts:
        g.setdefault(c.t1, set()).add(c.t2)
        g.setdefault(c.t2, set()).add(c.t1)
    return g

def prune_targets(
    target_peps: Set[str],
    conflicts: List[Conflict],
    max_remove_fraction: float,
    prefer_remove_short: bool
) -> Tuple[Set[str], List[Tuple[str, str, int]]]:
    g = build_conflict_graph(conflicts)
    removed: Set[str] = set()
    plan: List[Tuple[str, str, int]] = []

    max_remove = int(math.floor(len(target_peps) * max_remove_fraction))
    if max_remove < 0:
        max_remove = 0

    def degree(p: str) -> int:
        return len(g.get(p, set()) - removed)

    def key(p: str) -> Tuple[int, int, str]:
        return (-degree(p), (len(p) if prefer_remove_short else 0), p)

    while True:
        candidates = [p for p in g.keys() if p not in removed and degree(p) > 0]
        if not candidates:
            break
        if len(plan) >= max_remove:
            break
        p = min(candidates, key=key)
        d = degree(p)
        removed.add(p)
        plan.append((p, "remove_high_conflict_degree", d))

    kept = set(target_peps) - removed
    return kept, plan

def write_removal_plan(plan: List[Tuple[str, str, int]], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "removal_plan.tsv"
    with out_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["peptide", "reason", "degree_at_removal"])
        for pep, reason, deg in plan:
            w.writerow([pep, reason, deg])
    print(f"Wrote: {out_tsv}")


# ----------------------------
# Main
# ----------------------------

def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    target_list = read_peptide_fasta(args.target_pep_fasta, args.min_len, args.max_len, args.ignore_il_equivalence)
    target_peps = set(target_list)
    del target_list  # critical for large datasets

    if args.cmd == "qc":
        run_qc(target_peps, args.out_dir, args.mass_tol_ppm, args.mass_tol_da)
        return 0

    if args.cmd == "conflicts":
        conflicts = find_conflicts(
            target_peps=target_peps,
            k=args.k,
            min_shared_kmers=args.min_shared_kmers,
            jaccard_threshold=args.jaccard_threshold,
            lcs_threshold=args.lcs_threshold,
            max_candidates=args.max_candidates,
            max_conflicts_per_peptide=args.max_conflicts_per_peptide,
            allow_self=args.self_match,
            workers=args.workers,
            precompute_decoy_kmers=args.precompute_decoy_kmers,
        )
        write_conflicts(conflicts, args.out_dir)
        run_qc(target_peps, args.out_dir, mass_tol_ppm=0.0, mass_tol_da=0.0)
        return 0

    if args.cmd == "prune":
        conflicts = find_conflicts(
            target_peps=target_peps,
            k=args.k,
            min_shared_kmers=args.min_shared_kmers,
            jaccard_threshold=args.jaccard_threshold,
            lcs_threshold=args.lcs_threshold,
            max_candidates=args.max_candidates,
            max_conflicts_per_peptide=10_000,
            allow_self=False,
            workers=args.workers,
            precompute_decoy_kmers=args.precompute_decoy_kmers,
        )
        write_conflicts(conflicts, args.out_dir)

        kept_targets, plan = prune_targets(
            target_peps=target_peps,
            conflicts=conflicts,
            max_remove_fraction=args.max_remove_fraction,
            prefer_remove_short=bool(args.prefer_remove_short)
        )
        write_removal_plan(plan, args.out_dir)

        filtered_target = args.out_dir / "filtered_target.fa"
        write_peptide_fasta(sorted(kept_targets), filtered_target, prefix="T_")
        print(f"Wrote: {filtered_target}")

        if args.write_decoy_fasta:
            filtered_decoy = args.out_dir / "filtered_decoy.fa"
            kept_decoys = sorted({p[::-1] for p in kept_targets})
            write_peptide_fasta(kept_decoys, filtered_decoy, prefix="D_")
            print(f"Wrote: {filtered_decoy}")

        qc_dir_before = args.out_dir / "qc_before"
        qc_dir_after = args.out_dir / "qc_after"
        run_qc(target_peps, qc_dir_before, mass_tol_ppm=0.0, mass_tol_da=0.0)
        run_qc(kept_targets, qc_dir_after, mass_tol_ppm=0.0, mass_tol_da=0.0)

        return 0

    raise RuntimeError("Unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
