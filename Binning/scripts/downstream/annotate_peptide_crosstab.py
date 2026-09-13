#!/usr/bin/env python3
import argparse
import sys
import re
from pathlib import Path
import pandas as pd

def str2bool(v):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("yes", "true", "t", "1", "y"):
        return True
    if s in ("no", "false", "f", "0", "n"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected (True/False).")

def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate peptide–protein crosstab with Cluster/Gene/Function (from FASTA headers) and uniqueness metrics by PeptideFlanked."
    )
    p.add_argument("-i", "--crosstab-tsv", required=True, help="Input peptide–protein crosstab TSV.")
    p.add_argument("-f", "--fasta", required=True, help="Protein FASTA with Cluster/Gene/Function in headers.")
    p.add_argument("-o", "--output-tsv", required=True, help="Output annotated TSV path.")
    # Column names in the crosstab (keep defaults aligned with your generator):
    p.add_argument("--peptide-col", default="Peptide", help="Core peptide column (default: Peptide).")
    p.add_argument("--peptide-flanked-col", default="PeptideFlanked", help="Flanked peptide column (default: PeptideFlanked).")
    p.add_argument("--protein-col", default="Protein", help="Protein column (default: Protein).")
    # Output column names (fixed by your spec):
    p.add_argument("--cluster-col", default="Cluster", help="Name of the Cluster column to create (default: Cluster).")
    p.add_argument("--cluster-tag", default="Cluster=", help="FASTA header tag for cluster (default: Cluster=).")
    # Filtering option: Unique to Cluster
    p.add_argument(
    "-unique_to_cluster", "--unique-to-cluster",
    type=str2bool, default=False,
    help="True → keep only rows where 'Unique to Cluster' is True. False (default) → annotate only.")
    return p.parse_args()

# ---- Helpers ----
def get_gene_from_info(info: str) -> str:
    """
    Extract a gene symbol from the trailing info string (e.g., 'OS=... GN=TP53 PE=1 ...').
    Looks for GN= or Gene= (case-insensitive). Returns 'Unknown' if not found.
    """
    if not isinstance(info, str) or not info:
        return "Unknown"
    m = re.search(r'(?:^|\s|;|\|)(?:GN|Gene)=([^\s;|]+)', info, flags=re.IGNORECASE)
    return m.group(1) if m else "Unknown"

def _remove_cluster_tokens(text: str) -> str:
    """Strip any 'Cluster=...' tokens from a text blob, collapse whitespace, and trim separators."""
    if not isinstance(text, str):
        return ""
    cleaned = re.sub(r'(?:^|\s|;|\|)Cluster=[^\s;|]+', ' ', text, flags=re.IGNORECASE)
    cleaned = cleaned.strip(" ;|")
    cleaned = " ".join(cleaned.split())
    return cleaned

def get_protein_function(protein_id: str, proteins: dict) -> dict:
    """
    Given a protein_id and a proteins dict like {id: {'Function': <full header tail>}},
    split at 'OS' to separate the human-readable function from the info tail,
    remove any 'Cluster=*' tokens from function, then parse gene from the info tail.
    """
    protein_info = proteins.get(protein_id)
    if protein_info:
        function_string = protein_info.get('Function', '') or ''
        if "OS" in function_string:
            function_parts = function_string.split('OS', 1)  # first occurrence
            function = function_parts[0].strip()
            info = "OS" + function_parts[1].strip() if len(function_parts) > 1 else ''
        else:
            function = function_string.strip()
            info = ''
        # Remove any embedded Cluster= tags from the function text
        function = _remove_cluster_tokens(function)
        if not function:
            function = "Unknown"
        gene = get_gene_from_info(info)
        if not gene:
            gene = "Unknown"
        return {'function': function, 'gene': gene}
    else:
        return {'function': 'Unknown', 'gene': 'Unknown'}

# ---- FASTA header parser ----
def build_protein_info_map(fasta_path: Path, cluster_tag: str = "Cluster=") -> dict:
    """
    Parse FASTA headers and return:
      { protein_id -> {"Cluster": <str or None>, "Function": <full header tail string>} }

    Conventions:
      - protein_id = token up to first whitespace (e.g., 'tr|F9UP40|F9UP40_LACPL')
      - Cluster parsed via `cluster_tag` (default 'Cluster='); case-insensitive
      - Function kept as entire header tail (after protein_id), including OS=, GN=, etc.
    """
    prot2info = {}
    cluster_re = re.compile(r'(?:^|\s|;|\|)' + re.escape(cluster_tag) + r'([^\s;|]+)', flags=re.IGNORECASE)

    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n")
            protein_id = header.split()[0]            # accession-like token up to first whitespace
            rest = header[len(protein_id):].lstrip()  # full header tail (may include OS=, GN=, PE=, SV=, Cluster=, etc.)
            m = cluster_re.search(header)
            clust = m.group(1) if m else None
            prot2info[protein_id] = {"Cluster": clust, "Function": rest}

    return prot2info

def main():
    args = parse_args()
    in_path = Path(args.crosstab_tsv)
    fa_path = Path(args.fasta)
    out_path = Path(args.output_tsv)

    if not in_path.is_file():
        print(f"ERROR: Input crosstab not found: {in_path}", file=sys.stderr)
        sys.exit(1)
    if not fa_path.is_file():
        print(f"ERROR: FASTA not found: {fa_path}", file=sys.stderr)
        sys.exit(1)

    # Load crosstab (strings first; coerce sample cols later)
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    required = [args.peptide_col, args.peptide_flanked_col, args.protein_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: Missing required columns {missing} in {in_path}", file=sys.stderr)
        sys.exit(1)

    # Normalize keys
    for c in required:
        df[c] = df[c].astype(str).str.strip()

    # Identify sample columns (everything not reserved)
    reserved = set(required) | {args.cluster_col, "Unique", "Shared", "Unique to Cluster", "Gene", "Function"}
    sample_cols = [c for c in df.columns if c not in reserved]

    # Build Protein -> (Cluster, Function-tail) map
    prot2info = build_protein_info_map(fa_path, cluster_tag=args.cluster_tag)

    # Map Cluster directly
    df[args.cluster_col] = df[args.protein_col].map(lambda x: prot2info.get(str(x), {}).get("Cluster"))

    # Derive Function and Gene via helper
    info_series = df[args.protein_col].apply(lambda pid: get_protein_function(str(pid), prot2info))
    df["Function"] = info_series.apply(lambda d: d["function"])
    df["Gene"]     = info_series.apply(lambda d: d["gene"])

    # ---- Uniqueness annotations (by PeptideFlanked) ----
    # Unique: PeptideFlanked maps to exactly one Protein
    nprot = (
        df.groupby(args.peptide_flanked_col)[args.protein_col]
          .nunique(dropna=False)
          .rename("_nprot_")
    )
    df = df.merge(nprot, left_on=args.peptide_flanked_col, right_index=True, how="left")
    df["Unique"] = df["_nprot_"].eq(1)
    df.drop(columns=["_nprot_"], inplace=True)

    # Unique to Cluster: PeptideFlanked maps to exactly one Cluster
    nclus = (
        df.groupby(args.peptide_flanked_col)[args.cluster_col]
          .nunique(dropna=False)
          .rename("_nclus_")
    )
    df = df.merge(nclus, left_on=args.peptide_flanked_col, right_index=True, how="left")
    df["Unique to Cluster"] = df["_nclus_"].eq(1)
    df.drop(columns=["_nclus_"], inplace=True)

    # Shared = not Unique
    df["Shared"] = ~df["Unique"]

    # Optional gate: keep only peptides unique to their cluster
    if args.unique_to_cluster:
        df = df[df["Unique to Cluster"]].copy()

    # Coerce sample columns to numeric (FutureWarning-safe: no errors="ignore")
    for c in sample_cols:
        try:
            df[c] = pd.to_numeric(df[c])
        except Exception:
            # Non-numeric column — leave as-is
            pass

    # Final column order:
    # Peptide, PeptideFlanked, Protein, (Samples...), Cluster, Unique to Cluster, Unique, Shared, Gene, Function
    front = [args.peptide_col, args.peptide_flanked_col, args.protein_col]
    tail  = [args.cluster_col, "Unique to Cluster", "Unique", "Shared", "Gene", "Function"]
    middle = [c for c in df.columns if c not in set(front + tail)]
    df = df[front + middle + tail]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote annotations → {out_path}  (rows: {len(df)}, samples: {len(sample_cols)})")

if __name__ == "__main__":
    main()
