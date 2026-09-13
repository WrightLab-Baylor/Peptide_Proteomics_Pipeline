#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import pandas as pd
import numpy as np

# ---------- Helpers ----------
def _coerce_bool_series(s: pd.Series) -> pd.Series:
    return s.fillna(False).apply(lambda x: str(x).strip().lower() in ("true", "1", "t", "yes", "y"))

def _is_boolish_col(col: pd.Series) -> bool:
    # True for real bool dtype, 0/1 numeric columns, and common boolean-like strings
    if pd.api.types.is_bool_dtype(col):
        return True
    if pd.api.types.is_numeric_dtype(col) and set(pd.unique(col.dropna())) <= {0, 1}:
        return True
    uniq_vals = set(map(lambda v: str(v).strip().lower(), pd.unique(col.dropna())))
    return uniq_vals <= {"true","false","t","f","yes","no","y","n","0","1"}

def _infer_sample_cols(df: pd.DataFrame, exclude=()) -> list:
    exclude = set(exclude)
    numeric_like = []
    for c in df.columns:
        if c in exclude:
            continue
        # Skip boolean-like columns explicitly so they never get treated as intensities
        if _is_boolish_col(df[c]):
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            numeric_like.append(c)
    if not numeric_like:
        numeric_like = [c for c in df.columns if c not in exclude and not _is_boolish_col(df[c])]
    return numeric_like

def _pick_peptide_key(df: pd.DataFrame, prefer_flanked=True) -> str:
    if prefer_flanked and "PeptideFlanked" in df.columns:
        return "PeptideFlanked"
    elif "Peptide" in df.columns:
        return "Peptide"
    else:
        raise ValueError("Input must contain Peptide or PeptideFlanked")

def _read_table(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as e:
        raise SystemExit(f"[ERROR] Failed reading {path}: {e}")

def _filter_to_coverage_ids(rolled_df: pd.DataFrame, coverage_tsv: Path) -> pd.DataFrame:
    if not coverage_tsv.is_file():
        raise SystemExit(f"[ERROR] grouped coverage file not found: {coverage_tsv}")
    cov = _read_table(coverage_tsv)
    cov_ids = set(cov["Protein"].astype(str))
    out = rolled_df[rolled_df["Protein"].astype(str).isin(cov_ids)].copy()
    return out

# ---------- RRollup internals ----------
def _pick_reference(peptab: pd.DataFrame, sample_cols: list) -> int:
    miss = peptab[sample_cols].isna().sum(axis=1)
    totals = peptab[sample_cols].sum(axis=1, skipna=True)
    order = pd.DataFrame({'idx': peptab.index, 'miss': miss, 'tot': totals})
    return int(order.sort_values(['miss', 'tot'], ascending=[True, False]).iloc[0]['idx'])

def _median_ratio(ref_vec: np.ndarray, pep_vec: np.ndarray) -> float:
    mask = (~np.isnan(ref_vec)) & (~np.isnan(pep_vec)) & (pep_vec != 0)
    if not mask.any():
        return np.nan
    return float(np.median(ref_vec[mask] / pep_vec[mask]))

def _grubbs_filter(values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Two-sided, remove at most one outlier from a 1D array. Returns boolean mask of kept values."""
    x = values.astype(float)
    mask = ~np.isnan(x)
    n = int(mask.sum())
    if n < 3:
        return mask
    vals = x[mask]
    sd = np.nanstd(vals, ddof=1)
    if n < 3 or sd == 0:
        return mask
    z = np.abs(vals - np.nanmean(vals)) / sd
    i = int(np.argmax(z))
    G = float(z[i])
    from scipy import stats
    tcrit = stats.t.ppf(1 - alpha/(2*n), n - 2)
    Gcrit = ((n - 1)/np.sqrt(n)) * np.sqrt(tcrit**2 / (n - 2 + tcrit**2))
    if G > Gcrit:
        idx = np.flatnonzero(mask)[i]
        mask[idx] = False
    return mask

# ---------- Core: rollup from annotated peptide crosstab ----------
def rollup_from_annotated(
    pep_annot_tsv: Path,
    mode: str,
    use_flanked: bool = True,
    rollup: str = "sum",
    rrollup_summary: str | None = None,
    outlier_alpha: float | None = None,
) -> pd.DataFrame:
    """
    Roll up annotated peptide crosstab to protein level.

    mode : 'all_matches' | 'unique_only' | 'requires_unique'
    rollup : 'sum' | 'median' | 'mean' | 'rrollup'
    rrollup_summary : 'median' | 'mean' | None
      Only used if rollup == 'rrollup'. If None, defaults to 'median'.
    outlier_alpha : float or None (RRollup only)
    """
    df = _read_table(pep_annot_tsv)

    pep_key = _pick_peptide_key(df, prefer_flanked=use_flanked)
    have_plain = ("Peptide" in df.columns) and (pep_key != "Peptide")

    # Uniqueness flag ALWAYS from 'Unique' now (UTC is informational only)
    if "Unique" not in df.columns:
        df["Unique"] = False
    df["Unique"] = _coerce_bool_series(df["Unique"])

    # Annotations to carry
    gene_cols = [c for c in ("Gene", "gene") if c in df.columns]
    func_cols = [c for c in ("Function", "function") if c in df.columns]
    any_cols  = [c for c in ("Unique to Cluster",) if c in df.columns]  # informational only

    # Sample columns: exclude boolean-like
    meta_exclude = {pep_key, "Protein", "Unique", "Cluster"} | set(gene_cols) | set(func_cols) | set(any_cols)
    sample_cols = _infer_sample_cols(df, exclude=list(meta_exclude))

    # Slim working frame
    cols_needed = ["Protein", "Unique", pep_key] + (["Cluster"] if "Cluster" in df.columns else []) + gene_cols + func_cols + any_cols + sample_cols
    slim = df.reindex(columns=cols_needed).copy()
    slim = slim.rename(columns={pep_key: "Peptide_pref"})
    if have_plain:
        slim["Peptide_plain"] = df["Peptide"].astype(str)
    else:
        slim["Peptide_plain"] = ""

    # Apply mode filter
    if mode == "unique_only":
        working = slim[slim["Unique"] == True].copy()
    elif mode == "requires_unique":
        keep_ids = slim.groupby("Protein")["Unique"].any()
        working = slim[slim["Protein"].isin(keep_ids[keep_ids].index)].copy()

    elif mode == "all_matches":
        working = slim.copy()
    else:
        raise SystemExit(f"[ERROR] Unsupported mode: {mode}")

    if working.empty:
        # Build a minimal empty frame with expected columns
        return pd.DataFrame(columns=["Protein"] + sample_cols + ["Cluster","Gene","Function","Peptide Number",
                                                                 "Unique Peptide(s)","Shared Peptide(s)",
                                                                 "Unique to Cluster","Peptides","Flanked Peptides"])

    # Maps for annotation columns (carry through later)
    def _first_non_null(x):
        for v in x:
            if pd.notna(v) and str(v) != "":
                return v
        return ""
    def _any_boolish(x):
        return any(str(v).strip().lower() in ("true","1","t","yes","y") for v in x if pd.notna(v))

    cluster_map = working.groupby("Protein")["Cluster"].first() if "Cluster" in working.columns else pd.Series(dtype=object)
    gene_map    = {c: working.groupby("Protein")[c].apply(_first_non_null) for c in gene_cols}
    func_map    = {c: working.groupby("Protein")[c].apply(_first_non_null) for c in func_cols}
    any_maps    = {c: working.groupby("Protein")[c].apply(_any_boolish) for c in any_cols}

    # Presence flags (protein-level) based on 'Unique'
    uniq_map   = working.groupby("Protein")["Unique"].any()
    shared_map = working.groupby("Protein")["Unique"].apply(lambda s: (~s.astype(bool)).any())

    # --------- Non-RRollup rollups (sum/median/mean) ----------
    if rollup in ("sum", "median", "mean"):
        if rollup == "sum":
            agg_fn = {c: "sum" for c in sample_cols}
        elif rollup == "median":
            agg_fn = {c: "median" for c in sample_cols}
        else:  # mean
            agg_fn = {c: "mean" for c in sample_cols}

        grouped = working.groupby("Protein", as_index=False).agg({
            **agg_fn,
            "Peptide_pref": lambda x: ";".join(sorted(set(map(str, x)))),
            "Peptide_plain": lambda x: ";".join(sorted({s for s in map(str, x) if s})),
        })

        # Attach annotations
        if not cluster_map.empty:
            grouped["Cluster"] = grouped["Protein"].map(cluster_map)
        for c in gene_cols:
            grouped[c] = grouped["Protein"].map(gene_map[c])
        for c in func_cols:
            grouped[c] = grouped["Protein"].map(func_map[c])
        for c in any_cols:
            grouped[c] = grouped["Protein"].map(any_maps[c]).astype(bool)

        # Presence flags
        grouped["Unique Peptide(s)"] = grouped["Protein"].map(uniq_map).fillna(False).astype(bool)
        grouped["Shared Peptide(s)"] = grouped["Protein"].map(shared_map).fillna(False).astype(bool)

    # --------- RRollup ----------
    elif rollup == "rrollup":
        # Default rrollup summary if None/invalid
        if str(rrollup_summary).lower() not in ("median", "mean"):
            rrollup_summary = "median"

        prot_rows = []
        for prot, peptab in working.groupby("Protein", sort=False):
            if peptab.shape[0] == 1:
                vals = peptab.iloc[0][sample_cols].values.astype(float)
                prot_rows.append((prot, vals))
                continue

            ref_idx = _pick_reference(peptab, sample_cols)
            ref_row = peptab.loc[ref_idx]
            ref = ref_row[sample_cols].values.astype(float)

            scaled = []
            for _, row in peptab.iterrows():
                v = row[sample_cols].values.astype(float)
                if row.name == ref_idx:
                    v_scaled = v.copy()
                else:
                    sf = _median_ratio(ref, v)
                    v_scaled = v * sf if (np.isfinite(sf) and sf != 0) else np.full_like(v, np.nan)
                scaled.append(v_scaled)
            scaled = np.vstack(scaled)

            if outlier_alpha is not None:
                keep = np.ones_like(scaled, dtype=bool)
                for j in range(scaled.shape[1]):
                    kmask = _grubbs_filter(scaled[:, j], alpha=float(outlier_alpha))
                    keep[:, j] = kmask
                scaled = np.where(keep, scaled, np.nan)

            if str(rrollup_summary).lower() == "mean":
                prot_vals = np.nanmean(scaled, axis=0)
            else:
                prot_vals = np.nanmedian(scaled, axis=0)

            prot_rows.append((prot, prot_vals))

        if not prot_rows:
            grouped = pd.DataFrame(columns=["Protein"] + sample_cols)
        else:
            grouped = pd.DataFrame({"Protein": [p for p, _ in prot_rows]})
            for j, c in enumerate(sample_cols):
                grouped[c] = [vals[j] for _, vals in prot_rows]

        # Attach annotations & flags
        if not cluster_map.empty:
            grouped["Cluster"] = grouped["Protein"].map(cluster_map)
        for c in gene_cols:
            grouped[c] = grouped["Protein"].map(gene_map[c])
        for c in func_cols:
            grouped[c] = grouped["Protein"].map(func_map[c])
        for c in any_cols:
            grouped[c] = grouped["Protein"].map(any_maps[c]).astype(bool)

        grouped["Unique Peptide(s)"] = grouped["Protein"].map(uniq_map).fillna(False).astype(bool)
        grouped["Shared Peptide(s)"] = grouped["Protein"].map(shared_map).fillna(False).astype(bool)

        # Peptide lists for provenance
        peps_pref = working.groupby("Protein")["Peptide_pref"].apply(lambda x: ";".join(sorted(set(map(str, x)))))
        peps_plain = working.groupby("Protein")["Peptide_plain"].apply(lambda x: ";".join(sorted({s for s in map(str, x) if s})))
        grouped["Flanked Peptides"] = grouped["Protein"].map(peps_pref).fillna("")
        grouped["Peptides"] = grouped["Protein"].map(peps_plain).fillna("")

    else:
        raise SystemExit(f"[ERROR] Unsupported --rollup: {rollup}")

    # Normalize peptide list columns & counts (works for all rollups)
    if "Flanked Peptides" not in grouped.columns:
        grouped["Flanked Peptides"] = grouped["Protein"].map(
            working.groupby("Protein")["Peptide_pref"].apply(lambda x: ";".join(sorted(set(map(str, x)))))
        ).fillna("")
    if "Peptides" not in grouped.columns:
        grouped["Peptides"] = grouped["Protein"].map(
            working.groupby("Protein")["Peptide_plain"].apply(lambda x: ";".join(sorted({s for s in map(str, x) if s})))
        ).fillna("")

    base_list_col = "Flanked Peptides" if grouped["Flanked Peptides"].astype(str).str.len().gt(0).any() else "Peptides"
    grouped["Peptide Number"] = grouped[base_list_col].apply(
        lambda s: 0 if pd.isna(s) or str(s) == "" else len(set(str(s).split(";")))
    )

    # === Clean up: drop internal peptide columns so they don't appear at the end ===
    for col in ["Peptide_pref", "Peptide_plain"]:
        if col in grouped.columns:
            grouped = grouped.drop(columns=[col])

    return grouped


def parse_args():
    p = argparse.ArgumentParser(description="Protein rollup from annotated peptide crosstab.")
    p.add_argument("-i", "--input_tsv", required=True, dest="input_tsv",
                   help="Annotated peptide crosstab TSV (e.g., peptide_crosstab_annotated.tsv).")
    p.add_argument("-o", "--output_protein_tsv", required=True, dest="output_protein_tsv",
                   help="Path to write protein-level crosstab (e.g., protein_crosstab_annotated.tsv).")
    p.add_argument("--mode", choices=["unique_only", "all_matches", "requires_unique"],
                   default="all_matches", help="Peptide inclusion mode (default: all_matches).")
    p.add_argument("--rollup", choices=["sum","median","mean","rrollup"], default="sum",
                   help="Across-peptide rollup: sum/median/mean or rrollup (default: sum).")
    # Accept any string for rrollup_summary; shell may pass 'None' when rollup != rrollup
    p.add_argument("--rrollup_summary", default=None,
                   help="For --rollup rrollup only: 'median' or 'mean' (default: median if omitted). Ignored otherwise.")
    p.add_argument("--outlier_alpha", type=float, default=None,
                   help="For --rollup rrollup: optional Grubbs alpha (e.g., 0.05). Omit to disable.")
    p.add_argument("--coverage_tsv", default=None,
                   help="Path to grouped_coverage.tsv. If not provided, uses <input_dir>/map_files/grouped_coverage.tsv.")
    p.add_argument("--use_flanked", choices=["True","False"], default="True",
                   help="Use PeptideFlanked as peptide identity if available (default True).")
    return p.parse_args()


def _finalize_protein_column_order(df: pd.DataFrame) -> pd.DataFrame:
    """
    Final order:
    Protein | <naturally sorted sample intensities> | Cluster | Gene | Function | Peptide Number | Unique Peptide(s) | Shared Peptide(s) | Unique to Cluster | Peptides | Flanked Peptides
    Any other columns are appended last.
    """
    import re as _re

    def _natural_key(s: str):
        return [int(t) if t.isdigit() else t.lower() for t in _re.split(r'(\d+)', s)]

    meta = {"Protein","Cluster","Gene","Function","Peptide Number","Unique Peptide(s)","Shared Peptide(s)",
            "Unique to Cluster","Peptides","Flanked Peptides","Peptide_pref","Peptide_plain"}
    samples = sorted(
        [c for c in df.columns
         if c not in meta and pd.api.types.is_numeric_dtype(df[c]) and not _is_boolish_col(df[c])],
        key=_natural_key
    )

    desired = (["Protein"] + samples + ["Cluster","Gene","Function","Peptide Number",
               "Unique Peptide(s)","Shared Peptide(s)","Unique to Cluster","Peptides","Flanked Peptides"])

    final = [c for c in desired if c in df.columns]
    seen = set(final)
    final += [c for c in df.columns if c not in seen]
    return df[final]

def main():
    args = parse_args()
    use_flanked = (args.use_flanked == "True")

    input_path = Path(args.input_tsv).resolve()
    input_dir = input_path.parent

    # 1) Protein rollup from annotated peptide crosstab (RRollup optional)
    rolled = rollup_from_annotated(
        input_path,
        mode=args.mode,
        use_flanked=use_flanked,
        rollup=args.rollup,
        rrollup_summary=args.rrollup_summary,
        outlier_alpha=args.outlier_alpha,
    )

    # 2) Always filter to coverage IDs
    cov_path = Path(args.coverage_tsv).resolve() if args.coverage_tsv else (input_dir / "map_files" / "grouped_coverage.tsv")
    rolled = _filter_to_coverage_ids(rolled, cov_path)

    # 3) Finalize order and write protein crosstab
    rolled = _finalize_protein_column_order(rolled)
    out_path = Path(args.output_protein_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rolled.to_csv(out_path, sep='\t', index=False)
    print(f"✅ Protein crosstab written: {out_path}  (rows={rolled.shape[0]})")

    # 4) (Clustered crosstab intentionally omitted.)

if __name__ == "__main__":
    main()
