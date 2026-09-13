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

# Filter clusters by coverage: map covered Proteins -> Clusters using the peptide table
def _filter_clusters_by_coverage(rolled_df: pd.DataFrame, pep_df: pd.DataFrame, coverage_tsv: Path) -> pd.DataFrame:
    if not coverage_tsv.is_file():
        raise SystemExit(f"[ERROR] grouped coverage file not found: {coverage_tsv}")
    cov = _read_table(coverage_tsv)
    allowed_clusters = set()
    if "Cluster" in cov.columns and cov["Cluster"].notna().any():
        allowed_clusters = set(pd.to_numeric(cov["Cluster"], errors="coerce").dropna().astype(int))
    elif "Protein" in cov.columns:
        covered_proteins = set(cov["Protein"].astype(str))
        prot2cluster = pep_df.dropna(subset=["Protein", "Cluster"])[["Protein","Cluster"]].copy()
        prot2cluster["Protein"] = prot2cluster["Protein"].astype(str)
        prot2cluster["Cluster"] = pd.to_numeric(prot2cluster["Cluster"], errors="coerce")
        prot2cluster = prot2cluster.dropna(subset=["Cluster"])
        allowed_clusters = set(prot2cluster[prot2cluster["Protein"].isin(covered_proteins)]["Cluster"].astype(int))
    else:
        return rolled_df
    if "Cluster" in rolled_df.columns:
        out = rolled_df[pd.to_numeric(rolled_df["Cluster"], errors="coerce").isin(allowed_clusters)].copy()
    else:
        out = rolled_df.copy()
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

# ---------- Core: cluster rollup from annotated peptide crosstab ----------
def rollup_clusters_from_annotated(
    pep_annot_tsv: Path,
    mode: str,
    use_flanked: bool = True,
    rollup: str = "sum",
    rrollup_summary: str | None = None,
    outlier_alpha: float | None = None,
) -> pd.DataFrame:
    df = _read_table(pep_annot_tsv)

    if "Cluster" not in df.columns:
        raise SystemExit("[ERROR] Input file lacks 'Cluster' column required for cluster-level rollup.")

    pep_key = _pick_peptide_key(df, prefer_flanked=use_flanked)
    have_plain = ("Peptide" in df.columns) and (pep_key != "Peptide")

    if "Unique" not in df.columns:
        df["Unique"] = False
    df["Unique"] = _coerce_bool_series(df["Unique"])

    # Identify sample columns (exclude flags/metadata)
    meta_exclude = {pep_key, "Protein", "Unique", "Cluster", "Gene", "gene", "Function", "function", "Unique to Cluster"}
    sample_cols = _infer_sample_cols(df, exclude=list(meta_exclude))

    # Slim working frame
    cols_needed = ["Cluster", "Unique", pep_key] + (["Protein"] if "Protein" in df.columns else []) + ["Unique to Cluster"] + sample_cols
    slim = df.reindex(columns=[c for c in cols_needed if c in df.columns]).copy()
    slim = slim.rename(columns={pep_key: "Peptide_pref"})
    if have_plain:
        slim["Peptide_plain"] = df["Peptide"].astype(str)
    else:
        slim["Peptide_plain"] = ""

    # Mode filter
    if mode == "unique_only":
        working = slim[slim["Unique"] == True].copy()
    elif mode == "requires_unique":
        keep_ids = slim.groupby("Cluster")["Unique"].any()
        working = slim[slim["Cluster"].isin(keep_ids[keep_ids].index)].copy()

    elif mode == "all_matches":
        working = slim.copy()
    else:
        raise SystemExit(f"[ERROR] Unsupported mode: {mode}")

    if working.empty:
        cols = ["Protein"] + sample_cols + ["Cluster","Gene","Function","Peptide Number",
                                            "Unique Peptide(s)","Shared Peptide(s)",
                                            "Unique to Cluster","Peptides","Flanked Peptides"]
        return pd.DataFrame(columns=cols)

    # --------- Non-RRollup rollups (sum/median/mean) ----------
    if rollup in ("sum", "median", "mean"):
        if rollup == "sum":
            agg_fn = {c: "sum" for c in sample_cols}
        elif rollup == "median":
            agg_fn = {c: "median" for c in sample_cols}
        else:
            agg_fn = {c: "mean" for c in sample_cols}

        grouped = working.groupby("Cluster", as_index=False).agg({
            **agg_fn,
            "Peptide_pref": lambda x: ";".join(sorted(set(map(str, x)))),
            "Peptide_plain": lambda x: ";".join(sorted({s for s in map(str, x) if s})),
            "Unique to Cluster": "any" if "Unique to Cluster" in working.columns else (lambda x: False),
        })

    # --------- RRollup ----------
    elif rollup == "rrollup":
        if str(rrollup_summary).lower() not in ("median", "mean"):
            rrollup_summary = "median"

        rows = []
        for clust, peptab in working.groupby("Cluster", sort=False):
            if peptab.shape[0] == 1:
                vals = peptab.iloc[0][sample_cols].values.astype(float)
                rows.append((clust, vals))
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
                clust_vals = np.nanmean(scaled, axis=0)
            else:
                clust_vals = np.nanmedian(scaled, axis=0)

            rows.append((clust, clust_vals))

        if not rows:
            grouped = pd.DataFrame(columns=["Cluster"] + sample_cols)
        else:
            grouped = pd.DataFrame({"Cluster": [c for c, _ in rows]})
            for j, c in enumerate(sample_cols):
                grouped[c] = [vals[j] for _, vals in rows]

        # carry UTC info by any()
        if "Unique to Cluster" in working.columns:
            utc_any = working.groupby("Cluster")["Unique to Cluster"].any().reset_index()
            grouped = grouped.merge(utc_any, on="Cluster", how="left")
        else:
            grouped["Unique to Cluster"] = False

        # Peptide lists for provenance
        peps_pref = working.groupby("Cluster")["Peptide_pref"].apply(lambda x: ";".join(sorted(set(map(str, x)))))
        peps_plain = working.groupby("Cluster")["Peptide_plain"].apply(lambda x: ";".join(sorted({s for s in map(str, x) if s})))
        grouped["Flanked Peptides"] = grouped["Cluster"].map(peps_pref).fillna("")
        grouped["Peptides"] = grouped["Cluster"].map(peps_plain).fillna("")

    else:
        raise SystemExit(f"[ERROR] Unsupported --rollup: {rollup}")

    # Presence flags (cluster-level) based on 'Unique' peptide presence
    uniq_map   = working.groupby("Cluster")["Unique"].any()
    shared_map = working.groupby("Cluster")["Unique"].apply(lambda s: (~s.astype(bool)).any())
    grouped["Unique Peptide(s)"] = grouped["Cluster"].map(uniq_map).fillna(False).astype(bool)
    grouped["Shared Peptide(s)"] = grouped["Cluster"].map(shared_map).fillna(False).astype(bool)

    # Rename peptide list columns if needed
    if "Flanked Peptides" not in grouped.columns:
        grouped = grouped.rename(columns={"Peptide_pref": "Flanked Peptides"})
    if "Peptides" not in grouped.columns:
        grouped = grouped.rename(columns={"Peptide_plain": "Peptides"})

    # Peptide Number
    base_list_col = "Flanked Peptides" if grouped["Flanked Peptides"].astype(str).str.len().gt(0).any() else "Peptides"
    grouped["Peptide Number"] = grouped[base_list_col].apply(
        lambda s: 0 if pd.isna(s) or str(s) == "" else len(set(str(s).split(";")))
    )

    # Synthesize 'Protein' for shape compatibility
    grouped["Protein"] = grouped["Cluster"].apply(lambda x: f"Cluster {int(x)}" if pd.notna(x) and str(x) != "" else "Cluster NA")

    # Leave Gene/Function BLANK by design
    grouped["Gene"] = ""
    grouped["Function"] = ""

    # Clean internal helper peptide cols if any
    for col in ["Peptide_pref", "Peptide_plain"]:
        if col in grouped.columns:
            grouped = grouped.drop(columns=[col])

    return grouped


def parse_args():
    p = argparse.ArgumentParser(description="Cluster rollup from annotated peptide crosstab.")
    p.add_argument("-i", "--input_tsv", required=True, dest="input_tsv",
                   help="Annotated peptide crosstab TSV (e.g., peptide_crosstab_annotated.tsv).")
    p.add_argument("-o", "--output_cluster_tsv", required=True, dest="output_cluster_tsv",
                   help="Path to write cluster-level crosstab (e.g., clustered_crosstab_annotated.tsv).")
    p.add_argument("--mode", choices=["unique_only", "all_matches", "requires_unique"],
                   default="all_matches", help="Peptide inclusion mode (default: all_matches).")
    p.add_argument("--rollup", choices=["sum","median","mean","rrollup"], default="sum",
                   help="Across-peptide rollup: sum/median/mean or rrollup (default: sum).")
    p.add_argument("--rrollup_summary", default=None,
                   help="For --rollup rrollup only: 'median' or 'mean' (default: median if omitted). Ignored otherwise.")
    p.add_argument("--outlier_alpha", type=float, default=None,
                   help="For --rollup rrollup: optional Grubbs alpha (e.g., 0.05). Omit to disable.")
    p.add_argument("--coverage_tsv", default=None,
                   help="Path to grouped_coverage.tsv. If not provided, uses <input_dir>/map_files/grouped_coverage.tsv.")
    p.add_argument("--use_flanked", choices=["True","False"], default="True",
                   help="Use PeptideFlanked as peptide identity if available (default True).")
    return p.parse_args()


def _finalize_cluster_column_order(df: pd.DataFrame) -> pd.DataFrame:
    """
    Final order identical to protein crosstab:
    Protein | <naturally sorted sample intensities> | Cluster | Gene | Function | Peptide Number | Unique Peptide(s) | Shared Peptide(s) | Unique to Cluster | Peptides | Flanked Peptides
    """
    import re as _re

    def _natural_key(s: str):
        return [int(t) if t.isdigit() else t.lower() for t in _re.split(r'(\d+)', s)]

    meta = {"Protein","Cluster","Gene","Function","Peptide Number","Unique Peptide(s)","Shared Peptide(s)",
            "Unique to Cluster","Peptides","Flanked Peptides"}
    samples = sorted(
        [c for c in df.columns if c not in meta and pd.api.types.is_numeric_dtype(df[c]) and not _is_boolish_col(df[c])],
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
    pep_df = _read_table(input_path)

    # 1) Cluster rollup
    rolled = rollup_clusters_from_annotated(
        input_path,
        mode=args.mode,
        use_flanked=use_flanked,
        rollup=args.rollup,
        rrollup_summary=args.rrollup_summary,
        outlier_alpha=args.outlier_alpha,
    )

    # 2) Filter to coverage-linked clusters using peptide map
    cov_path = Path(args.coverage_tsv).resolve() if args.coverage_tsv else (input_dir / "map_files" / "grouped_coverage.tsv")
    rolled = _filter_clusters_by_coverage(rolled, pep_df, cov_path)

    # 3) Finalize order and write cluster crosstab
    rolled = _finalize_cluster_column_order(rolled)
    out_path = Path(args.output_cluster_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rolled.to_csv(out_path, sep='\t', index=False)
    print(f"✅ Cluster crosstab written: {out_path}  (rows={rolled.shape[0]})")

if __name__ == "__main__":
    main()
