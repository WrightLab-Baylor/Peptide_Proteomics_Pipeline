# Binning workflow

The Binning arm starts from existing partition-search outputs plus a matching
full-database search. It does **not** create the peptide partitions or run the
upstream MS-GF+ searches.

The main launcher, `scripts/pipeline/run_binning_pipeline.sh`, performs the
following stages for Individual-bin (`*_group_<N>`) and/or Pairwise-bin
(`*_pair_<N>_<M>`) directories:

1. Prepare SIC tables from partition-level `SICdir` inputs.
2. Cull partition SIC assignments against a full-database `*_withsyn.tsv`
   reference, with the validated SYN fallback behavior.
3. Estimate partition-level target-decoy FDR.
4. Merge protein mappings from PHRP SYN outputs.
5. Collect partition FDR values and audit scan/peptide conflicts.
6. Compare full-database and partitioned FDR values.
7. Summarize and plot FDR changes across the full-database FDR distribution.
8. Optionally derive an empirical correction curve, apply it to the full
   database results, and regenerate peptide/protein rollups.

For conditions represented by several numbered runs, the matched-formula
workflow can then apply each run's formula only to the source sample from which
that formula was derived:

```bash
bash scripts/pipeline/run_all_matched_formula_outputs.sh
```

By default, the mass matched runner and the result collator use `Binning/data`
as the data root. Explicit path overrides remain available for HPC or external
storage.

Finally, collate compact result tables with:

```bash
python3 scripts/analysis/collate_binning_results.py --include-fdr
```

See `../README.md` for setup and command examples, and `Results Guide.tsv` for
the definitions of collated output files and columns.

## Compatibility terminology

Public documentation uses **Individual bins** and **Pairwise bins**. Historical
on-disk names such as `Group_Downstream`, `*_group_<N>`, and
`Pairwise_Downstream` are intentionally retained because they are part of the
validated file interface.

Likewise, the collated `AnalysisFamily` value `Pairwise` is intentionally used
for the primary Individual-vs-Pairwise comparison, while `size_series` denotes
the Individual-bin-size analysis.


Size-series datasets are discovered dynamically from `<Condition>_Groups_<LABEL>` directories; labels are not restricted to the manuscript size series.
