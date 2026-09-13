# Binning

This directory contains the peptide-binning FDR workflow used to measure how
database partitioning changes peptide-spectrum-match confidence and to derive
empirical FDR corrections that can be applied back to the original
full-database search.

The workflow preserves the complete full-database PSM assignments during final
correction. Partition searches are used to estimate systematic FDR changes, not
to replace the original search results.

## Scope

The Binning arm begins with **existing partition-search outputs** and a matching
full-database search. It does not build peptide bins or run the upstream
RAW-to-MS-GF+ search itself.

For each partition analysis it can:

- prepare partition SIC tables;
- retain full-database-matching scan/peptide assignments;
- estimate partition-level FDR;
- compare partitioned and full-database FDR values;
- derive and apply an empirical correction curve;
- regenerate peptide and protein rollups;
- build matched-formula outputs across numbered samples; and
- collate compact scan-recovery, FDR, peptide, and protein summaries.

## Requirements

A Conda environment is provided:

```bash
conda env create -f environment.yml
conda activate peptide-binning-fdr
```

The Python workflow uses NumPy, pandas, SciPy, Matplotlib, and Seaborn. The shell
launchers also expect standard Unix tools such as `bash`, `find`, `sed`, `awk`,
`sort`, and `realpath`.

## Expected data layout

Repository-local runs normally live under `Binning/data/`:

```text
Binning/
├── data/
│   ├── MyCondition_1/
│   │   ├── MyCondition_Binning/
│   │   │   ├── ..._group_1/
│   │   │   ├── ..._group_2/
│   │   │   ├── ..._pair_1_2/
│   │   │   └── ...
│   │   ├── MyCondition_Groups_1.2M/
│   │   ├── MyCondition_Groups_600k/
│   │   └── ...
│   ├── MyCondition_2/
│   ├── MyCondition_Full/
│   │   ├── fdr_esti/
│   │   │   └── *_withsyn.tsv
│   │   ├── *.fasta
│   │   ├── peptide_crosstab_annotated.tsv
│   │   └── protein_crosstab_annotated.tsv
│   ├── Matched_Formula_Output/
│   └── results/
├── scripts/
└── docs/
```

Each partition directory is expected to contain the upstream files required by
the selected stages, including `SICdir` and the partition PHRP outputs under
`results/PHRPOut`.

Large datasets do not need to live inside the repository. The launchers accept
explicit paths, and the mass matched runner/collator accept an alternate data
root.

## Main workflow

The primary launcher is:

```bash
bash scripts/pipeline/run_binning_pipeline.sh \
    --dataset-dir data/MyCondition_1/MyCondition_Binning \
    --reference-withsyn data/MyCondition_Full/fdr_esti/sample_withsyn.tsv \
    --full-dir data/MyCondition_Full
```

The launcher automatically detects Individual (`*_group_*`) and Pairwise
(`*_pair_*_*`) partitions unless told otherwise. Use `--dry-run` to preview the
workflow without writing outputs.

For an Individual-bin-size series, run the same launcher against each
`MyCondition_Groups_<LABEL>` directory. The collator discovers these directories
dynamically; labels such as `1.2M`, `500k`, or custom experiment labels are all
supported. Familiar paper labels are only given a preferred display order.

See:

```bash
bash scripts/pipeline/run_binning_pipeline.sh --help
```

for stage controls and compatibility options.

## Matched-formula workflow

When several numbered runs represent different samples from the same condition,
the recommended matched workflow can apply each numbered-run formula only to
its corresponding source sample and combine the corrected samples for one
downstream rollup.

With the standard `Binning/data` layout:

```bash
bash scripts/pipeline/run_all_matched_formula_outputs.sh
```

Use `--dry-run` first if you want to inspect the discovered condition/design
matches. An external data root can be supplied with `--root PATH`.

The lower-level matched runner is
`scripts/pipeline/run_matched_correction_downstream.sh`.

## Collating results

To build compact result tables for every discovered numbered run:

```bash
python3 scripts/analysis/collate_binning_results.py --include-fdr
```

To collate one run only:

```bash
python3 scripts/analysis/collate_binning_results.py \
    --run MyCondition_1 \
    --include-fdr \
    --strict
```

The collator defaults to `Binning/data`, writes per-run outputs under
`data/results/<Run>/`, and automatically includes matched-formula summaries when
`data/Matched_Formula_Output/` is present. Use `--skip-matched` when you want to
collate only the numbered-run outputs.

The authoritative field-level description of these tables is:

**`docs/Results Guide.tsv`**

## Counting rules

The current collator uses biological identifiers rather than raw crosstab row
counts:

- peptides: distinct non-empty `Peptide`;
- `PeptideFlanked`: distinct non-empty values retained as a QC metric;
- proteins: distinct non-empty `Protein`;
- source-row counts: retained only for provenance/QC.

A peptide expanded across multiple protein-mapping rows is therefore counted
once as a peptide.

## Terminology and compatibility names

Public terminology:

- **Individual bins** — single peptide partitions;
- **Pairwise bins** — combinations of two Individual bins.

Validated on-disk names such as `*_group_<N>`, `Group_Downstream`,
`*_pair_<N>_<M>`, and `Pairwise_Downstream` are intentionally retained for
compatibility.

The collated `AnalysisFamily` value `Pairwise` intentionally denotes the primary
Individual-vs-Pairwise comparison. `size_series` denotes the Individual-bin-size
series.

## Directory guide

- `scripts/pipeline/` — user-facing orchestration shells.
- `scripts/binning/` — SIC preparation, classical-guided culling, and partition
  FDR collection.
- `scripts/binning/_internal/` — required implementation engines used by
  `classical_guided_cull.py`; normally do not invoke these directly.
- `scripts/fdr/` — FDR estimation and protein-mapping merge utilities.
- `scripts/correction/` — correction-curve derivation and application.
- `scripts/downstream/` — peptide filtering, mapping, and protein/cluster rollup.
- `scripts/analysis/` — production comparison, plotting, and collation tools.
- `docs/Workflow.md` — stage-by-stage workflow overview.
- `docs/Results Guide.tsv` — collated-output data dictionary.
- `examples/` — example layout and command notes; no large example data are
  bundled.

## Reproducibility notes

The current release preserves historical machine-facing directory names when
changing them would alter validated interfaces. User-facing descriptions use the
newer Individual/Pairwise terminology.

For maintainers or coding assistants, see `AI_README.md` before changing paths,
counting rules, correction semantics, or output schemas.
