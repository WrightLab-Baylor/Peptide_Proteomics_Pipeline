# Examples

No large example datasets are bundled with the repository. This folder documents
the expected layout and provides portable command examples.

## Minimal layout

```text
Binning/
├── data/
│   ├── MyCondition_1/
│   │   └── MyCondition_Binning/
│   │       ├── sample_group_1/
│   │       ├── sample_group_2/
│   │       └── sample_pair_1_2/
│   └── MyCondition_Full/
│       ├── fdr_esti/
│       │   └── sample_withsyn.tsv
│       └── database.fasta
└── scripts/
```

Each partition directory must already contain the upstream SIC/PHRP products
required by the stages you plan to run.

## Preview a run

From the `Binning` directory:

```bash
bash scripts/pipeline/run_binning_pipeline.sh \
    --dataset-dir data/MyCondition_1/MyCondition_Binning \
    --reference-withsyn data/MyCondition_Full/fdr_esti/sample_withsyn.tsv \
    --full-dir data/MyCondition_Full \
    --dry-run
```

Remove `--dry-run` after the discovered paths and stages look correct.

## Collate results

```bash
python3 scripts/analysis/collate_binning_results.py \
    --run MyCondition_1 \
    --include-fdr
```

For full output definitions, see `../docs/Results Guide.tsv`.

Small, de-identified example outputs may be added here in the future, but this
directory should not contain large raw or search-result datasets.
