# AI Helper / Maintainer Guide — Binning

Repository: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline
Workflow directory: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline/tree/main/Binning

Related workflow directories:

- FASTA Processing: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline/tree/main/FASTA_Processing
- Targeted Peptide Removal: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline/tree/main/Targeted_Peptide_Removal

**FASTA Processing is step 1** unless the required peptide FASTA and Binning
partition databases already exist. Before inventing missing preprocessing logic,
inspect `FASTA_Processing/README.md` and `FASTA_Processing/AI_README.md`.

This file is an anchor point for an AI assistant helping a user **run, understand,
troubleshoot, or maintain** the Binning workflow.

The intended user may be a new graduate student with little or no command-line
experience. Prefer concrete commands, explain what each command will do before it
is run, and preserve the validated scientific behavior described below.

---

## 1. First priority: help the user run Binning safely

When a user asks how to run Binning, do not begin by modifying code. First help
them identify:

1. the Binning software directory;
2. the Binning data root;
3. the numbered run they want to process, such as `Soil_1`;
4. the matching full-database run, such as `Soil_Full`;
5. the partitioned dataset directory, such as `Soil_1/Soil_Binning`;
6. the matching full-database `*_withsyn.tsv` reference.

The standard repository-local data root is:

```text
Binning/data/
```

A typical condition looks like:

```text
Binning/
├── scripts/
├── docs/
├── examples/
└── data/
    ├── Soil_1/
    │   └── Soil_Binning/
    ├── Soil_2/
    │   └── Soil_Binning/
    ├── ...
    └── Soil_Full/
        └── fdr_esti/
            └── <sample>_withsyn.tsv
```

For large HPC datasets, the data may live elsewhere. In that case, use explicit
paths instead of moving the data into the repository.

### For a novice user

Prefer to:

- give one command at a time;
- explain what directory they should be in;
- use absolute paths when there is any ambiguity;
- recommend `--dry-run` before destructive or expensive operations;
- avoid changing scripts merely to accommodate a local path;
- never assume a dataset or reference file when multiple candidates exist.

If the user is unsure about the directory structure, ask them to show a shallow
directory listing such as:

```bash
find data -maxdepth 3 -type d | sort
```

or, for an external data root:

```bash
find /path/to/binning_data -maxdepth 3 -type d | sort
```

Do not request a huge recursive listing unless needed.

---

## 2. Main user-facing workflow

The primary launcher is:

```text
scripts/pipeline/run_binning_pipeline.sh
```

It intentionally requires an explicit partitioned dataset directory and an
explicit full-database `*_withsyn.tsv` reference.

From the Binning software root, a dry run looks like:

```bash
bash scripts/pipeline/run_binning_pipeline.sh \
    --dataset-dir data/Soil_1/Soil_Binning \
    --reference-withsyn data/Soil_Full/fdr_esti/<sample>_withsyn.tsv \
    --dry-run
```

When the user also wants the correctional-formula stage, provide the matching
full-database run directory:

```bash
bash scripts/pipeline/run_binning_pipeline.sh \
    --dataset-dir data/Soil_1/Soil_Binning \
    --reference-withsyn data/Soil_Full/fdr_esti/<sample>_withsyn.tsv \
    --full-dir data/Soil_Full
```

Do not invent `<sample>`. Resolve the actual filename from the user's directory.

### What the main launcher does

The launcher coordinates:

1. SIC preparation;
2. classical-guided culling;
3. partition-level FDR estimation;
4. protein-mapping merge;
5. partition FDR collection and recovery audit;
6. full-database versus partitioned FDR comparison;
7. FDR-delta summary/plot generation;
8. optional correction-formula derivation, application, and downstream
   peptide/protein rollup.

The launcher can operate on Individual bins, Pairwise bins, or both, depending
on the partition folders present and the selected options.

When explaining failures, identify **which stage failed** before recommending a
fix.

---

## 3. Classical-guided culling

The public culling interface is:

```text
scripts/binning/classical_guided_cull.py
```

The files under:

```text
scripts/binning/_internal/
```

are required implementation engines. They are not obsolete duplicates and
normally should not be invoked directly.

User-facing terminology is:

- **Individual bins**
- **Pairwise bins**

Validated machine-facing names remain:

- `*_group_<N>`
- `*_pair_<N>_<M>`
- `Group_Downstream`
- `Pairwise_Downstream`

Do not rename those paths casually.

---

## 4. Correctional-formula workflows

### One analysis

The correctional downstream runner is:

```text
scripts/pipeline/run_correctional_formula_downstream.sh
```

It derives an empirical correction curve from the partition/full-database
comparison, applies it to the original full-database results, and runs the
standard peptide/protein downstream workflow.

The source full-database directory is treated as read-only.

### Matched multi-sample application

For a multi-sample condition, each numbered-series formula can be applied only
to the sample from which that formula was derived.

The single-analysis matched runner is:

```text
scripts/pipeline/run_matched_correction_downstream.sh
```

The mass-discovery wrapper is:

```text
scripts/pipeline/run_all_matched_formula_outputs.sh
```

Both default to:

```text
<Binning>/data
```

unless a different root is provided.

A useful first command is:

```bash
bash scripts/pipeline/run_all_matched_formula_outputs.sh --dry-run
```

The expected matched output root is:

```text
data/Matched_Formula_Output/
```

Do not describe matched application as applying one sample's formula globally
to unrelated samples.

---

## 5. Collating final results

The primary result collator is:

```text
scripts/analysis/collate_binning_results.py
```

It defaults to:

```text
<Binning>/data
```

and writes to:

```text
<Binning>/data/results/
```

For one numbered run:

```bash
python3 scripts/analysis/collate_binning_results.py \
    --run Soil_1 \
    --include-fdr \
    --strict
```

If the user wants to validate only the numbered-run output and intentionally
exclude matched-formula collation:

```bash
python3 scripts/analysis/collate_binning_results.py \
    --run Soil_1 \
    --include-fdr \
    --skip-matched \
    --strict
```

Important outputs include:

```text
results/<run>/binning_results_summary.tsv
results/<run>/main_partition_scan_recovery_summary.tsv
results/<run>/group_size_scan_recovery_summary.tsv
results/<run>/main_partition_peptide_protein_gains.tsv
results/<run>/group_size_peptide_protein_gains.tsv
```

Optional FDR summaries are produced with `--include-fdr`.

When `Matched_Formula_Output` exists, the collator can also write:

```text
results/matched_formula_peptide_protein_gains.tsv
results/matched_formula_input_manifest.tsv
results/matched_formula_results_guide.tsv
```

For field definitions and output meanings, consult:

```text
docs/Results Guide.tsv
```

---

## 6. How to help a novice interpret success

A successful run does not mean merely "the script exited."

Help the user verify:

- the expected output directories exist;
- the expected summary TSVs were written;
- there are no non-empty conflict-audit tables where conflicts are expected to
  be absent;
- corrected and baseline peptide/protein counts are plausible;
- warnings about missing analyses are understood rather than ignored.

For collation, the most useful compact output is usually:

```text
binning_results_summary.tsv
```

Key concepts:

- `RecoveredScans` describes recovery of full-database reference assignments
  across partitions;
- `PeptideGain` and `ProteinGain` are changes relative to the full-database
  baseline;
- positive gain means more entities recovered after correction;
- negative gain means fewer;
- FDR-change statistics use the sign convention documented below.

When comparing a rerun to a known-good result, prioritize:

```text
binning_results_summary.tsv
main_partition_peptide_protein_gains.tsv
group_size_peptide_protein_gains.tsv
```

Tiny last-digit floating-point formatting differences in FDR summaries can be
benign. Changes in peptide/protein identities or counts require investigation.

---

## 7. Troubleshooting approach for an AI helper

When a command fails:

1. identify the exact failing stage and error message;
2. confirm the referenced input path actually exists;
3. check whether the expected input file type is present;
4. distinguish a missing-input problem from a scientific/data conflict;
5. use `--dry-run` where available to inspect discovery logic;
6. avoid rewriting code until the failure is understood.

Useful checks include:

```bash
pwd
```

```bash
find data -maxdepth 3 -type d | sort
```

```bash
find data/Soil_Full/fdr_esti -maxdepth 1 -type f -name '*_withsyn.tsv' -print
```

```bash
bash -n scripts/pipeline/run_binning_pipeline.sh
```

```bash
python3 -m py_compile scripts/analysis/collate_binning_results.py
```

If there are multiple plausible reference files or dataset roots, do not guess.
Explain the ambiguity and help the user identify the correct one.

---

## 8. Scientific invariants

These are easy to accidentally "clean up" but are intentional parts of the
validated Binning workflow.

### Full-database results remain the final search space

Partition searches are used to estimate systematic FDR changes. Correction is
applied back to the original full-database `*_withsyn.tsv` results.

### Correction sign convention

```text
DeltaFDR = BinnedFDR - BaseFDR
```

Improvement is therefore negative.

Correction outputs are clipped to `[0, 1]`.

### Matched application

For multi-sample conditions, each numbered-run correction formula can be
applied only to the source sample from which that formula was derived.
`run_all_matched_formula_outputs.sh` orchestrates this behavior.

### Entity counting

- peptide count = distinct non-empty `Peptide`;
- `PeptideFlanked` is a parallel QC count;
- protein count = distinct non-empty `Protein`;
- crosstab source-row counts are audit/QC only.

Do **not** revert peptide counting to crosstab row counts.

### Classical-guided culling

Match on:

```text
(ScanNum, flank-stripped peptide)
```

while retaining PTM notation.

SIC rows take precedence over SYN fallback rows. SYN fallback contributes at
most one representative row per scan per partition and must contain a numeric
`SpecEValue` before FDR estimation.

---

## 9. Collator schema invariants

`scripts/analysis/collate_binning_results.py` is the source of truth for
collated outputs.

`docs/Results Guide.tsv` documents the current schema.

Important intentional labels:

```text
AnalysisFamily = Pairwise
```

for the primary Individual-vs-Pairwise comparison, and:

```text
AnalysisFamily = size_series
```

for the Individual-bin-size series.

Do **not** change `Pairwise` to `Main Partition`; that terminology was
deliberately rejected.

The current production peptide/protein gain tables must retain the distinct
entity-counting rules described above.

### Size-series labels

Size-series datasets are discovered dynamically from directories named:

```text
<Condition>_Groups_<LABEL>
```

`<LABEL>` is metadata, not a fixed parsing requirement. Familiar manuscript
labels such as `1.2M`, `600k`, `300k`, `150k`, and `75k` receive preferred
display ordering when present, but arbitrary labels are supported.

Do not hard-code the manuscript labels as the only valid size-series datasets.
`BinSize` carries the discovered label through collated outputs.

---

## 10. Root and path behavior

- `run_binning_pipeline.sh` intentionally requires explicit
  `--dataset-dir` and `--reference-withsyn` inputs.
- `collate_binning_results.py` defaults to `<Binning>/data`.
- `run_all_matched_formula_outputs.sh` defaults to `<Binning>/data`.
- `run_matched_correction_downstream.sh` defaults its data root to
  `<Binning>/data` but still requires a condition/design.
- Explicit root/path overrides are supported for external/HPC storage.

Keep orchestration-aware discovery in launchers. Low-level scientific utilities
should generally receive explicit paths rather than guess the data tree.

---

## 11. Release structure

The production analysis folder intentionally contains only generic workflow
utilities.

Dataset-specific manuscript audits, including historical Ocean-outlier or
Fecal/Ocean/Soil comparison scripts, are development/manuscript artifacts and
should not be restored to the public Binning arm without deliberate
generalization.

There is no active `config/` system in this arm. Do not tell a user to edit a
configuration file unless a future version actually implements one.

The `examples/` directory documents example layout and invocation, not a bundled
large test dataset.

---

## 12. Maintainer regression expectations

Before changing scientific behavior or output schemas:

1. run Python syntax checks and shell `bash -n`;
2. use `--dry-run` where supported;
3. regenerate a known run such as `Soil_1`;
4. compare:
   - `binning_results_summary.tsv`;
   - `main_partition_peptide_protein_gains.tsv`;
   - `group_size_peptide_protein_gains.tsv`;
5. investigate any scientific difference rather than normalizing it away.

Tiny last-digit floating-point formatting differences in FDR summaries can
occur without scientific differences; count or identity changes require
investigation.

---

## 13. Documentation responsibilities

If an output filename, column, counting rule, or meaning changes, update:

- `README.md` when the user workflow changes;
- `AI_README.md` when AI-helper guidance or scientific invariants change; and
- `docs/Results Guide.tsv` when collated outputs change.

When helping a user run the current release, prefer the existing public
interfaces and documentation over modifying scripts.
