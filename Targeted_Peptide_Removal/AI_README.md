# AI Helper / Maintainer Guide — Targeted Peptide Removal

Repository: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline
Workflow directory: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline/tree/main/Targeted_Peptide_Removal

Related workflow directories:

- FASTA Processing: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline/tree/main/FASTA_Processing
- Binning: https://github.com/WrightLab-Baylor/Peptide_Proteomics_Pipeline/tree/main/Binning

This file is an anchor point for an AI assistant helping a user **run,
understand, troubleshoot, or maintain** the Targeted Peptide Removal workflow.

The intended user may be a new graduate student with little or no command-line
experience. Guide the user through the validated public interfaces first.
Do not begin by rewriting scripts.

---

## 1. Understand the three-workflow repository

This repository does **not** begin from a protein FASTA.

Upstream FASTA processing is handled by the sibling `FASTA_Processing/` workflow
directory, which digests the protein FASTA into a peptide FASTA and produces the
peptide/protein digest map.
The Targeted Peptide Removal workflow begins with the resulting peptide FASTA.

The sibling `Binning/` workflow directory contains the Binning workflow.

When cross-workflow resources are needed, inspect all three workflow directories
in this repository before inventing missing functionality.

Conceptually:

```text
FASTA_Processing workflow directory
    |
    | protein FASTA -> peptide FASTA + peptide/protein digest map
    v
Targeted Peptide Removal
    |
    | evaluator -> culler -> threshold series -> downstream -> collator
    v
peptide/protein gain-loss summaries
```

Do not reimplement upstream digestion logic inside this repository merely to
make a local run convenient.

---

## 2. First priority: orient the user

When a user asks how to run this workflow, determine:

1. the Targeted Peptide Removal repository root;
2. the peptide FASTA produced by the upstream FASTA-processing workflow;
3. the condition name, such as `Fecal`, `Ocean`, or `Soil`;
4. where evaluator outputs should be stored;
5. whether the unmodified baseline run already exists;
6. whether threshold runs already contain search/downstream outputs;
7. whether the protein FASTA and `peptide_protein_digest_map.tsv` have been
   staged into the run directories.

The repository-local data root is:

```text
Targeted_Peptide_Removal/data/
```

For a novice user:

- give one command at a time;
- explain what it will do before asking them to run it;
- prefer `--dry-run` where available;
- use absolute paths when directory context is unclear;
- do not guess among multiple FASTA, digest-map, or run candidates;
- do not modify scripts just to accommodate a local filesystem path.

Useful orientation commands include:

```bash
pwd
```

```bash
find data -maxdepth 3 -type d | sort
```

and, when locating resources:

```bash
find /path/to/project -maxdepth 3 -type f \
    \( -name '*.fasta' -o -name '*.faa' -o -name 'peptide_protein_digest_map.tsv' \) \
    -print
```

---

## 3. Production workflow: evaluator

Use:

```text
scripts/removal/peptide_fasta_evaluator.py
```

For the targeted-removal workflow, the main production subcommand is:

```text
conflicts
```

Validated defaults:

```text
k = 3
min_shared_kmers = 3
jaccard_threshold = 0.5
lcs_threshold = 0.5
max_candidates = 500
```

Similarity thresholds below `0.5` are intentionally rejected. Do not advise a
user to bypass this floor casually: lower thresholds cause an impractical
candidate/memory explosion for large metaproteomic peptide databases.

Typical command:

```bash
python3 scripts/removal/peptide_fasta_evaluator.py conflicts \
    --target_pep_fasta /path/to/Condition_peptides.fasta \
    --out_dir /path/to/evaluator_output \
    --workers 32
```

Do not invent the worker count. Help the user choose one appropriate to the
machine.

The important output for the culling workflow is:

```text
conflicts.tsv
```

The evaluator also has `qc` and `prune` modes. The validated threshold-series
workflow uses `conflicts` followed by `conflict_culler_fast.py`.

Palindrome and target/decoy-overlap statistics are QC information. The current
workflow does not automatically remove peptides merely because they are
palindromic or have target/decoy overlap.

---

## 4. Production workflow: culling and threshold series

The conservative culler is:

```text
scripts/removal/conflict_culler_fast.py
```

It removes both endpoints of every conflict edge meeting the requested
Jaccard/LCS thresholds.

The preferred user-facing orchestration is:

```text
scripts/pipeline/generate_peptide_removal_series.sh
```

Default series:

```text
1.0,0.9,0.8,0.7,0.6,0.5
```

Typical command:

```bash
bash scripts/pipeline/generate_peptide_removal_series.sh \
    --condition Soil \
    --peptide-fasta /path/to/Soil_peptides.fasta \
    --conflicts-tsv /path/to/evaluator_output/conflicts.tsv
```

Default output root:

```text
<Targeted_Peptide_Removal>/data
```

Expected threshold layout:

```text
data/Soil_0.9/
├── data/
├── database/
│   └── Soil_0.9.fasta
└── culling/
```

The threshold generator intentionally does **not** copy the shared protein FASTA
or `peptide_protein_digest_map.tsv`.

Use the validated stager:

```text
scripts/pipeline/stage_shared_fasta_resources.sh
```

With the standard sibling-repository layout, a typical command is:

```bash
bash scripts/pipeline/stage_shared_fasta_resources.sh \
    --condition Soil \
    --dry-run
```

After reviewing the dry run, repeat without `--dry-run`.

The stager discovers the clustered protein FASTA and digest map from the
corresponding `FASTA_Processing/data/<Condition>/` directory and creates
relative symbolic links in existing baseline and threshold run directories.
It never overwrites regular files. Do not replace this with ad hoc copies unless
there is a specific reason to duplicate the resources.

Do not overwrite a threshold directory that already contains search/downstream
results merely to regenerate its peptide FASTA.

---

## 5. Baseline requirement

The collator requires an unmodified baseline condition alongside the threshold
runs.

Example:

```text
data/
├── Soil/
├── Soil_1.0/
├── Soil_0.9/
├── Soil_0.8/
├── Soil_0.7/
├── Soil_0.6/
└── Soil_0.5/
```

`generate_peptide_removal_series.sh` does not create `data/Soil/`.

Make sure the user understands that all gain/loss values are relative to this
unmodified baseline.

---

## 6. Shared-resource staging

Before search/downstream processing, use:

```text
scripts/pipeline/stage_shared_fasta_resources.sh
```

The helper is deliberately conservative:

- it operates only on existing `<Condition>` / `<Condition>_*` run directories;
- it creates relative symbolic links;
- it leaves regular files untouched;
- it can replace wrong/stale symlinks only when `--replace-links` is supplied;
- it supports explicit source/destination roots for non-sibling repository layouts.

Run it with `--dry-run` first when guiding a novice user.

---

## 7. Search/downstream handoff

This repository does not perform the upstream RAW/search/MASIC workflow that
creates `SICdir`.

Before using the downstream launcher, each baseline/threshold run should contain
the resources expected by:

```text
scripts/pipeline/run_peptide_removal_downstream.sh
```

Required run layout includes:

```text
<run>/
├── SICdir/
├── database/
├── peptide_protein_digest_map.tsv
└── <one top-level protein FASTA or FAA>
```

If these are absent, identify the missing upstream step instead of modifying
the downstream launcher to guess resources.

---

## 8. Production workflow: downstream

Use:

```text
scripts/pipeline/run_peptide_removal_downstream.sh
```

Always consider a dry run first:

```bash
bash scripts/pipeline/run_peptide_removal_downstream.sh \
    --run-dir data/Soil_0.5 \
    --dry-run
```

Validated defaults:

```text
peptide filter             MSMSScore
stringency                 8
second filter              FDR
second stringency          0.05
highest intensity only     True
unique to cluster          False
uniqueness mode            all_matches
minimum peptide count      1
peptide counting mode      per_sample
cluster sums               True
rollup type                sum
RRollup mode               None
```

The launcher resolves downstream scripts from this repository. Do not tell the
user to copy a `msgpypulse/` directory into every run.

### CD-HIT

The cluster generator checks for existing cluster information **before**
requiring CD-HIT.

Therefore, do not diagnose a missing `cd-hit` executable unless clustering is
actually needed.

If clustering is needed, CD-HIT may be:

- available as `cd-hit` on `PATH`; or
- supplied explicitly using `--cdhit-path`.

A future task for the FASTA_Processing workflow directory is to include
AI-facing installation guidance for CD-HIT.

---

## 9. Production workflow: collation

Use:

```text
scripts/analysis/collate_peptide_removal_results.py
```

With the standard layout:

```bash
python3 scripts/analysis/collate_peptide_removal_results.py
```

Default root:

```text
<Targeted_Peptide_Removal>/data
```

Default output:

```text
<Targeted_Peptide_Removal>/data/results
```

The collator discovers a complete baseline and its `Condition_<threshold>`
siblings.

### Critical counting invariants

The collator has been explicitly audited for distinct-entity counting.

- `Database_Peptides` = number of peptide FASTA records.
- identified peptides = distinct non-empty values in the `Peptide` column.
- identified proteins = distinct non-empty values in the `Protein` column.
- gained/lost peptide and protein values are set differences of those distinct
  entities.

**Do not replace these calculations with row counts.**

A peptide expanded across multiple protein mappings still counts as one peptide.

---

## 10. Scientific invariants

### Threshold floor

Similarity thresholds below `0.5` are intentionally unsupported.

### Conflict definition

The evaluator assesses reverse-target/target sequence similarity using the
configured k-mer Jaccard and LCS criteria.

Reciprocal filtering is not part of the current public workflow.

### Conservative removal

`conflict_culler_fast.py` removes both endpoints of every qualifying conflict
edge. Do not silently change this to a greedy/minimal-removal strategy without
scientific validation.

### No automatic palindrome removal

Do not automatically eliminate palindromic peptides or target/decoy-overlap
peptides simply because they are detected.

### Baseline comparison

Peptide/protein gains and losses are always calculated relative to the
unmodified baseline condition.

---

## 11. Troubleshooting approach

When a command fails:

1. identify the exact script and stage;
2. read the actual error;
3. confirm the referenced path exists;
4. inspect the minimum relevant directory level;
5. distinguish missing upstream resources from a code failure;
6. use `--dry-run` where supported;
7. avoid modifying scientific code until the failure is understood.

Useful checks:

```bash
bash -n scripts/pipeline/generate_peptide_removal_series.sh
```

```bash
bash -n scripts/pipeline/run_peptide_removal_downstream.sh
```

```bash
python3 -m py_compile scripts/removal/peptide_fasta_evaluator.py
```

```bash
python3 -m py_compile scripts/analysis/collate_peptide_removal_results.py
```

If a user supplies an ambiguous directory tree, help them locate the correct
files rather than guessing.

---

## 12. Release/maintenance notes

The downstream Python utilities are intentionally carried with this repository
so the peptide-removal branch can run independently.

A future repository cleanup may centralize shared code, but do not relocate or
deduplicate these utilities unless every launcher and dependency has been
updated and regression-tested.

There is no active per-repository `environment.yml` intended for the final
three-workflow release. A shared environment definition will be maintained
outside the three workflow repositories.

---

## 13. Regression expectations

Before changing evaluator, culler, downstream, or collator scientific behavior:

1. run Python compilation and shell `bash -n`;
2. exercise `--help` for public interfaces;
3. test evaluator → culler → threshold generation on a small fixture;
4. confirm thresholds below `0.5` are rejected;
5. dry-run downstream processing;
6. rerun a known real condition when possible;
7. compare distinct peptide/protein counts and gain/loss results against a
   known-good output.

Do not normalize away count or identity differences.

---

## 14. Documentation responsibilities

If user-facing behavior changes, update:

- `README.md`;
- `AI_README.md`;
- any future results/output guide for this repository.

When repository URLs are available, add them near the top of this file and
explicitly instruct the helper to inspect all three related repositories before
inventing missing upstream functionality.
