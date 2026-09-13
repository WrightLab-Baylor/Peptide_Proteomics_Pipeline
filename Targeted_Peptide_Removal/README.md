# Targeted Peptide Removal

Targeted Peptide Removal evaluates peptide FASTA databases for high target–decoy
sequence overlap, generates progressively culled peptide databases, runs the
standard peptide/protein downstream workflow on completed searches, and
summarizes peptide/protein gains and losses relative to an unmodified baseline.

This repository contains the **peptide-removal-specific** portion of the
workflow. Protein FASTA digestion and peptide FASTA generation are handled by a
separate FASTA-processing repository.

## Workflow at a glance

```text
protein FASTA
    |
    |  upstream FASTA-processing repository
    v
target peptide FASTA + peptide/protein digest map
    |
    v
peptide_fasta_evaluator.py
    |
    |  conflicts.tsv
    v
conflict_culler_fast.py
    |
    |  repeated by generate_peptide_removal_series.sh
    v
Condition_1.0 ... Condition_0.5 peptide FASTAs
    |
    |  search / SIC generation performed by the normal proteomics workflow
    v
run_peptide_removal_downstream.sh
    |
    v
peptide/protein crosstabs
    |
    v
collate_peptide_removal_results.py
    |
    v
condition-level peptide/protein gain/loss summaries
```

## Repository layout

```text
Targeted_Peptide_Removal/
├── README.md
├── AI_README.md
├── data/
└── scripts/
    ├── analysis/
    │   └── collate_peptide_removal_results.py
    ├── removal/
    │   ├── peptide_fasta_evaluator.py
    │   └── conflict_culler_fast.py
    ├── pipeline/
    │   ├── generate_peptide_removal_series.sh
    │   ├── stage_shared_fasta_resources.sh
    │   └── run_peptide_removal_downstream.sh
    └── downstream/
        └── downstream Python utilities
```

`data/` is the default local project root. Large datasets may live elsewhere;
all public launchers accept explicit paths where appropriate.

## Requirements

The workflow requires Python 3 and the Python packages used by the scripts,
including `pandas` and `numpy`.

CD-HIT (`cd-hit`) is required only when downstream processing needs to generate
protein clusters and no suitable cluster annotations already exist. If CD-HIT
is not on `PATH`, an explicit executable can be supplied to the downstream
launcher with `--cdhit-path`.

The shared environment for all three related workflow repositories will be
maintained outside this repository.

---

# 1. Evaluate the peptide FASTA

The primary evaluator is:

```text
scripts/removal/peptide_fasta_evaluator.py
```

For the production targeted-removal workflow, use the `conflicts` subcommand.
The validated defaults are:

```text
k-mer size               3
minimum shared k-mers    3
Jaccard threshold        0.5
LCS threshold            0.5
maximum candidates       500
```

Similarity thresholds below `0.5` are intentionally rejected because the
candidate space and memory requirements become impractical for large
metaproteomic peptide databases.

Example:

```bash
python3 scripts/removal/peptide_fasta_evaluator.py conflicts \
    --target_pep_fasta /path/to/Condition_peptides.fasta \
    --out_dir /path/to/evaluator_output \
    --workers 32
```

The main product used by the removal workflow is:

```text
evaluator_output/conflicts.tsv
```

The evaluator also provides `qc` and `prune` subcommands. The production
threshold-series workflow uses `conflicts` followed by
`conflict_culler_fast.py`; the integrated `prune` mode is retained as an
optional utility.

The evaluator reports target/decoy overlap and palindrome-related QC. It does
**not** automatically eliminate palindromic peptides or target–decoy overlaps
simply because they exist.

---

# 2. Generate the peptide-removal threshold series

The threshold-series launcher is:

```text
scripts/pipeline/generate_peptide_removal_series.sh
```

By default it generates:

```text
1.0, 0.9, 0.8, 0.7, 0.6, 0.5
```

where each value is the minimum Jaccard and LCS threshold passed to the
conservative conflict culler.

Example from the repository root:

```bash
bash scripts/pipeline/generate_peptide_removal_series.sh \
    --condition Soil \
    --peptide-fasta /path/to/Soil_peptides.fasta \
    --conflicts-tsv /path/to/evaluator_output/conflicts.tsv
```

The default output root is:

```text
Targeted_Peptide_Removal/data/
```

For example:

```text
data/
├── Soil_1.0/
│   ├── data/
│   ├── database/
│   │   └── Soil_1.0.fasta
│   └── culling/
├── Soil_0.9/
│   ├── data/
│   ├── database/
│   │   └── Soil_0.9.fasta
│   └── culling/
├── ...
└── Soil_0.5/
    ├── data/
    ├── database/
    │   └── Soil_0.5.fasta
    └── culling/
```

Each `culling/` directory contains the culler audit products, including:

```text
cull_report.tsv
removed_peptides.tsv
remaining_conflicts_after_cull.tsv
```

The culler uses a conservative rule: if a peptide participates in a qualifying
conflict edge, both endpoints are removed.

The threshold generator intentionally does not duplicate the shared protein FASTA
or `peptide_protein_digest_map.tsv`. After the baseline and threshold run
folders exist, stage those resources from the sibling FASTA Processing output
with:

```bash
bash scripts/pipeline/stage_shared_fasta_resources.sh \
    --condition Soil
```

With the standard sibling repository layout, the stager discovers:

```text
FASTA_Processing/data/Soil/protein/*_clustered.fasta
FASTA_Processing/data/Soil/peptide/peptide_protein_digest_map.tsv
```

and creates relative symbolic links in every existing `Soil` and `Soil_*` run
directory under `Targeted_Peptide_Removal/data/`.

Preview the links first with:

```text
--dry-run
```

If the repositories or data live elsewhere, use `--fasta-data-root`,
`--protein-fasta`, `--digest-map`, or `--data-root` explicitly. Regular files
are never overwritten.

---

# 3. Prepare the baseline and stage shared resources

The collator compares every threshold run against an **unmodified baseline**
directory with the condition name alone.

For example:

```text
data/
├── Soil/        # unmodified baseline
├── Soil_1.0/
├── Soil_0.9/
├── Soil_0.8/
├── Soil_0.7/
├── Soil_0.6/
└── Soil_0.5/
```

The threshold-series generator creates the threshold directories, but it does
not create the baseline directory. Create or stage the unmodified baseline run,
then run `stage_shared_fasta_resources.sh` so the baseline and threshold runs
share the same authoritative protein FASTA and digest map.

Before running the downstream launcher, each baseline/threshold run must have
the outputs and resources expected by the downstream workflow. In particular,
the launcher expects:

```text
<run>/
├── SICdir/
├── database/
├── peptide_protein_digest_map.tsv
└── <one top-level protein FASTA or FAA>
```

The RAW-data search and SIC-generation stages occur outside this repository.
Use the peptide FASTA in each run's `database/` directory as the search
database, then return here once the required SIC/search outputs are present.

---

# 4. Run searches and peptide/protein downstream processing

The public downstream launcher is:

```text
scripts/pipeline/run_peptide_removal_downstream.sh
```

Example:

```bash
bash scripts/pipeline/run_peptide_removal_downstream.sh \
    --run-dir data/Soil_0.5
```

A dry run is recommended first:

```bash
bash scripts/pipeline/run_peptide_removal_downstream.sh \
    --run-dir data/Soil_0.5 \
    --dry-run
```

Validated defaults are:

```text
peptide filter             MSMSScore
MSMSScore threshold        8
second filter              FDR
FDR threshold              0.05
highest intensity only     True
unique to cluster          False
uniqueness mode            all_matches
minimum peptide count      1
peptide counting mode      per_sample
cluster sums               True
rollup type                sum
RRollup mode               None
```

If CD-HIT is required and not on `PATH`:

```bash
bash scripts/pipeline/run_peptide_removal_downstream.sh \
    --run-dir data/Soil_0.5 \
    --cdhit-path /path/to/cd-hit
```

The cluster generator first checks whether suitable cluster information already
exists. CD-HIT is required only when clustering actually needs to be generated.

Run the same downstream workflow for the unmodified baseline.

---

# 5. Collate peptide-removal results

The collator is:

```text
scripts/analysis/collate_peptide_removal_results.py
```

With the standard repository layout:

```bash
python3 scripts/analysis/collate_peptide_removal_results.py
```

The default input root is:

```text
Targeted_Peptide_Removal/data/
```

and the default output directory is:

```text
Targeted_Peptide_Removal/data/results/
```

The collator discovers complete baseline conditions and their
`Condition_<threshold>` siblings automatically.

It measures:

- peptide database size;
- database peptides removed;
- fraction of the peptide database removed;
- identified peptides passing downstream filters;
- peptides gained and lost relative to baseline;
- net peptide change;
- proteins passing;
- proteins gained and lost relative to baseline;
- net protein change.

## Counting rules

These definitions are important:

- `Database_Peptides` = FASTA records in the peptide search database.
- identified peptide count = **distinct non-empty `Peptide` values**.
- identified protein count = **distinct non-empty `Protein` values**.
- gained/lost peptide and protein counts are **set differences of distinct
  entities**, not TSV row counts.

A peptide mapped to multiple proteins therefore still counts as one identified
peptide.

---

# Safety and reproducibility notes

- Similarity thresholds below `0.5` are intentionally unsupported.
- Do not use crosstab row counts as peptide counts.
- Keep an unmodified baseline run for every condition.
- The threshold generator will not overwrite run directories that appear to
  contain search/downstream results.
- Prefer `--dry-run` before rerunning downstream processing.
- Preserve culling reports and collated summaries as provenance for each run.
- The protein FASTA and peptide/protein digest map should come from the same
  upstream FASTA-processing workflow used to produce the peptide search
  database.

For an AI assistant or future maintainer, see `AI_README.md`.
