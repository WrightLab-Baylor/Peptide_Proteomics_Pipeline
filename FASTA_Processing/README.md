# FASTA Processing

FASTA Processing provides the shared upstream database-preparation utilities used
by the related **Binning** and **Targeted Peptide Removal** workflows.

It has two main jobs:

1. cluster and digest a protein FASTA into reusable peptide-level resources; and
2. partition a peptide FASTA into Binning-compatible Individual and Pairwise
   databases.

This repository intentionally stops at database preparation. Search, FDR,
downstream peptide/protein analysis, and result collation belong to the related
workflow repositories.

> Repository links for FASTA Processing, Binning, and Targeted Peptide Removal
> will be added after the final GitHub repositories are created.

## Repository layout

```text
FASTA_Processing/
├── README.md
├── AI_README.md
├── data/
└── scripts/
    ├── fasta/
    │   ├── fasta_cluster_generator.py
    │   └── peptide_fasta_builder.py
    ├── binning/
    │   └── peptide_bin_pair_builder.py
    └── pipeline/
        ├── prepare_peptide_fasta.sh
        └── prepare_binning_databases.sh
```

`data/` is the default local output root for FASTA processing.

A shared Conda environment for FASTA Processing, Binning, and Targeted Peptide
Removal is maintained outside the three individual repositories.

---

# 1. Prepare a peptide FASTA

The main user-facing launcher is:

```text
scripts/pipeline/prepare_peptide_fasta.sh
```

It performs:

```text
protein FASTA
    ↓
CD-HIT protein clustering
    ↓
cluster-annotated protein FASTA
    ↓
tryptic digestion
    ↓
peptide FASTA + peptide/protein digest map
```

The default settings are:

```text
protein clustering identity     0.95
CD-HIT word length              5
CD-HIT threads                  all available CPUs
minimum peptide length          6
maximum peptide length          50
target/decoy overlap removal    disabled
```

Digestion uses tryptic cleavage after K/R except when followed by P.

## Basic example

From the FASTA Processing repository root:

```bash
bash scripts/pipeline/prepare_peptide_fasta.sh \
    --condition Soil \
    --protein-fasta /path/to/Soil_proteins.fasta
```

By default, outputs are written to:

```text
data/Soil/
├── protein/
│   └── Soil_clustered.fasta
└── peptide/
    ├── Soil_clustered_peptide.fasta
    ├── peptide_protein_digest_map.tsv
    ├── digest_stats.txt
    └── removed_peptides.tsv
```

These outputs are intended to be reused by the downstream Binning and Targeted
Peptide Removal workflows.

## Useful overrides

Choose a different output root:

```bash
bash scripts/pipeline/prepare_peptide_fasta.sh \
    --condition Soil \
    --protein-fasta /path/to/Soil_proteins.fasta \
    --output-root /path/to/output
```

Preview the commands without running them:

```bash
bash scripts/pipeline/prepare_peptide_fasta.sh \
    --condition Soil \
    --protein-fasta /path/to/Soil_proteins.fasta \
    --dry-run
```

Use a specific CD-HIT executable:

```bash
bash scripts/pipeline/prepare_peptide_fasta.sh \
    --condition Soil \
    --protein-fasta /path/to/Soil_proteins.fasta \
    --cdhit-path /path/to/cd-hit
```

The cluster generator first checks whether the protein FASTA already contains
usable `Cluster=` annotations. If it does, CD-HIT is not required.

## Target/decoy-overlap removal

`peptide_fasta_builder.py` retains an optional target/decoy-overlap removal
mode:

```text
--remove-tda-overlaps
```

This is **disabled by default**.

The validated current workflow does not automatically remove peptides merely
because they are palindromic or share target/decoy sequence relationships.
Those relationships can instead be evaluated downstream by the Targeted Peptide
Removal workflow.

---

# 2. Prepare databases for Binning

The Binning database launcher is:

```text
scripts/pipeline/prepare_binning_databases.sh
```

It takes an existing peptide FASTA and writes datasets directly into a
Binning-compatible data tree.

The default destination is a sibling:

```text
../Binning/data
```

relative to the FASTA Processing repository. If the Binning data root is
elsewhere, provide it explicitly with:

```text
--binning-data-root /path/to/Binning/data
```

The launcher supports two related database designs:

1. a primary Individual/Pairwise partition; and
2. an Individual-bin size series.

## 2A. Primary Individual/Pairwise databases

The production default is:

```text
9 Individual bins
36 Pairwise bins
```

Example:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta
```

This writes under:

```text
Binning/data/Soil_1/Soil_Binning/
```

with Individual-bin folders such as:

```text
..._group_1
..._group_2
...
..._group_9
```

and Pairwise-bin folders such as:

```text
..._pair_1_2
..._pair_1_3
...
```

Each generated partition contains the `database/` and `data/` subdirectories
expected by the Binning workflow.

The condition name is normally inferred from the run name by removing a trailing
number, so `Soil_1` becomes `Soil`.

### Custom primary partition

Change the number of Individual bins:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Test_1 \
    --peptide-fasta /path/to/Test_peptides.fasta \
    --num-groups 12
```

Generate only Individual bins:

```text
--output-mode groups
```

Generate only Pairwise bins:

```text
--output-mode pairs
```

Generate both:

```text
--output-mode both
```

## 2B. Individual-bin size series

The default size-series contains five familiar manuscript labels:

```text
1.2M,600k,300k,150k,75k
```

but neither the **number of series entries** nor their labels is fixed. Binning discovers size-series
datasets generically from:

```text
<Condition>_Groups_<LABEL>
```

so custom labels such as `500k`, `PilotA`, or `Small` are valid.

There are **two independent choices** when defining a size series:

1. how the number of Individual bins is determined; and
2. what label is used for each output dataset.

### Method 1: target peptides per Individual bin

Let the launcher calculate the nearest whole number of groups from the peptide
FASTA size:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta \
    --build size-series
```

The default targets are:

```text
1.2M,600k,300k,150k,75k
```

Custom targets can be supplied with:

```text
--size-targets 2M,1M,500k
```

By default, those targets also become the dataset labels. Labels can instead be
set independently:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta \
    --build size-series \
    --size-targets 1M,500k,250k \
    --size-labels Large,Medium,Small
```

### Method 2: exact group counts

If the desired number of groups is already known, supply those counts directly:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta \
    --build size-series \
    --size-groups 11,22,44,89,177
```

If `--size-labels` is omitted, exact group counts are paired positionally with
`--size-targets`, preserving the familiar default labels:

```text
size target/label    exact groups
1.2M                 11
600k                 22
300k                 44
150k                 89
75k                  177
```

For a completely custom exact-group series, provide labels directly:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta \
    --build size-series \
    --size-groups 10,20,40 \
    --size-labels Small,Medium,Large
```

Exact group counts can also stand alone:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta \
    --build size-series \
    --size-groups 10,20,40
```

Without explicit targets or labels, those datasets are named:

```text
Soil_Groups_G10
Soil_Groups_G20
Soil_Groups_G40
```

This means a size series can contain any number of entries; five is only the
default manuscript series.

This writes:

```text
Soil_Groups_Small
Soil_Groups_Medium
Soil_Groups_Large
```

The Binning collator treats `<LABEL>` as descriptive metadata rather than a
required numeric size. Exact-group mode is useful when reproducing a predefined
partitioning scheme or when a custom series is easier to describe by group
count than by nominal peptides per bin.

## 2C. Build the primary and size-series databases together

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta \
    --build all
```

With exact size-series group counts:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta \
    --build all \
    --size-groups 11,22,44,89,177
```

---

# 3. Reproducibility and safety

## Shuffle seed

The bin builder performs a master shuffle before partitioning.

For reproducible partitions, provide:

```text
--seed 12345
```

The same seed is passed to every requested build.

## Dry run

Before writing Binning datasets:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta /path/to/Soil_peptides.fasta \
    --build all \
    --dry-run
```

The launcher refuses to overwrite a non-empty existing Binning dataset
directory.

## Shared partition logic

Individual and Pairwise databases are built from the same master shuffle and
Individual-bin boundaries. Pairwise bins therefore represent combinations of
the same underlying Individual bins rather than independent reshuffles.

---

# 4. How this repository connects to the other workflows

## Binning

Binning consumes peptide FASTA partitions created by
`prepare_binning_databases.sh`.

Search, classical-guided culling, partition-level FDR analysis, correctional
formula processing, downstream peptide/protein analysis, and result collation
are handled in the Binning repository.

## Targeted Peptide Removal

Targeted Peptide Removal consumes the clustered peptide FASTA,
`peptide_protein_digest_map.tsv`, and corresponding protein FASTA generated by
the FASTA-processing workflow.

Its evaluator, targeted conflict removal, downstream analysis, and result
collation are handled in the Targeted Peptide Removal repository.

---

# 5. Public scripts

## `scripts/fasta/fasta_cluster_generator.py`

Clusters proteins with CD-HIT and annotates retained FASTA headers with cluster
identifiers.

## `scripts/fasta/peptide_fasta_builder.py`

Digests the clustered protein FASTA into a peptide FASTA and writes the
peptide-to-protein digest map and digest statistics.

## `scripts/binning/peptide_bin_pair_builder.py`

Builds Individual and/or Pairwise peptide FASTA partitions from one master
shuffle.

## `scripts/pipeline/prepare_peptide_fasta.sh`

Recommended user-facing launcher for protein clustering and peptide digestion.

## `scripts/pipeline/prepare_binning_databases.sh`

Recommended user-facing launcher for primary and size-series Binning databases.

For AI-assisted operation or maintenance guidance, see `AI_README.md`.
