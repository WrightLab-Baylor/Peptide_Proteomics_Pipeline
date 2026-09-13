# AI Helper / Maintainer Guide — FASTA Processing

This file is an anchor point for an AI assistant helping a user **run,
understand, troubleshoot, or maintain** the FASTA Processing workflow.

The intended user may be a new graduate student with little or no command-line
experience. Prefer the public launchers, explain commands before asking the user
to run them, and do not rewrite scientific code merely to accommodate a local
path.

---

## 1. Understand this repository's role

FASTA Processing is the shared upstream repository for two related workflows:

```text
FASTA Processing
    ├── Binning
    └── Targeted Peptide Removal
```

Its responsibilities are limited to:

1. protein FASTA clustering;
2. protein-to-peptide digestion and digest-map generation;
3. generation of Individual/Pairwise peptide databases for Binning.

It does **not** perform MS-GF+ searching, FDR estimation, peptide/protein
downstream processing, correctional-formula analysis, or result collation.

Those operations belong to the Binning and Targeted Peptide Removal
repositories.

Once GitHub URLs are finalized, add links to all three repositories near the
top of this file. Future AI helpers should be instructed to inspect all three
before inventing missing functionality.

---

## 2. First priority: orient the user

When a user asks how to prepare a FASTA, identify:

1. the FASTA Processing repository root;
2. the input protein FASTA;
3. the condition/dataset name;
4. whether the input FASTA is already cluster-annotated;
5. whether the user needs:
   - only the reusable peptide FASTA resources;
   - Binning databases;
   - or both;
6. where the related Binning repository/data root lives if Binning databases are
   requested.

For a novice:

- give one command at a time;
- explain expected outputs;
- use `--dry-run` when supported;
- use absolute paths if the working directory is unclear;
- do not guess between multiple protein or peptide FASTAs;
- do not move or rename user inputs unless explicitly requested.

Useful orientation commands:

```bash
pwd
```

```bash
find data -maxdepth 3 -type f | sort
```

and, if a related repository is elsewhere:

```bash
find /path/to/project -maxdepth 3 -type d | sort
```

---

## 3. Primary FASTA-processing workflow

Use:

```text
scripts/pipeline/prepare_peptide_fasta.sh
```

Typical command:

```bash
bash scripts/pipeline/prepare_peptide_fasta.sh \
    --condition Soil \
    --protein-fasta /path/to/Soil_proteins.fasta
```

Default output:

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

Validated defaults:

```text
protein cluster identity     0.95
CD-HIT word length           5
CD-HIT threads               0 / all available CPUs
minimum peptide length       6
maximum peptide length       50
TDA-overlap removal          disabled
```

Digestion is tryptic: cleavage after K/R except before P.

Do not silently change these defaults when helping a user reproduce the
validated workflow.

---

## 4. CD-HIT behavior and installation

`scripts/fasta/fasta_cluster_generator.py` checks whether the FASTA already
contains usable `Cluster=` annotations **before** requiring CD-HIT.

Therefore:

- if cluster annotations are already present, do not require CD-HIT;
- if clustering is needed, use `cd-hit` from `PATH` or an explicit
  `--cdhit-path`.

### Preferred installation with the shared Conda environment

The shared environment includes `cd-hit` from Bioconda.

If helping a user create or update the environment manually, a typical Conda
installation is:

```bash
conda install -c conda-forge -c bioconda cd-hit
```

or include:

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - cd-hit
```

in the shared environment file.

After installation, verify:

```bash
which cd-hit
```

and:

```bash
cd-hit -h
```

If CD-HIT is installed somewhere that is not on `PATH`, pass the executable:

```text
--cdhit-path /path/to/cd-hit
```

Do not add personal workstation paths to the source code.

---

## 5. Peptide FASTA scientific invariants

The public builder is:

```text
scripts/fasta/peptide_fasta_builder.py
```

Key behavior:

- tryptic cleavage after K/R except before P;
- default peptide length range 6–50;
- peptide-to-protein mappings are written to
  `peptide_protein_digest_map.tsv`;
- target/decoy-overlap removal exists as an option but is disabled by default.

Do not automatically remove palindromic peptides or target/decoy sequence
relationships merely because they are detected.

Targeted overlap evaluation/removal belongs to the Targeted Peptide Removal
repository.

---

## 6. Preparing the primary Binning databases

Use:

```text
scripts/pipeline/prepare_binning_databases.sh
```

The underlying builder is:

```text
scripts/binning/peptide_bin_pair_builder.py
```

but most users should use the launcher.

Production default:

```text
9 Individual bins
36 Pairwise bins
```

Typical command:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta data/Soil/peptide/Soil_clustered_peptide.fasta
```

Default Binning destination:

```text
../Binning/data
```

If the related Binning repository is elsewhere:

```text
--binning-data-root /path/to/Binning/data
```

Expected dataset:

```text
Binning/data/Soil_1/Soil_Binning/
```

The launcher intentionally refuses to overwrite a non-empty existing dataset.

---

## 7. Binning size series

The launcher supports both calculated and exact group-count approaches.

The default series has five familiar manuscript labels:

```text
1.2M,600k,300k,150k,75k
```

Five entries are **not required**. These are defaults and familiar manuscript
labels, not a required naming or series-length contract. Current Binning discovery accepts any:

```text
<Condition>_Groups_<LABEL>
```

and carries `<LABEL>` through as `BinSize`. Familiar paper labels receive
preferred display ordering when present, but custom labels are valid.

### Calculated group counts

If `--size-groups` is absent, the launcher counts FASTA records and converts
each target peptides/bin value into the nearest whole number of groups.

Example:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta /path/to/Soil_peptides.fasta \
    --build size-series
```

### Exact group counts

When reproducing a predefined experiment, exact group counts can be supplied:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta /path/to/Soil_peptides.fasta \
    --build size-series \
    --size-groups 11,22,44,89,177
```

If `--size-labels` is omitted, `--size-groups` is positional relative to
`--size-targets` and those target values supply the labels.

Custom labels can be supplied independently:

```bash
bash scripts/pipeline/prepare_binning_databases.sh \
    --run-name Soil_1 \
    --peptide-fasta /path/to/Soil_peptides.fasta \
    --build size-series \
    --size-groups 10,20,40 \
    --size-labels Small,Medium,Large
```

`--size-labels` may also be used with calculated `--size-targets`.

Exact group counts can be supplied without targets or labels. In that case,
automatic labels are generated as `G<count>` (for example `G10`, `G20`,
`G40`). This allows arbitrary-length group-count series without placeholder
size targets.

Do not describe exact group-count mode as approximate. The supplied counts are
used directly, and labels are descriptive metadata rather than group-count
inputs.

---

## 8. Binning partition invariants

The bin builder uses one master shuffle.

Individual bins are slices of that master shuffle.

Pairwise bins are combinations of the same Individual bins.

Therefore:

- do not reshuffle independently for each pair;
- do not generate Pairwise bins from independently generated Individual
  partitions;
- preserve the requested seed when reproducibility matters.

The public terminology is:

- **Individual bins**
- **Pairwise bins**

Machine-facing output names retain:

```text
*_group_<N>
*_pair_<N>_<M>
```

because the Binning workflow consumes them.

---

## 9. Reproducibility

For a reproducible master shuffle:

```text
--seed 12345
```

The same seed is applied across primary and size-series builds in one launcher
invocation.

Before creating databases, recommend:

```text
--dry-run
```

to inspect paths and group counts.

When the user supplies exact size-series group counts:

- if `--size-labels` is present, verify that its entry count matches
  `--size-groups`;
- if explicit `--size-targets` are also present and are being used as labels,
  verify their count matches `--size-groups`;
- if only `--size-groups` is supplied, no five-entry target list is required;
  automatic `G<count>` labels are valid.

Do not force custom labels or series lengths back to the manuscript defaults
merely for collator compatibility; current Binning discovery is label-flexible.

---

## 10. Handoff to Binning

After partition generation, do not continue implementing Binning analysis here.

The next workflow belongs to the Binning repository and begins from dataset
directories such as:

```text
Binning/data/Soil_1/Soil_Binning/
```

and:

```text
Binning/data/Soil_1/Soil_Groups_1.2M/
```

A future AI helper should consult the Binning `README.md` and `AI_README.md`
before guiding the user through search/downstream/correction stages.

---

## 11. Handoff to Targeted Peptide Removal

Targeted Peptide Removal uses FASTA-processing outputs including:

```text
clustered protein FASTA
clustered peptide FASTA
peptide_protein_digest_map.tsv
```

Do not duplicate the digest code in that repository.

The Targeted Peptide Removal workflow owns:

- peptide similarity evaluation;
- conservative conflict culling;
- threshold-series removal databases;
- peptide/protein downstream analysis;
- peptide/protein gain/loss collation.

A planned helper will make staging/symlinking the shared protein FASTA and
digest map into peptide-removal threshold runs easier. Until that helper exists,
do not invent provenance for those resources.

---

## 12. Troubleshooting approach

When a command fails:

1. identify which public launcher failed;
2. read the actual error;
3. confirm the input file exists;
4. confirm whether CD-HIT is truly required;
5. inspect output/destination paths;
6. use `--dry-run` where available;
7. avoid changing scientific defaults until the operational problem is
   understood.

Useful checks:

```bash
bash -n scripts/pipeline/prepare_peptide_fasta.sh
```

```bash
bash -n scripts/pipeline/prepare_binning_databases.sh
```

```bash
python3 -m py_compile scripts/fasta/peptide_fasta_builder.py
```

```bash
python3 -m py_compile scripts/binning/peptide_bin_pair_builder.py
```

---

## 13. Maintainer invariants

Before changing scientific behavior:

1. syntax-check Python and shell files;
2. exercise public `--help` interfaces;
3. test an already-clustered FASTA without CD-HIT;
4. test clustering when CD-HIT is available;
5. test peptide digestion and digest-map generation;
6. test a 9-group/36-pair primary Binning build;
7. test size-series group-count calculation;
8. test exact `--size-groups` mode;
9. inspect output directory names against Binning expectations.

Do not normalize away changes in peptide identities, peptide-to-protein mappings,
partition membership, or group counts.

---

## 14. Documentation responsibilities

If behavior changes, update:

- `README.md`;
- `AI_README.md`;
- the related Binning documentation if output layouts or names change;
- the related Targeted Peptide Removal documentation if shared FASTA resources
  change.

When the GitHub repositories are created, add all three repository links and
tell future AI helpers to inspect all three related projects before creating
new cross-workflow code.
