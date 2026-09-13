# AI Helper Guide — Proteomics FDR Adjustment Workflows

This is the top-level handoff guide for AI assistants working with the three
related workflow forks:

```text
FASTA_Processing
Binning
Targeted_Peptide_Removal
```

## Critical workflow order

**FASTA_Processing is step 1 before either Binning or Targeted Peptide Removal.**

Do not guide a new user directly into Binning or Targeted Peptide Removal unless
the required FASTA-processing outputs already exist.

FASTA_Processing creates the shared upstream resources, including:

```text
clustered protein FASTA
peptide FASTA
peptide_protein_digest_map.tsv
```

and, for Binning, can also create Individual/Pairwise partition databases.

The normal flow is therefore:

```text
FASTA_Processing
    ├── Binning
    └── Targeted_Peptide_Removal
```

## Before helping with a workflow

Read the detailed AI and human README files for the relevant forks:

```text
FASTA_Processing/README.md
FASTA_Processing/AI_README.md

Binning/README.md
Binning/AI_README.md

Targeted_Peptide_Removal/README.md
Targeted_Peptide_Removal/AI_README.md
```

If a task crosses repository boundaries, inspect all relevant README files
before suggesting commands or code changes.

Do not duplicate upstream functionality inside a downstream fork merely because
an input is missing. Instead, identify which workflow should create that input
and guide the user through the correct preceding step.

## Shared environment

All three forks use the shared top-level:

```text
environment.yml
```

Prefer this shared environment over creating independent per-fork environments.

## Guidance style

The intended user may be a new graduate student with little command-line
experience. Prefer:

- one command at a time;
- explicit working directories and paths;
- `--dry-run` where supported;
- existing public launchers over direct low-level script invocation;
- identifying missing upstream inputs before modifying code;
- preserving the scientific defaults and invariants documented in each fork.

The fork-specific `AI_README.md` files contain the detailed operational and
maintenance constraints and should be treated as the primary source of truth
for each workflow.
