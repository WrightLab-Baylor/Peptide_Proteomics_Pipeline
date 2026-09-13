# Proteomics FDR Adjustment Workflows

This repository collection contains three related workflow forks:

1. **FASTA_Processing**
2. **Binning**
3. **Targeted_Peptide_Removal**

## Start here

**FASTA_Processing is step 1 for both downstream workflows.**

Use it first to:

- cluster the protein FASTA;
- digest proteins into the peptide FASTA;
- create `peptide_protein_digest_map.tsv`; and
- when needed, create the Individual/Pairwise databases used by Binning.

After FASTA preparation, continue with either:

```text
FASTA_Processing
├── Binning
└── Targeted_Peptide_Removal
```

### Binning

Use **Binning** for peptide-database partitioning analysis, classical-guided
culling, partition-level FDR comparison, correctional-formula processing, and
peptide/protein result collation.

See:

```text
Binning/README.md
```

for the complete workflow and usage instructions.

### Targeted Peptide Removal

Use **Targeted_Peptide_Removal** to evaluate peptide target/decoy sequence
similarity, generate targeted-removal database series, run peptide/protein
downstream processing, and collate peptide/protein gains and losses.

See:

```text
Targeted_Peptide_Removal/README.md
```

for the complete workflow and usage instructions.

### FASTA Processing

See:

```text
FASTA_Processing/README.md
```

for clustering, digestion, digest-map generation, and Binning database
preparation.

## Shared environment

The three workflows use one shared Conda environment:

```bash
conda env create -f environment.yml
conda activate proteomics-fdr-tools
```

The detailed README in each fork is the source of truth for that workflow's
commands, expected directory layout, inputs, outputs, and scientific defaults.
