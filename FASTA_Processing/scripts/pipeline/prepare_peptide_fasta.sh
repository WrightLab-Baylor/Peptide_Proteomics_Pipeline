#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CLUSTER_SCRIPT="${REPO_ROOT}/scripts/fasta/fasta_cluster_generator.py"
DIGEST_SCRIPT="${REPO_ROOT}/scripts/fasta/peptide_fasta_builder.py"

CONDITION=""
PROTEIN_FASTA=""
OUTPUT_ROOT="${REPO_ROOT}/data"
IDENTITY="0.95"
WORD_LENGTH="5"
THREADS="0"
MEMORY_MB="0"
MIN_LEN="6"
MAX_LEN="50"
CDHIT_PATH=""
REMOVE_TDA_OVERLAPS=0
DRY_RUN=0

usage() {
cat <<'EOF'
Usage:
  prepare_peptide_fasta.sh --condition NAME --protein-fasta FILE [options]

Cluster a protein FASTA with CD-HIT, annotate protein headers with cluster IDs,
then digest the clustered protein FASTA into a peptide FASTA and
peptide_protein_digest_map.tsv.

Required:
  --condition NAME           Condition/dataset name used for the output folder.
  --protein-fasta FILE       Input protein FASTA/FAA.

Output:
  By default, products are written under:
    <FASTA_Processing>/data/<Condition>/
      protein/<Condition>_clustered.fasta
      peptide/<Condition>_clustered_peptide.fasta
      peptide/peptide_protein_digest_map.tsv
      peptide/digest_stats.txt
      peptide/removed_peptides.tsv

Options:
  --output-root DIR          Output root (default: <FASTA_Processing>/data).
  --identity FLOAT           CD-HIT identity threshold (default: 0.95).
  --word-length INT          CD-HIT word length (default: 5).
  --threads INT              CD-HIT threads; 0 = all CPUs (default: 0).
  --memory-mb INT            CD-HIT memory limit in MB; 0 = unlimited (default: 0).
  --cdhit-path FILE          Explicit cd-hit executable.
  --min-len INT              Minimum peptide length (default: 6).
  --max-len INT              Maximum peptide length (default: 50).
  --remove-tda-overlaps      Enable optional target/decoy-overlap removal.
                             Disabled by default.
  --dry-run                  Print commands without executing them.
  -h, --help                 Show this help.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --condition) CONDITION="$2"; shift 2 ;;
        --protein-fasta) PROTEIN_FASTA="$2"; shift 2 ;;
        --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
        --identity) IDENTITY="$2"; shift 2 ;;
        --word-length) WORD_LENGTH="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --memory-mb) MEMORY_MB="$2"; shift 2 ;;
        --cdhit-path) CDHIT_PATH="$2"; shift 2 ;;
        --min-len) MIN_LEN="$2"; shift 2 ;;
        --max-len) MAX_LEN="$2"; shift 2 ;;
        --remove-tda-overlaps) REMOVE_TDA_OVERLAPS=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$CONDITION" ]] || { echo "ERROR: --condition is required." >&2; exit 2; }
[[ -n "$PROTEIN_FASTA" ]] || { echo "ERROR: --protein-fasta is required." >&2; exit 2; }
[[ -f "$PROTEIN_FASTA" ]] || { echo "ERROR: protein FASTA not found: $PROTEIN_FASTA" >&2; exit 2; }

CONDITION_DIR="${OUTPUT_ROOT}/${CONDITION}"
PROTEIN_DIR="${CONDITION_DIR}/protein"
PEPTIDE_DIR="${CONDITION_DIR}/peptide"
CLUSTERED_FASTA="${PROTEIN_DIR}/${CONDITION}_clustered.fasta"

cluster_cmd=(
    python3 "$CLUSTER_SCRIPT"
    --input "$PROTEIN_FASTA"
    --output "$CLUSTERED_FASTA"
    --identity "$IDENTITY"
    --word-length "$WORD_LENGTH"
    --threads "$THREADS"
    --memory "$MEMORY_MB"
)
if [[ -n "$CDHIT_PATH" ]]; then
    cluster_cmd+=(--cdhit-path "$CDHIT_PATH")
fi

digest_cmd=(
    python3 "$DIGEST_SCRIPT"
    --input_fasta "$CLUSTERED_FASTA"
    --out_dir "$PEPTIDE_DIR"
    --min_len "$MIN_LEN"
    --max_len "$MAX_LEN"
)
if [[ "$REMOVE_TDA_OVERLAPS" -eq 1 ]]; then
    digest_cmd+=(--remove_tda_overlaps)
fi

printf 'Protein clustering:\n  '
printf '%q ' "${cluster_cmd[@]}"
printf '\nPeptide digestion:\n  '
printf '%q ' "${digest_cmd[@]}"
printf '\n'

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "Dry run only; no commands executed."
    exit 0
fi

mkdir -p "$PROTEIN_DIR" "$PEPTIDE_DIR"

"${cluster_cmd[@]}"
"${digest_cmd[@]}"

echo
echo "FASTA preprocessing complete."
echo "Condition root:      $CONDITION_DIR"
echo "Clustered protein:   $CLUSTERED_FASTA"
echo "Peptide directory:   $PEPTIDE_DIR"
