#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
cat <<'USAGE'
Usage:
  stage_shared_fasta_resources.sh --condition NAME [options]

Purpose:
  Stage shared FASTA-processing resources into existing Targeted Peptide
  Removal baseline/threshold run directories using relative symbolic links.

Required:
  --condition NAME           Condition/run family name, e.g. Fecal, Ocean, Soil.

Source options:
  --fasta-data-root DIR      FASTA_Processing data root.
                             Default: sibling FASTA_Processing/data
  --protein-fasta FILE       Explicit clustered protein FASTA/FAA source.
  --digest-map FILE          Explicit peptide_protein_digest_map.tsv source.

Destination options:
  --data-root DIR            Targeted Peptide Removal data root.
                             Default: <Targeted_Peptide_Removal>/data
  --replace-links            Replace stale/wrong symbolic links only.
                             Regular files are never overwritten.
  --dry-run                  Print planned actions without modifying files.
  -h, --help                 Show this help.

Default source discovery:
  <FASTA-data-root>/<Condition>/protein/*_clustered.fasta
  <FASTA-data-root>/<Condition>/peptide/peptide_protein_digest_map.tsv

Destination discovery:
  Existing directories named exactly <Condition> or <Condition>_* under
  --data-root are staged. The script does not create missing run directories.

Each run receives:
  <run>/peptide_protein_digest_map.tsv
  <run>/<protein FASTA basename>

Links are relative, so sibling repositories can be moved together without
baking machine-specific absolute paths into the run directories.
USAGE
}

die(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd -P)"
COLLECTION_ROOT="$(cd -- "$PROJECT_ROOT/.." && pwd -P)"

CONDITION=""
DATA_ROOT="$PROJECT_ROOT/data"
FASTA_DATA_ROOT="$COLLECTION_ROOT/FASTA_Processing/data"
PROTEIN_FASTA=""
DIGEST_MAP=""
REPLACE_LINKS=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --condition) CONDITION="${2:?Missing value for --condition}"; shift 2 ;;
    --data-root) DATA_ROOT="${2:?Missing value for --data-root}"; shift 2 ;;
    --fasta-data-root) FASTA_DATA_ROOT="${2:?Missing value for --fasta-data-root}"; shift 2 ;;
    --protein-fasta) PROTEIN_FASTA="${2:?Missing value for --protein-fasta}"; shift 2 ;;
    --digest-map) DIGEST_MAP="${2:?Missing value for --digest-map}"; shift 2 ;;
    --replace-links) REPLACE_LINKS=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$CONDITION" ]] || die "--condition is required."
[[ "$CONDITION" != *'/'* ]] || die "--condition must not contain '/'."
[[ -d "$DATA_ROOT" ]] || die "Targeted Peptide Removal data root not found: $DATA_ROOT"

SOURCE_ROOT="$FASTA_DATA_ROOT/$CONDITION"

if [[ -z "$PROTEIN_FASTA" ]]; then
  [[ -d "$SOURCE_ROOT/protein" ]] || die "Protein resource directory not found: $SOURCE_ROOT/protein"
  mapfile -t protein_candidates < <(find "$SOURCE_ROOT/protein" -maxdepth 1 -type f \
    \( -name '*_clustered.fasta' -o -name '*_clustered.fa' -o -name '*_clustered.faa' \) | sort)
  [[ ${#protein_candidates[@]} -eq 1 ]] || die "Expected exactly one clustered protein FASTA/FAA in $SOURCE_ROOT/protein; found ${#protein_candidates[@]}. Use --protein-fasta explicitly."
  PROTEIN_FASTA="${protein_candidates[0]}"
fi

if [[ -z "$DIGEST_MAP" ]]; then
  DIGEST_MAP="$SOURCE_ROOT/peptide/peptide_protein_digest_map.tsv"
fi

[[ -f "$PROTEIN_FASTA" ]] || die "Protein FASTA not found: $PROTEIN_FASTA"
[[ -f "$DIGEST_MAP" ]] || die "Digest map not found: $DIGEST_MAP"
PROTEIN_FASTA="$(cd -- "$(dirname -- "$PROTEIN_FASTA")" && pwd -P)/$(basename -- "$PROTEIN_FASTA")"
DIGEST_MAP="$(cd -- "$(dirname -- "$DIGEST_MAP")" && pwd -P)/$(basename -- "$DIGEST_MAP")"

mapfile -t RUN_DIRS < <(find "$DATA_ROOT" -mindepth 1 -maxdepth 1 -type d \
  \( -name "$CONDITION" -o -name "${CONDITION}_*" \) ! -name results ! -name results_regression | sort)
[[ ${#RUN_DIRS[@]} -gt 0 ]] || die "No existing run directories found for $CONDITION under $DATA_ROOT"

relative_target() {
  python3 - "$1" "$2" <<'PY'
import os, sys
source, dest_dir = sys.argv[1:]
print(os.path.relpath(source, dest_dir))
PY
}

stage_link() {
  local source="$1"
  local dest="$2"
  local label="$3"
  local dest_dir
  dest_dir="$(dirname -- "$dest")"

  if [[ -L "$dest" ]]; then
    local current resolved_source resolved_dest
    current="$(readlink -- "$dest")"
    resolved_source="$(python3 - "$source" <<'PY'
import os, sys
print(os.path.realpath(sys.argv[1]))
PY
)"
    resolved_dest="$(python3 - "$dest" <<'PY'
import os, sys
print(os.path.realpath(sys.argv[1]))
PY
)"
    if [[ "$resolved_dest" == "$resolved_source" ]]; then
      printf '  OK      %-12s %s\n' "$label" "$dest"
      return
    fi
    if [[ "$REPLACE_LINKS" -ne 1 ]]; then
      die "Existing symlink points elsewhere: $dest -> $current (use --replace-links to replace symlinks only)"
    fi
    printf '  REPLACE %-12s %s\n' "$label" "$dest"
    [[ "$DRY_RUN" -eq 1 ]] || rm -- "$dest"
  elif [[ -e "$dest" ]]; then
    printf '  EXISTS  %-12s %s (regular file left unchanged)\n' "$label" "$dest"
    return
  fi

  local rel
  rel="$(relative_target "$source" "$dest_dir")"
  printf '  LINK    %-12s %s -> %s\n' "$label" "$dest" "$rel"
  if [[ "$DRY_RUN" -eq 0 ]]; then
    ln -s -- "$rel" "$dest"
  fi
}

log "Source protein FASTA: $PROTEIN_FASTA"
log "Source digest map:    $DIGEST_MAP"
log "Runs discovered:      ${#RUN_DIRS[@]}"

protein_name="$(basename -- "$PROTEIN_FASTA")"

for run_dir in "${RUN_DIRS[@]}"; do
  echo
  echo "Run: $run_dir"

  mapfile -t existing_fastas < <(find "$run_dir" -maxdepth 1 -type f \
      \( -name '*.fasta' -o -name '*.fa' -o -name '*.faa' \) -print | sort)
  mapfile -t existing_fasta_links < <(find "$run_dir" -maxdepth 1 -type l \
      \( -name '*.fasta' -o -name '*.fa' -o -name '*.faa' \) -print | sort)
  if (( ${#existing_fastas[@]} + ${#existing_fasta_links[@]} > 1 )); then
    die "Run already contains multiple top-level protein FASTA/FAA candidates: $run_dir"
  fi
  if (( ${#existing_fastas[@]} + ${#existing_fasta_links[@]} == 1 )); then
    existing="${existing_fastas[0]:-${existing_fasta_links[0]}}"
    if [[ "$(basename -- "$existing")" != "$protein_name" ]]; then
      die "Run already contains a differently named top-level protein FASTA/FAA: $existing"
    fi
  fi

  stage_link "$PROTEIN_FASTA" "$run_dir/$protein_name" "protein"
  stage_link "$DIGEST_MAP" "$run_dir/peptide_protein_digest_map.tsv" "digest-map"
done

echo
if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry run only; no links were created."
else
  echo "Shared FASTA resources staged successfully."
fi
