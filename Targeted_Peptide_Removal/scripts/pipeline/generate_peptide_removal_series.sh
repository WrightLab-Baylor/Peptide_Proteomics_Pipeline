#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
cat <<'USAGE'
Usage:
  generate_peptide_removal_series.sh \
    --condition NAME \
    --peptide-fasta FILE \
    --conflicts-tsv FILE \
    [options]

Purpose:
  Generate the standard targeted peptide-removal threshold series from one
  evaluator conflicts.tsv. Each threshold receives its own run directory under
  the Targeted_Peptide_Removal data root.

Required:
  --condition NAME       Condition/run family name, e.g. Fecal, Ocean, Soil.
  --peptide-fasta FILE   Unculled target peptide FASTA.
  --conflicts-tsv FILE   conflicts.tsv from peptide_fasta_evaluator.py.

Options:
  --data-root PATH       Output data root. Default: <Targeted_Peptide_Removal>/data
  --thresholds LIST      Comma-separated thresholds. Default: 1.0,0.9,0.8,0.7,0.6,0.5
  --overwrite            Rebuild generator-managed files only when the run has
                         no search/downstream outputs.
  --dry-run              Print actions without creating or modifying files.
  -h, --help             Show this help.

Generated layout:
  <data-root>/<Condition>_<threshold>/
  ├── data/
  ├── database/
  │   └── <Condition>_<threshold>.fasta
  └── culling/
      ├── cull_report.tsv
      ├── removed_peptides.tsv
      └── remaining_conflicts_after_cull.tsv

Notes:
  The empty data/ directory is created for acquisition/search inputs. Protein
  FASTA and peptide_protein_digest_map.tsv resources needed later are not guessed
  or copied by this launcher; stage them from the upstream FASTA-processing run.
USAGE
}

die(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd -P)"
CULLER="$PROJECT_ROOT/scripts/removal/conflict_culler_fast.py"

CONDITION=""
PEPTIDE_FASTA=""
CONFLICTS_TSV=""
DATA_ROOT="$PROJECT_ROOT/data"
THRESHOLDS_CSV="1.0,0.9,0.8,0.7,0.6,0.5"
OVERWRITE="False"
DRY_RUN="False"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --condition) CONDITION="${2:?Missing value for --condition}"; shift 2 ;;
    --peptide-fasta) PEPTIDE_FASTA="${2:?Missing value for --peptide-fasta}"; shift 2 ;;
    --conflicts-tsv) CONFLICTS_TSV="${2:?Missing value for --conflicts-tsv}"; shift 2 ;;
    --data-root) DATA_ROOT="${2:?Missing value for --data-root}"; shift 2 ;;
    --thresholds) THRESHOLDS_CSV="${2:?Missing value for --thresholds}"; shift 2 ;;
    --overwrite) OVERWRITE="True"; shift ;;
    --dry-run) DRY_RUN="True"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$CONDITION" ]] || die "--condition is required."
[[ "$CONDITION" != *'/'* ]] || die "--condition must not contain '/'."
[[ -n "$PEPTIDE_FASTA" ]] || die "--peptide-fasta is required."
[[ -n "$CONFLICTS_TSV" ]] || die "--conflicts-tsv is required."
[[ -f "$PEPTIDE_FASTA" ]] || die "Peptide FASTA not found: $PEPTIDE_FASTA"
[[ -f "$CONFLICTS_TSV" ]] || die "conflicts.tsv not found: $CONFLICTS_TSV"
[[ -f "$CULLER" ]] || die "Conflict culler not found: $CULLER"

if [[ "$DRY_RUN" == "True" ]]; then
  printf 'Would create data root: %q\n' "$DATA_ROOT"
else
  mkdir -p "$DATA_ROOT"
fi

IFS=',' read -r -a THRESHOLDS <<< "$THRESHOLDS_CSV"
((${#THRESHOLDS[@]} > 0)) || die "No thresholds supplied."

SUMMARY="$DATA_ROOT/${CONDITION}_peptide_removal_series_summary.tsv"
if [[ "$DRY_RUN" == "False" ]]; then
  printf 'Condition\tThreshold\tRunDirectory\tDatabaseFASTA\tTotalRecords\tQualifyingEdges\tRemovedPeptides\tKeptRecords\tRemainingEdges\n' > "$SUMMARY"
fi

for threshold in "${THRESHOLDS[@]}"; do
  threshold="${threshold//[[:space:]]/}"
  [[ -n "$threshold" ]] || continue

  python3 -c 'import sys; x=float(sys.argv[1]); assert 0.5 <= x <= 1.0, "threshold must be between 0.5 and 1.0"' "$threshold" || die "Invalid threshold: $threshold"

  label="${CONDITION}_${threshold}"
  run_dir="$DATA_ROOT/$label"
  cull_dir="$run_dir/culling"
  db_dir="$run_dir/database"
  raw_dir="$run_dir/data"
  db_fasta="$db_dir/${label}.fasta"

  if [[ -d "$run_dir" ]]; then
    processed="False"
    for p in SICdir SICs fdr_esti fdr_filtered peptide_crosstab.tsv peptide_crosstab_annotated.tsv protein_crosstab_annotated.tsv; do
      [[ -e "$run_dir/$p" ]] && processed="True"
    done
    [[ "$processed" == "False" ]] || die "Refusing to overwrite processed run directory: $run_dir"
    [[ "$OVERWRITE" == "True" ]] || die "Run directory already exists: $run_dir (use --overwrite only before search/downstream processing)"
  fi

  log "Generating $label"
  if [[ "$DRY_RUN" == "True" ]]; then
    printf 'Would create: %s %s %s\n' "$raw_dir" "$db_dir" "$cull_dir"
    printf 'Would run: python3 %q --peptide_fasta %q --conflicts_tsv %q --out_dir %q --min_jaccard %q --min_lcs %q\n' \
      "$CULLER" "$PEPTIDE_FASTA" "$CONFLICTS_TSV" "$cull_dir" "$threshold" "$threshold"
    printf 'Would copy: %q -> %q\n' "$cull_dir/filtered_target.fa" "$db_fasta"
    continue
  fi

  rm -rf -- "$cull_dir" "$db_dir"
  mkdir -p "$raw_dir" "$db_dir" "$cull_dir"

  python3 "$CULLER" \
    --peptide_fasta "$PEPTIDE_FASTA" \
    --conflicts_tsv "$CONFLICTS_TSV" \
    --out_dir "$cull_dir" \
    --min_jaccard "$threshold" \
    --min_lcs "$threshold"

  [[ -f "$cull_dir/filtered_target.fa" ]] || die "Culler did not produce filtered_target.fa for $label"
  cp -p -- "$cull_dir/filtered_target.fa" "$db_fasta"

  report="$cull_dir/cull_report.tsv"
  [[ -f "$report" ]] || die "Missing cull report: $report"
  metric(){ awk -F'\t' -v key="$1" '$1==key {print $2; exit}' "$report" | tr -d '\r'; }
  total_records="$(metric total_records_in_fasta)"
  qualifying_edges="$(metric qualifying_edges_seen)"
  removed_peptides="$(metric removed_peptides_unique)"
  kept_records="$(metric kept_records_in_fasta)"
  remaining_edges="$(metric edges_remaining_after_cull)"
  [[ "${remaining_edges:-}" == "0" ]] || die "$label retains qualifying conflicts: ${remaining_edges:-unknown}"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$CONDITION" "$threshold" "$run_dir" "$db_fasta" "$total_records" \
    "$qualifying_edges" "$removed_peptides" "$kept_records" "$remaining_edges" >> "$SUMMARY"
done

if [[ "$DRY_RUN" == "True" ]]; then
  log "Dry run complete; no files were created."
else
  log "Peptide-removal series complete: $SUMMARY"
fi
