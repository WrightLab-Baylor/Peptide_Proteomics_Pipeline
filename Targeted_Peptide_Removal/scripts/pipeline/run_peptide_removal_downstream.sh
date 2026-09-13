#!/usr/bin/env bash
set -Eeuo pipefail

usage(){ cat <<'USAGE'
Usage:
  run_peptide_removal_downstream.sh --run-dir PATH [options]

Purpose:
  Run the validated peptide-FASTA downstream workflow for one baseline or
  peptide-removal run directory.

Required:
  --run-dir PATH
      Run directory containing SICdir/, database/, peptide_protein_digest_map.tsv,
      and one top-level protein FASTA/FAA.

Validated defaults:
  --peptide-filter MSMSScore
  --stringency 8
  --peptide-filter2 FDR
  --stringency2 0.05
  --highest-intensity-only True
  --unique-to-cluster False
  --uniqueness-mode all_matches
  --minimum-peptide-count 1
  --peptide-counting-mode per_sample
  --cluster-sums True
  --rollup-type sum
  --rrollup-mode None

Other options:
  --cdhit-path PATH      Explicit cd-hit executable if it is not on PATH.
  --dry-run              Print commands without executing them.
  -h, --help             Show this help.
USAGE
}

die(){ printf 'ERROR: %s\n' "$*" >&2; exit 1; }
log(){ printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
quote(){ printf '%q ' "$@"; }

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd -P)"
TOOLS="$PROJECT_ROOT/scripts/downstream"

RUN_DIR=""
PEPTIDE_FILTER="MSMSScore"
STRINGENCY="8"
PEPTIDE_FILTER2="FDR"
STRINGENCY2="0.05"
HIGHEST_INTENSITY_ONLY="True"
UNIQUE_TO_CLUSTER="False"
UNIQUENESS_MODE="all_matches"
MINIMUM_PEPTIDE_COUNT="1"
PEPTIDE_COUNTING_MODE="per_sample"
CLUSTER_SUMS="True"
ROLLUP_TYPE="sum"
RROLLUP_MODE="None"
CDHIT_PATH=""
DRY_RUN="False"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="${2:?Missing value for --run-dir}"; shift 2 ;;
    --peptide-filter) PEPTIDE_FILTER="${2:?Missing value}"; shift 2 ;;
    --stringency) STRINGENCY="${2:?Missing value}"; shift 2 ;;
    --peptide-filter2) PEPTIDE_FILTER2="${2:?Missing value}"; shift 2 ;;
    --stringency2) STRINGENCY2="${2:?Missing value}"; shift 2 ;;
    --highest-intensity-only) HIGHEST_INTENSITY_ONLY="${2:?Missing value}"; shift 2 ;;
    --unique-to-cluster) UNIQUE_TO_CLUSTER="${2:?Missing value}"; shift 2 ;;
    --uniqueness-mode) UNIQUENESS_MODE="${2:?Missing value}"; shift 2 ;;
    --minimum-peptide-count) MINIMUM_PEPTIDE_COUNT="${2:?Missing value}"; shift 2 ;;
    --peptide-counting-mode) PEPTIDE_COUNTING_MODE="${2:?Missing value}"; shift 2 ;;
    --cluster-sums) CLUSTER_SUMS="${2:?Missing value}"; shift 2 ;;
    --rollup-type) ROLLUP_TYPE="${2:?Missing value}"; shift 2 ;;
    --rrollup-mode) RROLLUP_MODE="${2:?Missing value}"; shift 2 ;;
    --cdhit-path) CDHIT_PATH="${2:?Missing value}"; shift 2 ;;
    --dry-run) DRY_RUN="True"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$RUN_DIR" ]] || die "--run-dir is required."
[[ -d "$RUN_DIR" ]] || die "Run directory not found: $RUN_DIR"
RUN_DIR="$(cd -- "$RUN_DIR" && pwd -P)"

required_scripts=(prepare_SICs.py fdr_estimator.py peptide_map_expander.py true_pep_filter.py peptide_crosstab_generator.py database_cleaner.py fasta_cluster_generator_peptide.py peptide_crosstab_annotator.py peptide_protein_mapbuilder.py protein_rollup.py cluster_crosstab_generator.py)
for s in "${required_scripts[@]}"; do [[ -f "$TOOLS/$s" ]] || die "Missing downstream tool: $TOOLS/$s"; done
[[ -d "$RUN_DIR/SICdir" ]] || die "Missing SICdir/: $RUN_DIR/SICdir"
[[ -d "$RUN_DIR/database" ]] || die "Missing database/: $RUN_DIR/database"
[[ -f "$RUN_DIR/peptide_protein_digest_map.tsv" ]] || die "Missing peptide_protein_digest_map.tsv in $RUN_DIR"

mapfile -t protein_fastas < <(find "$RUN_DIR" -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.faa' \) -print | sort)
[[ ${#protein_fastas[@]} -eq 1 ]] || die "Expected exactly one top-level protein FASTA/FAA in $RUN_DIR; found ${#protein_fastas[@]}"

run(){
  log "RUN: $(quote "$@")"
  [[ "$DRY_RUN" == "True" ]] || "$@"
}

log "Peptide-removal downstream"
log "Run directory: $RUN_DIR"
log "Tools:         $TOOLS"
log "Dry run:       $DRY_RUN"

cd "$RUN_DIR"

log "Step 1/11: Preparing SICs"
run python3 "$TOOLS/prepare_SICs.py" -i SICdir/ -o SICs

log "Step 2/11: Estimating FDR"
run python3 "$TOOLS/fdr_estimator.py" -i SICs -o fdr_esti

log "Step 3/11: Expanding peptide-protein mappings"
run python3 "$TOOLS/peptide_map_expander.py" -w "$RUN_DIR"

log "Step 4/11: Filtering peptide-spectrum matches"
run python3 "$TOOLS/true_pep_filter.py" -i fdr_esti -o fdr_filtered -f "$PEPTIDE_FILTER" -t "$STRINGENCY" --filter2 "$PEPTIDE_FILTER2" --threshold2 "$STRINGENCY2"

log "Step 5/11: Generating peptide crosstab"
run python3 "$TOOLS/peptide_crosstab_generator.py" -i fdr_filtered -o peptide_crosstab.tsv --only-highest-intensity "$HIGHEST_INTENSITY_ONLY"

log "Step 6/11: Cleaning temporary database files"
run python3 "$TOOLS/database_cleaner.py"

log "Step 7/11: Ensuring protein FASTA has cluster annotations"
cluster_cmd=(python3 "$TOOLS/fasta_cluster_generator_peptide.py")
[[ -n "$CDHIT_PATH" ]] && cluster_cmd+=(--cdhit-path "$CDHIT_PATH")
run "${cluster_cmd[@]}"

mapfile -t protein_fastas_after < <(find "$RUN_DIR" -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.faa' \) -print | sort)
if [[ "$DRY_RUN" == "False" ]]; then
  [[ ${#protein_fastas_after[@]} -eq 1 ]] || die "Expected exactly one top-level protein FASTA/FAA after clustering; found ${#protein_fastas_after[@]}"
  PROTEIN_FASTA="${protein_fastas_after[0]}"
else
  PROTEIN_FASTA="${protein_fastas[0]}"
fi

log "Step 8/11: Annotating peptide crosstab"
run python3 "$TOOLS/peptide_crosstab_annotator.py" -i peptide_crosstab.tsv -f "$PROTEIN_FASTA" -o peptide_crosstab_annotated.tsv -unique_to_cluster "$UNIQUE_TO_CLUSTER"

log "Step 9/11: Building peptide-protein map"
run python3 "$TOOLS/peptide_protein_mapbuilder.py" -i peptide_crosstab_annotated.tsv -db . -p "$MINIMUM_PEPTIDE_COUNT" --mode "$UNIQUENESS_MODE" -f "$PEPTIDE_COUNTING_MODE"

log "Step 10/11: Generating protein rollup"
run python3 "$TOOLS/protein_rollup.py" -i peptide_crosstab_annotated.tsv -o protein_crosstab_annotated.tsv --mode "$UNIQUENESS_MODE" --rollup "$ROLLUP_TYPE" --rrollup "$RROLLUP_MODE"

log "Step 11/11: Generating cluster rollup when requested"
if [[ "${CLUSTER_SUMS,,}" == "true" ]]; then
  run python3 "$TOOLS/cluster_crosstab_generator.py" -i peptide_crosstab_annotated.tsv -o clustered_crosstab_annotated.tsv --mode "$UNIQUENESS_MODE" --rollup "$ROLLUP_TYPE" --rrollup "$RROLLUP_MODE"
else
  log "Cluster rollup disabled."
fi

log "Peptide-removal downstream completed successfully."
