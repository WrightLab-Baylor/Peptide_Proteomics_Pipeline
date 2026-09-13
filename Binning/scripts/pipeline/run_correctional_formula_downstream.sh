#!/usr/bin/env bash
#
# Derive an empirical FDR correction curve, apply it to a classical *_Full run,
# and generate peptide/protein rollup outputs inside the related Group/Pairwise
# downstream directory.
#
# The source *_Full directory is treated as read-only. Original fdr_esti files
# are never overwritten.

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_VERSION="1.0.2"

usage() {
    cat <<'USAGE'
Usage:
  run_correctional_formula_downstream.sh \
    --downstream-dir PATH \
    --full-dir PATH \
    [options]

Required:
  --downstream-dir PATH   Existing *_Downstream directory containing the
                          classical-vs-partitioned audit table.
  --full-dir PATH         Related classical *_Full directory containing
                          fdr_esti/*_withsyn.tsv and one top-level FASTA/FAA.

Normally inferred:
  --shared-tsv PATH       Shared comparison TSV. Default: uniquely locate
                          classic_vs_*/classic_vs_binned_scan_peptide_fdr_wide.tsv
                          below --downstream-dir.
  --condition-name NAME   Formula condition. Default: basename of the parent
                          partition directory containing *_Downstream.
  --project-root PATH     Binning software root. Default: inferred from this
                          shell location, with the current working directory
                          retained as a compatibility fallback.

Output locations below --downstream-dir:
  correctional_formula/<condition-name>/
  corrected_fdr_esti/
  protein_rollup/
  logs/correctional_formula_pipeline/

Downstream parameter defaults:
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
  --formula-dir PATH      Override formula output directory.
  --corrected-dir PATH    Override corrected fdr_esti output directory.
  --rollup-dir PATH       Override protein rollup output directory.
  --thresholds LIST       Application diagnostic cutoffs. Default: 0.01,0.05
  --overwrite             Remove and regenerate managed output directories.
  --skip-derive           Reuse an existing correction_formula.json.
  --skip-apply            Reuse existing corrected *_withsyn.tsv files.
  --dry-run               Print commands without executing them.
  -h, --help              Show this help.
USAGE
}

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

quote_cmd() {
    printf '%q ' "$@"
    printf '\n'
}

run_cmd() {
    log "RUN: $(quote_cmd "$@")"
    if [[ "$DRY_RUN" == "True" ]]; then
        return 0
    fi
    "$@"
}

require_file() {
    [[ -f "$1" ]] || die "Required file not found: $1"
}

require_dir() {
    [[ -d "$1" ]] || die "Required directory not found: $1"
}

normalize_bool() {
    case "${1,,}" in
        true|t|1|yes|y) printf 'True\n' ;;
        false|f|0|no|n) printf 'False\n' ;;
        *) die "Expected a boolean value, received: $1" ;;
    esac
}

# Required/inferred paths.
DOWNSTREAM_DIR=""
FULL_DIR=""
SHARED_TSV=""
CONDITION_NAME=""
PROJECT_ROOT=""
FORMULA_DIR=""
CORRECTED_DIR=""
ROLLUP_DIR=""

# Standard defaults requested for this workflow.
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
THRESHOLDS="0.01,0.05"

OVERWRITE="False"
SKIP_DERIVE="False"
SKIP_APPLY="False"
DRY_RUN="False"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --downstream-dir) DOWNSTREAM_DIR="${2:?Missing value for --downstream-dir}"; shift 2 ;;
        --full-dir) FULL_DIR="${2:?Missing value for --full-dir}"; shift 2 ;;
        --shared-tsv) SHARED_TSV="${2:?Missing value for --shared-tsv}"; shift 2 ;;
        --condition-name) CONDITION_NAME="${2:?Missing value for --condition-name}"; shift 2 ;;
        --project-root) PROJECT_ROOT="${2:?Missing value for --project-root}"; shift 2 ;;
        --formula-dir) FORMULA_DIR="${2:?Missing value for --formula-dir}"; shift 2 ;;
        --corrected-dir) CORRECTED_DIR="${2:?Missing value for --corrected-dir}"; shift 2 ;;
        --rollup-dir) ROLLUP_DIR="${2:?Missing value for --rollup-dir}"; shift 2 ;;
        --thresholds) THRESHOLDS="${2:?Missing value for --thresholds}"; shift 2 ;;
        --peptide-filter) PEPTIDE_FILTER="${2:?Missing value}"; shift 2 ;;
        --stringency) STRINGENCY="${2:?Missing value}"; shift 2 ;;
        --peptide-filter2) PEPTIDE_FILTER2="${2:?Missing value}"; shift 2 ;;
        --stringency2) STRINGENCY2="${2:?Missing value}"; shift 2 ;;
        --highest-intensity-only) HIGHEST_INTENSITY_ONLY="$(normalize_bool "${2:?Missing value}")"; shift 2 ;;
        --unique-to-cluster) UNIQUE_TO_CLUSTER="$(normalize_bool "${2:?Missing value}")"; shift 2 ;;
        --uniqueness-mode) UNIQUENESS_MODE="${2:?Missing value}"; shift 2 ;;
        --minimum-peptide-count) MINIMUM_PEPTIDE_COUNT="${2:?Missing value}"; shift 2 ;;
        --peptide-counting-mode) PEPTIDE_COUNTING_MODE="${2:?Missing value}"; shift 2 ;;
        --cluster-sums) CLUSTER_SUMS="$(normalize_bool "${2:?Missing value}")"; shift 2 ;;
        --rollup-type) ROLLUP_TYPE="${2:?Missing value}"; shift 2 ;;
        --rrollup-mode) RROLLUP_MODE="${2:?Missing value}"; shift 2 ;;
        --overwrite) OVERWRITE="True"; shift ;;
        --skip-derive) SKIP_DERIVE="True"; shift ;;
        --skip-apply) SKIP_APPLY="True"; shift ;;
        --dry-run) DRY_RUN="True"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -n "$DOWNSTREAM_DIR" ]] || die "--downstream-dir is required."
[[ -n "$FULL_DIR" ]] || die "--full-dir is required."

DOWNSTREAM_DIR="$(realpath -m "$DOWNSTREAM_DIR")"
FULL_DIR="$(realpath -m "$FULL_DIR")"
require_dir "$DOWNSTREAM_DIR"
require_dir "$FULL_DIR"
require_dir "$FULL_DIR/fdr_esti"

# Infer the Binning software root from this shell when possible; retain the
# current working directory as a compatibility fallback for development trees.
if [[ -z "$PROJECT_ROOT" ]]; then
    SHELL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -d "$SHELL_DIR/../correction" && -d "$SHELL_DIR/../downstream" ]]; then
        PROJECT_ROOT="$(cd "$SHELL_DIR/../.." && pwd -P)"
    elif [[ -d "$PWD/scripts/correction" && -d "$PWD/scripts/downstream" ]]; then
        PROJECT_ROOT="$PWD"
    else
        die "Could not infer Binning software root. Supply --project-root."
    fi
fi
PROJECT_ROOT="$(realpath -m "$PROJECT_ROOT")"

DERIVE_SCRIPT="$PROJECT_ROOT/scripts/correction/derive_fdr_correction.py"
APPLY_SCRIPT="$PROJECT_ROOT/scripts/correction/apply_fdr_correction.py"
DOWNSTREAM_TOOLS_DIR="$PROJECT_ROOT/scripts/downstream"
require_file "$DERIVE_SCRIPT"
require_file "$APPLY_SCRIPT"
require_dir "$DOWNSTREAM_TOOLS_DIR"

# Infer the shared comparison table. This works for both names such as
# classic_vs_group_fdr_audit and classic_vs_binned_fdr_audit.
if [[ -z "$SHARED_TSV" ]]; then
    mapfile -t SHARED_CANDIDATES < <(
        find "$DOWNSTREAM_DIR" -maxdepth 3 -type f \
            -path '*/classic_vs_*_fdr_audit/classic_vs_binned_scan_peptide_fdr_wide.tsv' \
            -print | sort
    )
    [[ ${#SHARED_CANDIDATES[@]} -eq 1 ]] || {
        printf 'Shared TSV candidates found: %s\n' "${#SHARED_CANDIDATES[@]}" >&2
        printf '  %s\n' "${SHARED_CANDIDATES[@]:-<none>}" >&2
        die "Expected exactly one shared audit TSV. Supply --shared-tsv explicitly."
    }
    SHARED_TSV="${SHARED_CANDIDATES[0]}"
else
    SHARED_TSV="$(realpath -m "$SHARED_TSV")"
fi
require_file "$SHARED_TSV"

if [[ -z "$CONDITION_NAME" ]]; then
    CONDITION_NAME="$(basename "$(dirname "$DOWNSTREAM_DIR")")"
fi
[[ "$CONDITION_NAME" != *'/'* ]] || die "--condition-name must not contain '/'."

FORMULA_DIR="${FORMULA_DIR:-$DOWNSTREAM_DIR/correctional_formula/$CONDITION_NAME}"
CORRECTED_DIR="${CORRECTED_DIR:-$DOWNSTREAM_DIR/corrected_fdr_esti}"
ROLLUP_DIR="${ROLLUP_DIR:-$DOWNSTREAM_DIR/protein_rollup}"
FORMULA_DIR="$(realpath -m "$FORMULA_DIR")"
CORRECTED_DIR="$(realpath -m "$CORRECTED_DIR")"
ROLLUP_DIR="$(realpath -m "$ROLLUP_DIR")"
FORMULA_JSON="$FORMULA_DIR/correction_formula.json"

# Locate the single top-level FASTA/FAA from the classical full run. Using it
# by absolute path avoids copying or modifying a large source database.
mapfile -t DB_FILES < <(
    find "$FULL_DIR" -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.faa' \) \
        ! -name '._*' -print | sort
)
[[ ${#DB_FILES[@]} -eq 1 ]] || {
    printf 'Top-level FASTA/FAA candidates found: %s\n' "${#DB_FILES[@]}" >&2
    printf '  %s\n' "${DB_FILES[@]:-<none>}" >&2
    die "Expected exactly one top-level FASTA/FAA in --full-dir."
}
DB_FILE="${DB_FILES[0]}"

# Validate enumerated settings early.
case "$UNIQUENESS_MODE" in all_matches|unique_only|requires_unique) ;; *) die "Invalid --uniqueness-mode: $UNIQUENESS_MODE" ;; esac
case "$PEPTIDE_COUNTING_MODE" in per_sample|global) ;; *) die "Invalid --peptide-counting-mode: $PEPTIDE_COUNTING_MODE" ;; esac

LOG_DIR="$DOWNSTREAM_DIR/logs/correctional_formula_pipeline"
if [[ "$DRY_RUN" == "False" ]]; then
    mkdir -p "$LOG_DIR"
    LOG_FILE="$LOG_DIR/run_$(date '+%Y%m%d_%H%M%S').log"
    exec > >(tee -a "$LOG_FILE") 2>&1
else
    LOG_FILE="<dry-run>"
fi

log "Correctional downstream pipeline v$SCRIPT_VERSION"
log "Project root:      $PROJECT_ROOT"
log "Downstream dir:    $DOWNSTREAM_DIR"
log "Classical full:    $FULL_DIR"
log "Shared TSV:        $SHARED_TSV"
log "Condition:         $CONDITION_NAME"
log "Source FASTA:      $DB_FILE"
log "Formula dir:       $FORMULA_DIR"
log "Corrected dir:     $CORRECTED_DIR"
log "Protein rollup:    $ROLLUP_DIR"
log "Log:               $LOG_FILE"
log "Parameters:        $PEPTIDE_FILTER, $STRINGENCY, $PEPTIDE_FILTER2, $STRINGENCY2, $HIGHEST_INTENSITY_ONLY, $UNIQUE_TO_CLUSTER, $UNIQUENESS_MODE, $MINIMUM_PEPTIDE_COUNT, $PEPTIDE_COUNTING_MODE, $CLUSTER_SUMS, $ROLLUP_TYPE, $RROLLUP_MODE"

# Managed outputs are protected unless --overwrite is explicit. Formula and
# corrected directories are independently reusable with --skip-*.
if [[ "$OVERWRITE" == "True" ]]; then
    if [[ "$SKIP_DERIVE" == "False" ]]; then
        run_cmd rm -rf -- "$FORMULA_DIR"
    fi
    if [[ "$SKIP_APPLY" == "False" ]]; then
        run_cmd rm -rf -- "$CORRECTED_DIR"
    fi
    run_cmd rm -rf -- "$ROLLUP_DIR"
else
    if [[ "$SKIP_DERIVE" == "False" && -e "$FORMULA_DIR" ]]; then
        die "Formula directory already exists: $FORMULA_DIR (use --overwrite or --skip-derive)"
    fi
    if [[ "$SKIP_APPLY" == "False" && -e "$CORRECTED_DIR" ]]; then
        die "Corrected directory already exists: $CORRECTED_DIR (use --overwrite or --skip-apply)"
    fi
    [[ ! -e "$ROLLUP_DIR" ]] || die "Protein rollup directory already exists: $ROLLUP_DIR (use --overwrite)"
fi

if [[ "$SKIP_DERIVE" == "False" ]]; then
    log "Step 1/8: Deriving the FDR correction formula"
    run_cmd python3 "$DERIVE_SCRIPT" \
        --shared "$SHARED_TSV" \
        --outdir "$FORMULA_DIR" \
        --condition_name "$CONDITION_NAME"
else
    log "Step 1/8: Reusing existing formula"
fi
if [[ "$DRY_RUN" == "False" ]]; then
    require_file "$FORMULA_JSON"
fi

if [[ "$SKIP_APPLY" == "False" ]]; then
    log "Step 2/8: Applying the FDR correction to classical *_withsyn.tsv files"
    run_cmd python3 "$APPLY_SCRIPT" \
        --formula "$FORMULA_JSON" \
        --full_dir "$FULL_DIR" \
        --outdir "$CORRECTED_DIR" \
        --thresholds "$THRESHOLDS"
else
    log "Step 2/8: Reusing existing corrected FDR tables"
fi
if [[ "$DRY_RUN" == "False" ]]; then
    require_dir "$CORRECTED_DIR"
    CORRECTED_FILES=("$CORRECTED_DIR"/*_withsyn.tsv)
    [[ ${#CORRECTED_FILES[@]} -gt 0 ]] || die "No corrected *_withsyn.tsv files found in $CORRECTED_DIR"
fi

run_cmd mkdir -p "$ROLLUP_DIR"

# The following tools write fixed output names. Execute them inside the isolated
# rollup folder while passing source resources by absolute path.
if [[ "$DRY_RUN" == "False" ]]; then
    pushd "$ROLLUP_DIR" >/dev/null
else
    log "DRY RUN working directory for Steps 3-8: $ROLLUP_DIR"
fi

log "Step 3/8: Filtering corrected peptide-spectrum matches"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/filter_peptides.py" \
    -i "$CORRECTED_DIR" \
    -o fdr_filtered \
    -f "$PEPTIDE_FILTER" \
    -t "$STRINGENCY" \
    --filter2 "$PEPTIDE_FILTER2" \
    --threshold2 "$STRINGENCY2"

log "Step 4/8: Generating the peptide crosstab"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/generate_peptide_crosstab.py" \
    -i fdr_filtered \
    -o peptide_crosstab.tsv \
    --only-highest-intensity "$HIGHEST_INTENSITY_ONLY"

log "Step 5/8: Annotating the peptide crosstab"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/annotate_peptide_crosstab.py" \
    -i peptide_crosstab.tsv \
    -f "$DB_FILE" \
    -o peptide_crosstab_annotated.tsv \
    -unique_to_cluster "$UNIQUE_TO_CLUSTER"

log "Step 6/8: Building the peptide-protein map"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/build_peptide_protein_map.py" \
    -i peptide_crosstab_annotated.tsv \
    -db "$FULL_DIR" \
    -p "$MINIMUM_PEPTIDE_COUNT" \
    --mode "$UNIQUENESS_MODE" \
    -f "$PEPTIDE_COUNTING_MODE"

log "Step 7/8: Generating the protein rollup"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/protein_rollup.py" \
    -i peptide_crosstab_annotated.tsv \
    -o protein_crosstab_annotated.tsv \
    --mode "$UNIQUENESS_MODE" \
    --rollup "$ROLLUP_TYPE" \
    --rrollup "$RROLLUP_MODE"

log "Step 8/8: Generating the cluster rollup when requested"
if [[ "$CLUSTER_SUMS" == "True" ]]; then
    run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/cluster_rollup.py" \
        -i peptide_crosstab_annotated.tsv \
        -o clustered_crosstab_annotated.tsv \
        --mode "$UNIQUENESS_MODE" \
        --rollup "$ROLLUP_TYPE" \
        --rrollup "$RROLLUP_MODE"
else
    log "Cluster rollup disabled."
fi

# Record the exact resolved configuration beside the outputs.
if [[ "$DRY_RUN" == "False" ]]; then
    cat > correctional_rollup_run_metadata.tsv <<META
metric\tvalue
shell_version\t$SCRIPT_VERSION
generated_local\t$(date --iso-8601=seconds)
project_root\t$PROJECT_ROOT
downstream_dir\t$DOWNSTREAM_DIR
full_dir\t$FULL_DIR
shared_tsv\t$SHARED_TSV
condition_name\t$CONDITION_NAME
formula_json\t$FORMULA_JSON
corrected_fdr_dir\t$CORRECTED_DIR
source_fasta\t$DB_FILE
peptide_filter\t$PEPTIDE_FILTER
stringency\t$STRINGENCY
peptide_filter2\t$PEPTIDE_FILTER2
stringency2\t$STRINGENCY2
highest_intensity_only\t$HIGHEST_INTENSITY_ONLY
unique_to_cluster\t$UNIQUE_TO_CLUSTER
uniqueness_mode\t$UNIQUENESS_MODE
minimum_peptide_count\t$MINIMUM_PEPTIDE_COUNT
peptide_counting_mode\t$PEPTIDE_COUNTING_MODE
cluster_sums\t$CLUSTER_SUMS
rollup_type\t$ROLLUP_TYPE
rrollup_mode\t$RROLLUP_MODE
META
fi

if [[ "$DRY_RUN" == "False" ]]; then
    popd >/dev/null
fi

log "All correctional-formula and downstream rollup steps completed successfully."
log "Formula package:    $FORMULA_DIR"
log "Corrected PSMs:     $CORRECTED_DIR"
log "Rollup products:    $ROLLUP_DIR"
