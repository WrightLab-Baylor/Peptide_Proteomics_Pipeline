#!/usr/bin/env bash
#
# Apply each numbered-series correction formula only to the sample from which
# that formula was derived, combine the matched corrected PSM tables, and run
# the standard peptide/protein downstream workflow once on the matched set.
#
# Expected Binning data layout:
#   <root>/<Condition>_1/
#   <root>/<Condition>_2/
#   ...
#   <root>/<Condition>_Full/
#
# Sample assignment is automatic. Each numbered series is matched to its full-
# database *_withsyn.tsv by the unique RAW/mzML stem under shared_data. The
# resolved mapping is written to matched_formula_manifest.tsv as provenance.

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_VERSION="1.1.1"

usage() {
    cat <<'USAGE'
Usage:
  run_matched_correction_downstream.sh \
    [--root PATH] \
    --condition NAME \
    --design individual|pairwise|size \
    [--bin-size LABEL] \
    [options]

Purpose:
  Build a matched correction output in which each numbered-series formula is
  applied only to its corresponding source sample, then run the standard
  peptide/protein downstream workflow on the combined corrected files.

Required:
  --condition NAME     Condition family, e.g. Fecal, Ocean, Soil.
  --design TYPE        individual, pairwise, or size.

Data location:
  --root PATH          Binning data root containing <Condition>_<integer> and
                       <Condition>_Full directories. Default: <Binning>/data.

Design-specific:
  --bin-size LABEL     Required only for --design size, e.g. 1.2M or 150k.

Normally inferred:
  --project-root PATH  Binning software root. Default: inferred from this
                       script when installed in scripts/pipeline/.
  --output-dir PATH    Analysis output directory. Default:
                       <root>/Matched_Formula_Output/<Condition>/<Analysis>

Automatic matching safeguards:
  1. Discover exactly one RAW/mzML sample stem per numbered series from
     shared_data (recursive fallback permitted).
  2. Optionally cross-check the unique PHRP *_syn.txt basename when available.
  3. Resolve exactly one matching <Condition>_Full/fdr_esti/*_withsyn.tsv.
  4. Resolve exactly one correction_formula.json for that series/design.
  5. Refuse ambiguous sample, reference, or formula matches before processing.
  6. Record every resolved path in matched_formula_manifest.tsv.

Output structure:
  <output-dir>/
    formula_applications/Series_<N>/
    corrected_fdr_esti/
    protein_rollup/
    logs/
    matched_formula_manifest.tsv
    matched_correction_application_summary.tsv
    matched_formula_run_metadata.tsv

Validated downstream defaults:
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
  --thresholds 0.01,0.05

Other options:
  --overwrite          Rebuild managed outputs in <output-dir>.
  --dry-run            Resolve and print planned matches/commands only.
  -h, --help           Show this help.
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

normalize_bool() {
    case "${1,,}" in
        true|t|1|yes|y) printf 'True\n' ;;
        false|f|0|no|n) printf 'False\n' ;;
        *) die "Expected a boolean value, received: $1" ;;
    esac
}

normalize_name() {
    printf '%s' "$1" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9]+//g'
}

require_dir() {
    [[ -d "$1" ]] || die "Required directory not found: $1"
}

require_file() {
    [[ -f "$1" ]] || die "Required file not found: $1"
}

# Canonical source-sample identity: unique RAW/mzML stem under shared_data.
discover_sample_stem() {
    local series_root="$1"
    local -a search_roots=()
    local direct="$series_root/shared_data"

    if [[ -d "$direct" ]]; then
        search_roots+=("$direct")
    else
        while IFS= read -r p; do
            search_roots+=("$p")
        done < <(find "$series_root" -type d -name shared_data -print 2>/dev/null | sort)
    fi

    ((${#search_roots[@]} > 0)) || {
        printf 'No shared_data directory found under %s\n' "$series_root" >&2
        return 1
    }

    local -a stems=()
    local root path base
    for root in "${search_roots[@]}"; do
        while IFS= read -r path; do
            base="$(basename "$path")"
            case "$base" in
                *.raw|*.RAW) stems+=("${base%.*}") ;;
                *.mzML|*.mzml) stems+=("${base%.*}") ;;
            esac
        done < <(
            find "$root" -maxdepth 1 \
                \( -iname '*.raw' -o -iname '*.mzml' \) -print 2>/dev/null | sort
        )
    done

    mapfile -t stems < <(printf '%s\n' "${stems[@]:-}" | sed '/^$/d' | sort -u)
    if [[ ${#stems[@]} -ne 1 ]]; then
        printf 'Expected exactly one unique RAW/mzML stem under %s; found %s\n' \
            "$series_root" "${#stems[@]}" >&2
        printf '  %s\n' "${stems[@]:-<none>}" >&2
        return 1
    fi
    printf '%s\n' "${stems[0]}"
}

# Optional independent sanity check from PHRP outputs. Zero matches are allowed
# because some cleaned/exported trees may omit partition-level PHRP files.
discover_syn_stem() {
    local analysis_root="$1"
    local -a stems=()
    local path base
    while IFS= read -r path; do
        base="$(basename "$path")"
        stems+=("${base%_syn.txt}")
    done < <(
        find "$analysis_root" -type f -path '*/results/PHRPOut/*_syn.txt' -print 2>/dev/null | sort
    )
    mapfile -t stems < <(printf '%s\n' "${stems[@]:-}" | sed '/^$/d' | sort -u)

    if [[ ${#stems[@]} -eq 0 ]]; then
        printf 'NA\n'
        return 0
    fi
    if [[ ${#stems[@]} -ne 1 ]]; then
        printf 'Expected zero or one unique SYN basename under %s; found %s\n' \
            "$analysis_root" "${#stems[@]}" >&2
        printf '  %s\n' "${stems[@]}" >&2
        return 1
    fi
    printf '%s\n' "${stems[0]}"
}

resolve_full_reference() {
    local full_dir="$1"
    local sample_stem="$2"
    local fdr_dir="$full_dir/fdr_esti"
    require_dir "$fdr_dir"

    local exact="$fdr_dir/${sample_stem}_withsyn.tsv"
    if [[ -f "$exact" ]]; then
        printf '%s\n' "$exact"
        return 0
    fi

    local target_norm
    target_norm="$(normalize_name "$sample_stem")"
    local -a exact_norm=()
    local -a loose_norm=()
    local path base stem norm
    while IFS= read -r path; do
        base="$(basename "$path")"
        stem="${base%_withsyn.tsv}"
        norm="$(normalize_name "$stem")"
        if [[ "$norm" == "$target_norm" ]]; then
            exact_norm+=("$path")
        elif [[ "$norm" == *"$target_norm"* || "$target_norm" == *"$norm"* ]]; then
            loose_norm+=("$path")
        fi
    done < <(find "$fdr_dir" -maxdepth 1 -type f -name '*_withsyn.tsv' -print | sort)

    if [[ ${#exact_norm[@]} -eq 1 ]]; then
        printf '%s\n' "${exact_norm[0]}"
        return 0
    fi
    if [[ ${#exact_norm[@]} -eq 0 && ${#loose_norm[@]} -eq 1 ]]; then
        printf '%s\n' "${loose_norm[0]}"
        return 0
    fi

    printf 'Could not uniquely resolve sample stem %q in %s\n' "$sample_stem" "$fdr_dir" >&2
    printf '  exact-normalized matches: %s\n' "${#exact_norm[@]}" >&2
    printf '  loose-normalized matches: %s\n' "${#loose_norm[@]}" >&2
    printf '  %s\n' "${exact_norm[@]:-}" "${loose_norm[@]:-}" >&2
    return 1
}

resolve_formula_json() {
    local downstream_dir="$1"
    local formula_root="$downstream_dir/correctional_formula"
    [[ -d "$formula_root" ]] || {
        printf 'Correctional formula directory not found: %s\n' "$formula_root" >&2
        return 1
    }

    local -a formulas=()
    while IFS= read -r p; do formulas+=("$p"); done < <(
        find "$formula_root" -type f -name correction_formula.json -print | sort
    )
    if [[ ${#formulas[@]} -ne 1 ]]; then
        printf 'Expected exactly one correction_formula.json under %s; found %s\n' \
            "$formula_root" "${#formulas[@]}" >&2
        printf '  %s\n' "${formulas[@]:-<none>}" >&2
        return 1
    fi
    printf '%s\n' "${formulas[0]}"
}

ROOT=""
CONDITION=""
DESIGN=""
BIN_SIZE=""
PROJECT_ROOT=""
OUTPUT_DIR=""

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
DRY_RUN="False"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root) ROOT="${2:?Missing value for --root}"; shift 2 ;;
        --condition) CONDITION="${2:?Missing value for --condition}"; shift 2 ;;
        --design) DESIGN="${2:?Missing value for --design}"; shift 2 ;;
        --bin-size) BIN_SIZE="${2:?Missing value for --bin-size}"; shift 2 ;;
        --project-root) PROJECT_ROOT="${2:?Missing value for --project-root}"; shift 2 ;;
        --output-dir) OUTPUT_DIR="${2:?Missing value for --output-dir}"; shift 2 ;;
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
        --dry-run) DRY_RUN="True"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -n "$CONDITION" ]] || die "--condition is required."
[[ "$DESIGN" =~ ^(individual|pairwise|size)$ ]] || die "--design must be individual, pairwise, or size."
if [[ "$DESIGN" == "size" ]]; then
    [[ -n "$BIN_SIZE" ]] || die "--bin-size is required with --design size."
else
    [[ -z "$BIN_SIZE" ]] || die "--bin-size is only valid with --design size."
fi

case "$UNIQUENESS_MODE" in all_matches|unique_only|requires_unique) ;; *) die "Invalid --uniqueness-mode: $UNIQUENESS_MODE" ;; esac
case "$PEPTIDE_COUNTING_MODE" in per_sample|global) ;; *) die "Invalid --peptide-counting-mode: $PEPTIDE_COUNTING_MODE" ;; esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -z "$PROJECT_ROOT" ]]; then
    candidate="$(realpath -m "$SCRIPT_DIR/../..")"
    if [[ -f "$candidate/scripts/correction/apply_fdr_correction.py" && -d "$candidate/scripts/downstream" ]]; then
        PROJECT_ROOT="$candidate"
    else
        die "Could not infer Binning software root. Supply --project-root."
    fi
fi
PROJECT_ROOT="$(realpath -m "$PROJECT_ROOT")"

ROOT="${ROOT:-$PROJECT_ROOT/data}"
ROOT="$(realpath -m "$ROOT")"
require_dir "$ROOT"

APPLY_SCRIPT="$PROJECT_ROOT/scripts/correction/apply_fdr_correction.py"
DOWNSTREAM_TOOLS_DIR="$PROJECT_ROOT/scripts/downstream"
require_file "$APPLY_SCRIPT"
require_dir "$DOWNSTREAM_TOOLS_DIR"

FULL_DIR="$ROOT/${CONDITION}_Full"
require_dir "$FULL_DIR"
require_dir "$FULL_DIR/fdr_esti"

case "$DESIGN" in
    individual)
        ANALYSIS_LABEL="Individual"
        ;;
    pairwise)
        ANALYSIS_LABEL="Pairwise"
        ;;
    size)
        ANALYSIS_LABEL="$BIN_SIZE"
        ;;
esac

OUTPUT_DIR="${OUTPUT_DIR:-$ROOT/Matched_Formula_Output/$CONDITION/$ANALYSIS_LABEL}"
OUTPUT_DIR="$(realpath -m "$OUTPUT_DIR")"
APPLICATION_ROOT="$OUTPUT_DIR/formula_applications"
CORRECTED_DIR="$OUTPUT_DIR/corrected_fdr_esti"
ROLLUP_DIR="$OUTPUT_DIR/protein_rollup"
LOG_DIR="$OUTPUT_DIR/logs"
MANIFEST="$OUTPUT_DIR/matched_formula_manifest.tsv"
APP_SUMMARY="$OUTPUT_DIR/matched_correction_application_summary.tsv"
RUN_METADATA="$OUTPUT_DIR/matched_formula_run_metadata.tsv"

# Discover numbered series generically and sort by numeric suffix.
declare -a SERIES_ROWS=()
for candidate in "$ROOT"/"${CONDITION}_"*; do
    [[ -d "$candidate" ]] || continue
    base="$(basename "$candidate")"
    suffix="${base#${CONDITION}_}"
    [[ "$suffix" =~ ^[1-9][0-9]*$ ]] || continue
    SERIES_ROWS+=("$suffix"$'\t'"$candidate")
done
((${#SERIES_ROWS[@]} > 0)) || die "No numbered series directories discovered for $CONDITION under $ROOT"
mapfile -t SERIES_ROWS < <(printf '%s\n' "${SERIES_ROWS[@]}" | sort -t $'\t' -k1,1n)

# Locate exactly one top-level FASTA/FAA for the classical full run.
mapfile -t DB_FILES < <(
    find "$FULL_DIR" -maxdepth 1 \( -type f -o -type l \) \
        \( -iname '*.fasta' -o -iname '*.faa' \) ! -name '._*' -print | sort
)
[[ ${#DB_FILES[@]} -eq 1 ]] || {
    printf 'Top-level FASTA/FAA candidates found: %s\n' "${#DB_FILES[@]}" >&2
    printf '  %s\n' "${DB_FILES[@]:-<none>}" >&2
    die "Expected exactly one top-level FASTA/FAA in $FULL_DIR"
}
DB_FILE="${DB_FILES[0]}"

# Resolve every match before creating scientific outputs. This makes matching a
# true preflight check rather than allowing a partial mixed-formula result.
declare -a RESOLVED_ROWS=()
declare -A SEEN_REFERENCE=()

log "Matched correction downstream v$SCRIPT_VERSION — preflight"
log "Root:          $ROOT"
log "Condition:     $CONDITION"
log "Design:        $DESIGN${BIN_SIZE:+ ($BIN_SIZE)}"
log "Series count:  ${#SERIES_ROWS[@]}"
log "Full run:      $FULL_DIR"
log "Output:        $OUTPUT_DIR"

for packed in "${SERIES_ROWS[@]}"; do
    IFS=$'\t' read -r series series_root <<< "$packed"

    case "$DESIGN" in
        individual)
            analysis_root="$series_root/${CONDITION}_Binning"
            downstream_dir="$analysis_root/Group_Downstream"
            ;;
        pairwise)
            analysis_root="$series_root/${CONDITION}_Binning"
            downstream_dir="$analysis_root/Pairwise_Downstream"
            ;;
        size)
            analysis_root="$series_root/${CONDITION}_Groups_${BIN_SIZE}"
            downstream_dir="$analysis_root/Group_Downstream"
            ;;
    esac

    require_dir "$analysis_root"
    require_dir "$downstream_dir"

    sample_stem="$(discover_sample_stem "$series_root")" || die "Sample discovery failed for Series $series"
    syn_stem="$(discover_syn_stem "$analysis_root")" || die "SYN cross-check failed for Series $series"
    if [[ "$syn_stem" != "NA" ]]; then
        sample_norm="$(normalize_name "$sample_stem")"
        syn_norm="$(normalize_name "$syn_stem")"
        [[ "$sample_norm" == "$syn_norm" ]] || die \
            "Series $series sample mismatch: shared_data=$sample_stem, PHRP=$syn_stem"
        syn_status="PASS"
    else
        syn_status="NOT_AVAILABLE"
    fi

    reference="$(resolve_full_reference "$FULL_DIR" "$sample_stem")" || die "Full-reference resolution failed for Series $series"
    formula_json="$(resolve_formula_json "$downstream_dir")" || die "Formula resolution failed for Series $series"

    ref_base="$(basename "$reference")"
    if [[ -n "${SEEN_REFERENCE[$ref_base]:-}" ]]; then
        die "Two numbered series resolved to the same full-run sample: $ref_base (Series ${SEEN_REFERENCE[$ref_base]} and $series)"
    fi
    SEEN_REFERENCE[$ref_base]="$series"

    RESOLVED_ROWS+=("$series"$'\t'"$series_root"$'\t'"$analysis_root"$'\t'"$downstream_dir"$'\t'"$sample_stem"$'\t'"$syn_stem"$'\t'"$syn_status"$'\t'"$reference"$'\t'"$formula_json")
    log "MATCH Series $series: $sample_stem -> $(basename "$reference") -> $formula_json"
done

# Protect existing managed outputs unless --overwrite is explicit.
managed=("$APPLICATION_ROOT" "$CORRECTED_DIR" "$ROLLUP_DIR" "$MANIFEST" "$APP_SUMMARY" "$RUN_METADATA")
if [[ "$DRY_RUN" == "False" ]]; then
    if [[ "$OVERWRITE" == "True" ]]; then
        rm -rf -- "$APPLICATION_ROOT" "$CORRECTED_DIR" "$ROLLUP_DIR"
        rm -f -- "$MANIFEST" "$APP_SUMMARY" "$RUN_METADATA"
    else
        for path in "${managed[@]}"; do
            [[ ! -e "$path" ]] || die "Managed output already exists: $path (use --overwrite)"
        done
    fi
    mkdir -p "$APPLICATION_ROOT" "$CORRECTED_DIR" "$ROLLUP_DIR" "$LOG_DIR"
    LOG_FILE="$LOG_DIR/run_$(date '+%Y%m%d_%H%M%S').log"
    exec > >(tee -a "$LOG_FILE") 2>&1
else
    LOG_FILE="<dry-run>"
fi

log "Preflight passed for ${#RESOLVED_ROWS[@]} matched series."

if [[ "$DRY_RUN" == "False" ]]; then
    printf '%s\n' $'Series\tSeriesDirectory\tAnalysisDirectory\tDownstreamDirectory\tAcquisitionStem\tPHRPSynStem\tPHRPCrossCheck\tFullReferenceWithSyn\tFormulaJSON\tCorrectedWithSyn\tMatchStatus' > "$MANIFEST"
fi

# Apply each formula to exactly one source file, retaining full per-sample
# application diagnostics in formula_applications/Series_<N>/.
for packed in "${RESOLVED_ROWS[@]}"; do
    IFS=$'\t' read -r series series_root analysis_root downstream_dir sample_stem syn_stem syn_status reference formula_json <<< "$packed"
    app_dir="$APPLICATION_ROOT/Series_${series}"
    corrected_file="$app_dir/$(basename "$reference")"
    combined_file="$CORRECTED_DIR/$(basename "$reference")"

    log "Applying Series $series formula only to $sample_stem"
    cmd=(
        python3 "$APPLY_SCRIPT"
        --formula "$formula_json"
        --full_dir "$FULL_DIR"
        --outdir "$app_dir"
        --pattern "$(basename "$reference")"
        --thresholds "$THRESHOLDS"
    )
    run_cmd "${cmd[@]}"

    if [[ "$DRY_RUN" == "False" ]]; then
        require_file "$corrected_file"
        [[ ! -e "$combined_file" ]] || die "Combined corrected filename collision: $combined_file"
        cp -p -- "$corrected_file" "$combined_file"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
            "$series" "$series_root" "$analysis_root" "$downstream_dir" \
            "$sample_stem" "$syn_stem" "$syn_status" "$reference" \
            "$formula_json" "$combined_file" >> "$MANIFEST"
    else
        log "DRY RUN would copy: $corrected_file -> $combined_file"
    fi
done

if [[ "$DRY_RUN" == "True" ]]; then
    log "Dry run complete; no files were created."
    exit 0
fi

# Aggregate per-sample application summaries without changing their originals.
python3 - "$APPLICATION_ROOT" "$MANIFEST" "$APP_SUMMARY" <<'PY'
from pathlib import Path
import sys
import pandas as pd

app_root = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
out_path = Path(sys.argv[3])
manifest = pd.read_csv(manifest_path, sep="\t")
frames = []
for row in manifest.itertuples(index=False):
    p = app_root / f"Series_{row.Series}" / "correction_application_summary.tsv"
    if not p.is_file():
        raise SystemExit(f"ERROR: missing application summary: {p}")
    df = pd.read_csv(p, sep="\t")
    df.insert(0, "Series", int(row.Series))
    df.insert(1, "MatchedSample", str(row.AcquisitionStem))
    df.insert(2, "FormulaJSON", str(row.FormulaJSON))
    frames.append(df)
pd.concat(frames, ignore_index=True).to_csv(out_path, sep="\t", index=False)
PY

# Standard post-correction downstream workflow. These commands intentionally
# mirror run_correctional_formula_downstream.sh so matched correction differs
# only in which formula is applied to each source *_withsyn.tsv.
pushd "$ROLLUP_DIR" >/dev/null

log "Downstream Step 1/6: Filtering corrected peptide-spectrum matches"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/filter_peptides.py" \
    -i "$CORRECTED_DIR" \
    -o fdr_filtered \
    -f "$PEPTIDE_FILTER" \
    -t "$STRINGENCY" \
    --filter2 "$PEPTIDE_FILTER2" \
    --threshold2 "$STRINGENCY2"

log "Downstream Step 2/6: Generating peptide crosstab"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/generate_peptide_crosstab.py" \
    -i fdr_filtered \
    -o peptide_crosstab.tsv \
    --only-highest-intensity "$HIGHEST_INTENSITY_ONLY"

log "Downstream Step 3/6: Annotating peptide crosstab"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/annotate_peptide_crosstab.py" \
    -i peptide_crosstab.tsv \
    -f "$DB_FILE" \
    -o peptide_crosstab_annotated.tsv \
    -unique_to_cluster "$UNIQUE_TO_CLUSTER"

log "Downstream Step 4/6: Building peptide-protein map"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/build_peptide_protein_map.py" \
    -i peptide_crosstab_annotated.tsv \
    -db "$FULL_DIR" \
    -p "$MINIMUM_PEPTIDE_COUNT" \
    --mode "$UNIQUENESS_MODE" \
    -f "$PEPTIDE_COUNTING_MODE"

log "Downstream Step 5/6: Generating protein rollup"
run_cmd python3 "$DOWNSTREAM_TOOLS_DIR/protein_rollup.py" \
    -i peptide_crosstab_annotated.tsv \
    -o protein_crosstab_annotated.tsv \
    --mode "$UNIQUENESS_MODE" \
    --rollup "$ROLLUP_TYPE" \
    --rrollup "$RROLLUP_MODE"

log "Downstream Step 6/6: Generating cluster rollup when requested"
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

cat > matched_formula_rollup_metadata.tsv <<META
metric\tvalue
shell_version\t$SCRIPT_VERSION
generated_local\t$(date --iso-8601=seconds)
root\t$ROOT
project_root\t$PROJECT_ROOT
condition\t$CONDITION
design\t$DESIGN
bin_size\t$BIN_SIZE
series_count\t${#RESOLVED_ROWS[@]}
full_dir\t$FULL_DIR
source_fasta\t$DB_FILE
corrected_fdr_dir\t$CORRECTED_DIR
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

popd >/dev/null

cat > "$RUN_METADATA" <<META
metric\tvalue
shell_version\t$SCRIPT_VERSION
generated_local\t$(date --iso-8601=seconds)
root\t$ROOT
project_root\t$PROJECT_ROOT
condition\t$CONDITION
design\t$DESIGN
analysis_label\t$ANALYSIS_LABEL
bin_size\t$BIN_SIZE
series_count\t${#RESOLVED_ROWS[@]}
full_dir\t$FULL_DIR
output_dir\t$OUTPUT_DIR
formula_applications\t$APPLICATION_ROOT
corrected_fdr_dir\t$CORRECTED_DIR
protein_rollup_dir\t$ROLLUP_DIR
match_manifest\t$MANIFEST
application_summary\t$APP_SUMMARY
log\t$LOG_FILE
META

log "Matched correction downstream completed successfully."
log "Match manifest:   $MANIFEST"
log "Corrected PSMs:   $CORRECTED_DIR"
log "Rollup products:  $ROLLUP_DIR"
log "Application QC:   $APP_SUMMARY"
