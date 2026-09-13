#!/usr/bin/env bash
#
# Discover all condition families and available correction designs under a
# Binning data root, then generate Matched_Formula_Output by calling
# run_matched_correction_downstream.sh for each eligible analysis.

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_VERSION="1.1.1"

usage() {
    cat <<'USAGE'
Usage:
  run_all_matched_formula_outputs.sh [--root PATH] [options]

Data location:
  --root PATH            Binning data root containing <Condition>_<integer> and
                         matching <Condition>_Full directories.
                         Default: <Binning>/data.

Normally inferred:
  --project-root PATH    Binning software root. Default: inferred from this
                         script when installed in scripts/pipeline/.
  --output-root PATH     Default: <root>/Matched_Formula_Output

Discovery behavior:
  - Finds every top-level <Condition>_<integer> directory that has a matching
    <Condition>_Full directory.
  - Discovers main Individual and Pairwise correction outputs when present.
  - Discovers available <Condition>_Groups_<LABEL> correction outputs without
    hard-coding dataset names or the number of numbered series.
  - The child runner performs automatic sample/formula matching and writes the
    provenance manifest for each analysis.

Options:
  --include REGEX        Only process target labels matching this regex.
  --exclude REGEX        Skip target labels matching this regex.
  --overwrite            Rebuild existing matched outputs.
  --fail-fast            Stop after the first failed analysis.
  --dry-run              Print discovered commands without running them.
  -h, --help             Show this help.

Examples:
  # Use the repository-local Binning/data directory.
  bash run_all_matched_formula_outputs.sh --dry-run

  # Or point to an external/HPC data root.
  bash run_all_matched_formula_outputs.sh \
      --root /path/to/binning_data \
      --include '^Fecal_' \
      --dry-run
USAGE
}

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

normalize_label() {
    local value="$1"
    value="${value// /_}"
    value="${value//\//_}"
    printf '%s' "$value" | sed -E 's/[^A-Za-z0-9._-]+/_/g; s/^_+//; s/_+$//'
}

ROOT=""
PROJECT_ROOT=""
OUTPUT_ROOT=""
INCLUDE_REGEX=""
EXCLUDE_REGEX=""
OVERWRITE="False"
FAIL_FAST="False"
DRY_RUN="False"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --root) ROOT="${2:?Missing value for --root}"; shift 2 ;;
        --project-root) PROJECT_ROOT="${2:?Missing value for --project-root}"; shift 2 ;;
        --output-root) OUTPUT_ROOT="${2:?Missing value for --output-root}"; shift 2 ;;
        --include) INCLUDE_REGEX="${2:?Missing value for --include}"; shift 2 ;;
        --exclude) EXCLUDE_REGEX="${2:?Missing value for --exclude}"; shift 2 ;;
        --overwrite) OVERWRITE="True"; shift ;;
        --fail-fast) FAIL_FAST="True"; shift ;;
        --dry-run) DRY_RUN="True"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -z "$PROJECT_ROOT" ]]; then
    candidate="$(realpath -m "$SCRIPT_DIR/../..")"
    if [[ -d "$candidate/scripts/pipeline" ]]; then
        PROJECT_ROOT="$candidate"
    else
        die "Could not infer Binning software root. Supply --project-root."
    fi
fi
PROJECT_ROOT="$(realpath -m "$PROJECT_ROOT")"

ROOT="${ROOT:-$PROJECT_ROOT/data}"
ROOT="$(realpath -m "$ROOT")"
[[ -d "$ROOT" ]] || die "Binning data root not found: $ROOT"
CHILD="$PROJECT_ROOT/scripts/pipeline/run_matched_correction_downstream.sh"
[[ -f "$CHILD" ]] || die "Matched correction child runner not found: $CHILD"

OUTPUT_ROOT="${OUTPUT_ROOT:-$ROOT/Matched_Formula_Output}"
OUTPUT_ROOT="$(realpath -m "$OUTPUT_ROOT")"

# Discover condition families from numbered top-level directories.
declare -A CONDITIONS_SEEN=()
for path in "$ROOT"/*; do
    [[ -d "$path" ]] || continue
    base="$(basename "$path")"
    if [[ "$base" =~ ^(.+)_([1-9][0-9]*)$ ]]; then
        condition="${BASH_REMATCH[1]}"
        [[ -d "$ROOT/${condition}_Full" ]] || continue
        CONDITIONS_SEEN["$condition"]=1
    fi
done

mapfile -t CONDITIONS < <(printf '%s\n' "${!CONDITIONS_SEEN[@]}" | sort -V)
((${#CONDITIONS[@]} > 0)) || die "No numbered condition families with matching *_Full directories were discovered."

# Target arrays.
declare -a TARGET_CONDITION=()
declare -a TARGET_DESIGN=()
declare -a TARGET_SIZE=()
declare -a TARGET_LABEL=()
declare -a TARGET_OUTPUT=()

add_target() {
    local condition="$1"
    local design="$2"
    local size="$3"
    local label output

    case "$design" in
        individual) label="${condition}_Individual"; output="$OUTPUT_ROOT/$condition/Individual" ;;
        pairwise) label="${condition}_Pairwise"; output="$OUTPUT_ROOT/$condition/Pairwise" ;;
        size) label="${condition}_${size}"; output="$OUTPUT_ROOT/$condition/$size" ;;
        *) die "Internal error: unknown design $design" ;;
    esac
    label="$(normalize_label "$label")"

    if [[ -n "$INCLUDE_REGEX" ]] && ! [[ "$label" =~ $INCLUDE_REGEX ]]; then
        return 0
    fi
    if [[ -n "$EXCLUDE_REGEX" ]] && [[ "$label" =~ $EXCLUDE_REGEX ]]; then
        return 0
    fi

    TARGET_CONDITION+=("$condition")
    TARGET_DESIGN+=("$design")
    TARGET_SIZE+=("$size")
    TARGET_LABEL+=("$label")
    TARGET_OUTPUT+=("$output")
}

for condition in "${CONDITIONS[@]}"; do
    # Use the first numbered series only for design discovery; the child runner
    # validates that every numbered series has the corresponding analysis.
    declare -a series_candidates=()
    for path in "$ROOT"/"${condition}_"*; do
        [[ -d "$path" ]] || continue
        base="$(basename "$path")"
        suffix="${base#${condition}_}"
        [[ "$suffix" =~ ^[1-9][0-9]*$ ]] || continue
        series_candidates+=("$suffix"$'\t'"$path")
    done
    mapfile -t series_candidates < <(printf '%s\n' "${series_candidates[@]}" | sort -t $'\t' -k1,1n)
    ((${#series_candidates[@]} > 0)) || continue
    IFS=$'\t' read -r _ first_series <<< "${series_candidates[0]}"

    main_root="$first_series/${condition}_Binning"
    if [[ -d "$main_root/Group_Downstream/correctional_formula" ]] && \
       find "$main_root/Group_Downstream/correctional_formula" -type f -name correction_formula.json -print -quit | grep -q .; then
        add_target "$condition" individual ""
    fi
    if [[ -d "$main_root/Pairwise_Downstream/correctional_formula" ]] && \
       find "$main_root/Pairwise_Downstream/correctional_formula" -type f -name correction_formula.json -print -quit | grep -q .; then
        add_target "$condition" pairwise ""
    fi

    # Discover all available size labels generically.
    declare -a discovered_sizes=()
    for size_root in "$first_series"/"${condition}_Groups_"*; do
        [[ -d "$size_root" ]] || continue
        size_base="$(basename "$size_root")"
        size="${size_base#${condition}_Groups_}"
        [[ -n "$size" ]] || continue
        formula_root="$size_root/Group_Downstream/correctional_formula"
        [[ -d "$formula_root" ]] || continue
        find "$formula_root" -type f -name correction_formula.json -print -quit | grep -q . || continue
        discovered_sizes+=("$size")
    done
    mapfile -t discovered_sizes < <(printf '%s\n' "${discovered_sizes[@]:-}" | sed '/^$/d' | sort -u)

    # Familiar manuscript labels are ordered first when present; all additional labels are discovered generically.
    declare -A emitted_size=()
    for size in "1.2M" "600k" "300k" "150k" "75k"; do
        for discovered in "${discovered_sizes[@]}"; do
            if [[ "$discovered" == "$size" ]]; then
                add_target "$condition" size "$size"
                emitted_size["$size"]=1
                break
            fi
        done
    done
    for size in "${discovered_sizes[@]}"; do
        [[ -n "${emitted_size[$size]:-}" ]] && continue
        add_target "$condition" size "$size"
    done

done

((${#TARGET_LABEL[@]} > 0)) || die "No eligible matched-correction analyses were discovered after filtering."

STAMP="$(date '+%Y%m%d_%H%M%S')"
BATCH_LOG_DIR="$OUTPUT_ROOT/_batch_logs"
MASTER_LOG="$BATCH_LOG_DIR/matched_mass_run_${STAMP}.txt"
SUMMARY_TSV="$BATCH_LOG_DIR/matched_mass_run_${STAMP}_summary.tsv"

if [[ "$DRY_RUN" == "False" ]]; then
    mkdir -p "$BATCH_LOG_DIR"
    exec > >(tee -a "$MASTER_LOG") 2>&1
    printf '%s\n' $'Target\tCondition\tDesign\tBinSize\tStatus\tExitCode\tElapsedSeconds\tOutputDirectory\tTranscript' > "$SUMMARY_TSV"
fi

log "Matched formula mass runner v$SCRIPT_VERSION"
log "Root:          $ROOT"
log "Repository:    $PROJECT_ROOT"
log "Output root:   $OUTPUT_ROOT"
log "Conditions:    ${CONDITIONS[*]}"
log "Targets:       ${#TARGET_LABEL[@]}"
log "Overwrite:     $OVERWRITE"
log "Fail fast:     $FAIL_FAST"
log "Dry run:       $DRY_RUN"

SUCCESS=0
FAILED=0

for i in "${!TARGET_LABEL[@]}"; do
    condition="${TARGET_CONDITION[$i]}"
    design="${TARGET_DESIGN[$i]}"
    size="${TARGET_SIZE[$i]}"
    label="${TARGET_LABEL[$i]}"
    output="${TARGET_OUTPUT[$i]}"
    transcript="$BATCH_LOG_DIR/${label}_${STAMP}.txt"

    cmd=(
        bash "$CHILD"
        --root "$ROOT"
        --project-root "$PROJECT_ROOT"
        --condition "$condition"
        --design "$design"
        --output-dir "$output"
    )
    [[ -n "$size" ]] && cmd+=(--bin-size "$size")
    [[ "$OVERWRITE" == "True" ]] && cmd+=(--overwrite)
    [[ "$DRY_RUN" == "True" ]] && cmd+=(--dry-run)

    log "======================================================================"
    log "START [$((i + 1))/${#TARGET_LABEL[@]}] $label"
    printf 'RUN:'
    printf ' %q' "${cmd[@]}"
    printf '\n'

    if [[ "$DRY_RUN" == "True" ]]; then
        "${cmd[@]}"
        continue
    fi

    start_epoch="$(date +%s)"
    set +e
    "${cmd[@]}" 2>&1 | tee "$transcript"
    status=${PIPESTATUS[0]}
    set -e
    elapsed="$(( $(date +%s) - start_epoch ))"

    if [[ $status -eq 0 ]]; then
        run_status="SUCCESS"
        ((SUCCESS+=1))
    else
        run_status="FAILED"
        ((FAILED+=1))
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$label" "$condition" "$design" "$size" "$run_status" "$status" \
        "$elapsed" "$output" "$transcript" >> "$SUMMARY_TSV"

    log "DONE [$label]: $run_status (exit=$status, ${elapsed}s)"
    if [[ $status -ne 0 && "$FAIL_FAST" == "True" ]]; then
        log "Stopping because --fail-fast was requested."
        exit "$status"
    fi
done

if [[ "$DRY_RUN" == "True" ]]; then
    log "Dry run complete; no matched outputs were created."
    exit 0
fi

log "======================================================================"
log "Matched formula mass run complete."
log "Successful targets: $SUCCESS"
log "Failed targets:     $FAILED"
log "Master log:         $MASTER_LOG"
log "Summary TSV:        $SUMMARY_TSV"

[[ $FAILED -eq 0 ]]
