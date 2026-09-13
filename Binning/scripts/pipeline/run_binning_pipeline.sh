#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
cat <<'USAGE'
Usage:
  run_binning_pipeline.sh \
    --dataset-dir PATH \
    --reference-withsyn FILE \
    [options]

Purpose:
  Run the Binning partition workflow from one portable entry point while
  preserving the established scientific behavior and output structure.

Stages:
  1. Prepare SICs
  2. Classical-guided SIC culling
  3. Partition-level FDR estimation
  4. SYN protein-mapping merge
  5. Partition FDR collection and conflict audit
  6. Full-database vs binned FDR comparison
  7. FDR-delta interval analysis and plots
  8. Derive/apply the correction formula and generate corrected rollups

Required:
  --dataset-dir PATH
      Dataset directory containing *_pair_i_j and/or *_group_i partitions.

  --reference-withsyn FILE
      Full-database/classical *_withsyn.tsv reference.

Correction stage:
  --full-dir PATH
      Matching classical/full run directory. When supplied, Stage 8 is enabled
      by default. The directory must contain fdr_esti/*_withsyn.tsv and the
      top-level FASTA/FAA expected by the correctional downstream runner.

  --skip-correction
      Do not run Stage 8 even when --full-dir is supplied.

  --correction-only
      Run only Stage 8 using already-generated Stage 5-7 outputs. This avoids
      rerunning SIC preparation, culling, partition FDR estimation, or
      collection when only correctional processing needs to be regenerated.

  --overwrite-correction
      Permit the correction runner to rebuild its managed formula, corrected
      FDR, and rollup output directories.

  --skip-correction-derive
      Reuse an existing correction_formula.json.

  --skip-correction-apply
      Reuse existing corrected *_withsyn.tsv files.

Partition selection:
  --partition-type auto|pair|group|both   Default: auto
  --pair-glob GLOB                       Default: *_pair_*_*
  --group-glob GLOB                      Default: *_group_*
  --skip-pair
  --skip-group

Tool locations:
  --tools-dir PATH
      Compatibility override for development/testing. When omitted, Binning,
      FDR, and analysis tools are resolved from <Binning>/scripts/.

  --correction-runner FILE
      Default: run_correctional_formula_downstream.sh beside this launcher.

  --project-root PATH
      Binning software root. Default: inferred from this script location.

Stage control:
  --skip-prepare-sics
  --skip-cull-fdr-merge
  --skip-collect-compare-plot
  --skip-collector
  --skip-plots
  --audit-only
      Run SIC preparation if requested, then the culling audit only.
      Do not write culling/FDR/merge/downstream outputs.

Preparation controls:
  --clean-mac-hidden
  --force-prepare
  --sicdir-name NAME                     Default: SICdir
  --sics-name NAME                       Default: SICs

Culling/FDR/merge controls:
  --sic-dir-name NAME                    Default: SICs
  --sic-glob GLOB                        Default: *_PlusSICStats.tsv
  --syn-dir-name NAME                    Default: results/PHRPOut
  --syn-glob GLOB                        Default: *_syn.txt
  --pair-culled-dir-name NAME            Default: SICs_classical_guided_pair
  --group-culled-dir-name NAME           Default: SICs_classical_guided_group
  --pair-fdr-dir-name NAME               Default: fdr_esti_pair
  --group-fdr-dir-name NAME              Default: fdr_esti_group
  --overwrite-cull
  --force-fdr
  --force-merge
  --skip-merge
  --no-syn-fallback

Collection/comparison/plot controls:
  --pair-downstream-dir-name NAME        Default: Pairwise_Downstream
  --group-downstream-dir-name NAME       Default: Group_Downstream
  --pattern GLOB                         Default: *_withsyn.tsv
  --thresholds LIST                      Default: 0.01,0.05
  --include-unshared

General:
  --dry-run
      Preview commands and filesystem decisions without creating/modifying
      workflow outputs.

  -h, --help

Examples:
  # From the Binning software root, preview a run stored under Binning/data.
  bash scripts/pipeline/run_binning_pipeline.sh \
      --dataset-dir data/Soil_1/Soil_Binning \
      --reference-withsyn data/Soil_Full/fdr_esti/sample_withsyn.tsv \
      --dry-run

  # Run both detected partition types and enable correctional processing.
  bash scripts/pipeline/run_binning_pipeline.sh \
      --dataset-dir data/Soil_1/Soil_Binning \
      --reference-withsyn data/Soil_Full/fdr_esti/sample_withsyn.tsv \
      --full-dir data/Soil_Full

  # Re-run only formula derivation/application/rollup from existing downstream.
  bash scripts/pipeline/run_binning_pipeline.sh \
      --dataset-dir data/Soil_1/Soil_Binning \
      --reference-withsyn data/Soil_Full/fdr_esti/sample_withsyn.tsv \
      --full-dir data/Soil_Full \
      --correction-only \
      --overwrite-correction
USAGE
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

DATASET_DIR=""
REFERENCE=""
FULL_DIR=""
REPO_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd -P)"
BINNING_DIR="$REPO_ROOT/scripts/binning"
FDR_DIR="$REPO_ROOT/scripts/fdr"
ANALYSIS_DIR="$REPO_ROOT/scripts/analysis"
DOWNSTREAM_DIR="$REPO_ROOT/scripts/downstream"

# Compatibility overrides remain available for development/testing.
TOOLS_DIR=""
CORRECTION_RUNNER=""
PROJECT_ROOT="$REPO_ROOT"

PARTITION_TYPE="auto"
PAIR_GLOB="*_pair_*_*"
GROUP_GLOB="*_group_*"
SKIP_PAIR="False"
SKIP_GROUP="False"

SKIP_PREP="False"
SKIP_CULL="False"
SKIP_DOWN="False"
SKIP_COLLECTOR="False"
SKIP_PLOTS="False"
AUDIT_ONLY="False"

SKIP_CORRECTION="False"
CORRECTION_ONLY="False"
OVERWRITE_CORRECTION="False"
SKIP_CORRECTION_DERIVE="False"
SKIP_CORRECTION_APPLY="False"

CLEAN="False"
FORCE_PREP="False"
SICDIR_NAME="SICdir"
SICS_NAME="SICs"

CULL_SIC_DIR="SICs"
SIC_GLOB="*_PlusSICStats.tsv"
SYN_DIR="results/PHRPOut"
SYN_GLOB="*_syn.txt"
PAIR_CULLED="SICs_classical_guided_pair"
GROUP_CULLED="SICs_classical_guided_group"
PAIR_FDR="fdr_esti_pair"
GROUP_FDR="fdr_esti_group"
OVERWRITE_CULL="False"
FORCE_FDR="False"
FORCE_MERGE="False"
SKIP_MERGE="False"
NO_SYN="False"

PAIR_DOWN="Pairwise_Downstream"
GROUP_DOWN="Group_Downstream"
PATTERN="*_withsyn.tsv"
THRESHOLDS="0.01,0.05"
INCLUDE_UNSHARED="False"

DRY_RUN="False"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dataset-dir) DATASET_DIR="$2"; shift 2 ;;
        --reference-withsyn) REFERENCE="$2"; shift 2 ;;
        --full-dir) FULL_DIR="$2"; shift 2 ;;
        --tools-dir|--binning-tools-dir) TOOLS_DIR="$2"; shift 2 ;;
        --correction-runner) CORRECTION_RUNNER="$2"; shift 2 ;;
        --project-root) PROJECT_ROOT="$2"; shift 2 ;;

        --partition-type) PARTITION_TYPE="$2"; shift 2 ;;
        --pair-glob) PAIR_GLOB="$2"; shift 2 ;;
        --group-glob) GROUP_GLOB="$2"; shift 2 ;;
        --skip-pair) SKIP_PAIR="True"; shift ;;
        --skip-group) SKIP_GROUP="True"; shift ;;

        --skip-prepare-sics) SKIP_PREP="True"; shift ;;
        --skip-cull-fdr-merge) SKIP_CULL="True"; shift ;;
        --skip-collect-compare-plot) SKIP_DOWN="True"; shift ;;
        --skip-collector) SKIP_COLLECTOR="True"; shift ;;
        --skip-plots) SKIP_PLOTS="True"; shift ;;
        --audit-only) AUDIT_ONLY="True"; shift ;;
        --skip-correction) SKIP_CORRECTION="True"; shift ;;
        --correction-only) CORRECTION_ONLY="True"; shift ;;
        --overwrite-correction) OVERWRITE_CORRECTION="True"; shift ;;
        --skip-correction-derive) SKIP_CORRECTION_DERIVE="True"; shift ;;
        --skip-correction-apply) SKIP_CORRECTION_APPLY="True"; shift ;;

        --clean-mac-hidden) CLEAN="True"; shift ;;
        --force-prepare) FORCE_PREP="True"; shift ;;
        --sicdir-name) SICDIR_NAME="$2"; shift 2 ;;
        --sics-name) SICS_NAME="$2"; shift 2 ;;

        --sic-dir-name) CULL_SIC_DIR="$2"; shift 2 ;;
        --sic-glob) SIC_GLOB="$2"; shift 2 ;;
        --syn-dir-name) SYN_DIR="$2"; shift 2 ;;
        --syn-glob) SYN_GLOB="$2"; shift 2 ;;
        --pair-culled-dir-name) PAIR_CULLED="$2"; shift 2 ;;
        --group-culled-dir-name) GROUP_CULLED="$2"; shift 2 ;;
        --pair-fdr-dir-name) PAIR_FDR="$2"; shift 2 ;;
        --group-fdr-dir-name) GROUP_FDR="$2"; shift 2 ;;
        --overwrite-cull) OVERWRITE_CULL="True"; shift ;;
        --force-fdr) FORCE_FDR="True"; shift ;;
        --force-merge) FORCE_MERGE="True"; shift ;;
        --skip-merge) SKIP_MERGE="True"; shift ;;
        --no-syn-fallback) NO_SYN="True"; shift ;;

        --pair-downstream-dir-name) PAIR_DOWN="$2"; shift 2 ;;
        --group-downstream-dir-name) GROUP_DOWN="$2"; shift 2 ;;
        --pattern) PATTERN="$2"; shift 2 ;;
        --thresholds) THRESHOLDS="$2"; shift 2 ;;
        --include-unshared) INCLUDE_UNSHARED="True"; shift ;;

        --dry-run) DRY_RUN="True"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ "$PARTITION_TYPE" =~ ^(auto|pair|group|both)$ ]] || {
    echo "ERROR: invalid --partition-type: $PARTITION_TYPE" >&2
    exit 2
}

[[ -n "$DATASET_DIR" && -n "$REFERENCE" ]] || {
    echo "ERROR: --dataset-dir and --reference-withsyn are required" >&2
    exit 2
}

if [[ -n "$TOOLS_DIR" ]]; then
    BINNING_DIR="${TOOLS_DIR%/}"
    FDR_DIR="${TOOLS_DIR%/}"
    ANALYSIS_DIR="${TOOLS_DIR%/}"
fi
DATASET_DIR="${DATASET_DIR%/}"
FULL_DIR="${FULL_DIR%/}"

if [[ -z "$CORRECTION_RUNNER" ]]; then
    CORRECTION_RUNNER="$SCRIPT_DIR/run_correctional_formula_downstream.sh"
fi
CORRECTION_RUNNER="${CORRECTION_RUNNER%/}"

if [[ -z "$PROJECT_ROOT" ]]; then
    PROJECT_ROOT="$REPO_ROOT"
fi
PROJECT_ROOT="${PROJECT_ROOT%/}"

# Correction-only intentionally bypasses Stages 1-7 while retaining partition
# discovery so pair/group formula contexts can still be resolved.
if [[ "$CORRECTION_ONLY" == "True" ]]; then
    SKIP_PREP="True"
    SKIP_CULL="True"
    SKIP_DOWN="True"
    SKIP_CORRECTION="False"
fi

[[ -d "$DATASET_DIR" ]] || { echo "ERROR: missing dataset directory: $DATASET_DIR" >&2; exit 1; }
[[ -f "$REFERENCE" ]] || { echo "ERROR: missing reference withsyn: $REFERENCE" >&2; exit 1; }
[[ -d "$BINNING_DIR" ]] || { echo "ERROR: missing binning directory: $BINNING_DIR" >&2; exit 1; }
[[ -d "$FDR_DIR" ]] || { echo "ERROR: missing FDR directory: $FDR_DIR" >&2; exit 1; }
[[ -d "$ANALYSIS_DIR" ]] || { echo "ERROR: missing analysis directory: $ANALYSIS_DIR" >&2; exit 1; }
[[ -d "$DOWNSTREAM_DIR" ]] || { echo "ERROR: missing downstream directory: $DOWNSTREAM_DIR" >&2; exit 1; }

CORRECTION_ENABLED="False"
if [[ "$SKIP_CORRECTION" != "True" && -n "$FULL_DIR" ]]; then
    CORRECTION_ENABLED="True"
fi
if [[ "$CORRECTION_ONLY" == "True" && -z "$FULL_DIR" ]]; then
    echo "ERROR: --correction-only requires --full-dir" >&2
    exit 2
fi
if [[ "$CORRECTION_ENABLED" == "True" ]]; then
    [[ -d "$FULL_DIR" ]] || { echo "ERROR: missing full directory: $FULL_DIR" >&2; exit 1; }
    [[ -f "$CORRECTION_RUNNER" ]] || { echo "ERROR: missing correction runner: $CORRECTION_RUNNER" >&2; exit 1; }
    [[ -n "$PROJECT_ROOT" && -d "$PROJECT_ROOT" ]] || {
        echo "ERROR: could not infer --project-root for correction stage; provide it explicitly." >&2
        exit 1
    }
fi

PREPARE_SCRIPT="$BINNING_DIR/prepare_sics.py"
CULLER="$BINNING_DIR/classical_guided_cull.py"
FDR_SCRIPT="$FDR_DIR/estimate_fdr.py"
MERGER="$FDR_DIR/merge_protein_mappings.py"
COLLECTOR="$BINNING_DIR/collect_partition_fdr.py"
COMPARE="$ANALYSIS_DIR/compare_full_vs_binned_fdr.py"
PLOT="$ANALYSIS_DIR/plot_fdr_delta.py"

required_scripts=(
    "$PREPARE_SCRIPT"
    "$CULLER"
    "$FDR_SCRIPT"
    "$MERGER"
    "$COLLECTOR"
    "$COMPARE"
    "$PLOT"
)
for script in "${required_scripts[@]}"; do
    [[ -f "$script" ]] || { echo "ERROR: required script missing: $script" >&2; exit 1; }
done

# Discover and validate partition folders exactly once.
shopt -s nullglob
PAIR_CANDIDATES=("$DATASET_DIR"/$PAIR_GLOB)
GROUP_CANDIDATES=("$DATASET_DIR"/$GROUP_GLOB)
shopt -u nullglob

PAIRS=()
for folder in "${PAIR_CANDIDATES[@]}"; do
    [[ -d "$folder" ]] || continue
    name="$(basename "$folder")"
    [[ "$name" =~ _pair_[0-9]+_[0-9]+$ ]] && PAIRS+=("$folder")
done

GROUP_FOLDERS=()
for folder in "${GROUP_CANDIDATES[@]}"; do
    [[ -d "$folder" ]] || continue
    name="$(basename "$folder")"
    [[ "$name" =~ _group_[0-9]+$ ]] && GROUP_FOLDERS+=("$folder")
done

ACTIVE=()
case "$PARTITION_TYPE" in
    auto)
        [[ "$SKIP_PAIR" != "True" ]] && ((${#PAIRS[@]})) && ACTIVE+=(pair)
        [[ "$SKIP_GROUP" != "True" ]] && ((${#GROUP_FOLDERS[@]})) && ACTIVE+=(group)
        ;;
    pair)
        [[ "$SKIP_PAIR" != "True" ]] || { echo "ERROR: pair requested but --skip-pair was supplied" >&2; exit 2; }
        ((${#PAIRS[@]})) || { echo "ERROR: pair requested but no valid pair folders found" >&2; exit 1; }
        ACTIVE+=(pair)
        ;;
    group)
        [[ "$SKIP_GROUP" != "True" ]] || { echo "ERROR: group requested but --skip-group was supplied" >&2; exit 2; }
        ((${#GROUP_FOLDERS[@]})) || { echo "ERROR: group requested but no valid group folders found" >&2; exit 1; }
        ACTIVE+=(group)
        ;;
    both)
        [[ "$SKIP_PAIR" != "True" ]] || { echo "ERROR: both requested but --skip-pair was supplied" >&2; exit 2; }
        [[ "$SKIP_GROUP" != "True" ]] || { echo "ERROR: both requested but --skip-group was supplied" >&2; exit 2; }
        ((${#PAIRS[@]})) || { echo "ERROR: both requested but no valid pair folders found" >&2; exit 1; }
        ((${#GROUP_FOLDERS[@]})) || { echo "ERROR: both requested but no valid group folders found" >&2; exit 1; }
        ACTIVE+=(pair group)
        ;;
esac

((${#ACTIVE[@]})) || {
    echo "ERROR: all detected partition types were skipped or absent" >&2
    exit 1
}

printf 'Dataset:                  %s\n' "$DATASET_DIR"
printf 'Reference withsyn:        %s\n' "$REFERENCE"
printf 'Repository root:           %s\n' "$REPO_ROOT"
printf 'Detected pair folders:    %d\n' "${#PAIRS[@]}"
printf 'Detected group folders:   %d\n' "${#GROUP_FOLDERS[@]}"
printf 'Active partition types:   %s\n' "${ACTIVE[*]}"
printf 'Dry run:                  %s\n' "$DRY_RUN"
printf 'Audit only:               %s\n' "$AUDIT_ONLY"
printf 'Stages: prepare=%s cull/FDR/merge=%s downstream=%s\n' \
    "$([[ "$SKIP_PREP" == "True" ]] && echo skip || echo run)" \
    "$([[ "$SKIP_CULL" == "True" ]] && echo skip || echo run)" \
    "$([[ "$SKIP_DOWN" == "True" ]] && echo skip || echo run)"

run_or_preview() {
    if [[ "$DRY_RUN" == "True" ]]; then
        printf 'Would run:'
        printf ' %q' "$@"
        printf '\n'
    else
        "$@"
    fi
}

mode_folders() {
    local mode="$1"
    if [[ "$mode" == "pair" ]]; then
        printf '%s\n' "${PAIRS[@]}"
    else
        printf '%s\n' "${GROUP_FOLDERS[@]}"
    fi
}

# -----------------------------------------------------------------------------
# Stage 1: Prepare SICs
# -----------------------------------------------------------------------------
prepare_mode() {
    local mode="$1"
    local -a folders=()
    if [[ "$mode" == "pair" ]]; then folders=("${PAIRS[@]}"); else folders=("${GROUP_FOLDERS[@]}"); fi

    local log_dir="$DATASET_DIR/${mode}_prepare_sics_logs"
    local summary=""
    if [[ "$DRY_RUN" != "True" ]]; then
        mkdir -p "$log_dir"
        summary="$log_dir/prepare_sics_summary_$(date +%Y%m%d_%H%M%S).tsv"
        printf 'partition_folder\tstatus\tmessage\n' > "$summary"
    fi

    echo
    echo "=== Stage 1: prepare SICs [$mode] (${#folders[@]} partitions) ==="

    local folder name input_dir output_dir log_file
    for folder in "${folders[@]}"; do
        name="$(basename "$folder")"
        input_dir="$folder/$SICDIR_NAME"
        output_dir="$folder/$SICS_NAME"
        log_file="$log_dir/${name}.log"

        echo "[$mode] $name"

        if [[ ! -d "$input_dir" ]]; then
            echo "  WARNING: missing $input_dir"
            [[ -z "$summary" ]] || printf '%s\tskipped\tmissing input directory\n' "$name" >> "$summary"
            continue
        fi

        if [[ "$FORCE_PREP" != "True" && -d "$output_dir" ]] && \
           find "$output_dir" -type f -name '*.tsv' -print -quit | grep -q .; then
            echo "  existing prepared SICs; skipped"
            [[ -z "$summary" ]] || printf '%s\tskipped_existing\tprepared SICs already exist\n' "$name" >> "$summary"
            continue
        fi

        if [[ "$DRY_RUN" == "True" ]]; then
            [[ "$CLEAN" == "True" ]] && echo "  would remove ._* sidecar files"
            printf '  would run: python3 %q -i %q -o %q\n' "$PREPARE_SCRIPT" "$input_dir" "$output_dir"
            continue
        fi

        : > "$log_file"
        if [[ "$CLEAN" == "True" ]]; then
            find "$folder" -type f -name '._*' -print -delete | tee -a "$log_file"
        fi
        python3 "$PREPARE_SCRIPT" -i "$input_dir" -o "$output_dir" 2>&1 | tee -a "$log_file"
        printf '%s\tprocessed\tSICs prepared\n' "$name" >> "$summary"
    done

    [[ -z "$summary" ]] || echo "Summary: $summary"
}

# -----------------------------------------------------------------------------
# Stages 2-4: Cull -> FDR -> merge
# -----------------------------------------------------------------------------
cull_fdr_merge_mode() {
    local mode="$1"
    local -a folders=()
    local culled fdr report_dir

    if [[ "$mode" == "pair" ]]; then
        folders=("${PAIRS[@]}")
        culled="$PAIR_CULLED"
        fdr="$PAIR_FDR"
        report_dir="$DATASET_DIR/Pairwise_Downstream/classical_guided_culling"
    else
        folders=("${GROUP_FOLDERS[@]}")
        culled="$GROUP_CULLED"
        fdr="$GROUP_FDR"
        report_dir="$DATASET_DIR/Group_Downstream/classical_guided_culling"
    fi

    echo
    echo "=== Stages 2-4: cull/FDR/merge [$mode] ==="
    echo "Partitions: ${#folders[@]}"
    echo "Culled dir: $culled"
    echo "FDR dir:    $fdr"

    local -a cull_cmd=(
        python3 "$CULLER"
        --partition-type "$mode"
        --dataset-dir "$DATASET_DIR"
        --reference-withsyn "$REFERENCE"
        --output-prefix "$mode"
        --report-dir "$report_dir"
        --sic-dir-name "$CULL_SIC_DIR"
        --sic-glob "$SIC_GLOB"
        --syn-dir-name "$SYN_DIR"
        --syn-glob "$SYN_GLOB"
        --out-dir-name "$culled"
    )
    [[ "$NO_SYN" == "True" ]] && cull_cmd+=(--no-syn-fallback)
    [[ "$OVERWRITE_CULL" == "True" ]] && cull_cmd+=(--overwrite)
    [[ "$AUDIT_ONLY" == "True" ]] && cull_cmd+=(--dry-run)

    if [[ "$DRY_RUN" == "True" ]]; then
        printf 'Would run:'; printf ' %q' "${cull_cmd[@]}"; printf '\n'

        # Preview the same per-partition FDR and SYN-merge commands that a real
        # run would execute after culling. This is reporting only: no workflow
        # files are created or modified in dry-run mode.
        local folder name input output
        for folder in "${folders[@]}"; do
            name="$(basename "$folder")"
            input="$folder/$culled"
            output="$folder/$fdr"

            echo "[$mode] $name"

            if [[ "$FORCE_FDR" != "True" && -d "$output" ]] && \
               find "$output" -type f -name '*_fdrstats.tsv' -print -quit | grep -q .; then
                echo "  would skip FDR: existing *_fdrstats.tsv found"
            else
                printf '  would run: python3 %q -i %q -o %q\n' "$FDR_SCRIPT" "$input" "$output"
            fi

            if [[ "$SKIP_MERGE" == "True" ]]; then
                echo "  would skip SYN merge: requested"
                continue
            fi

            if [[ "$FORCE_MERGE" != "True" && -d "$output" ]] && \
               find "$output" -type f -name '*_withsyn.tsv' -print -quit | grep -q .; then
                echo "  would skip SYN merge: existing *_withsyn.tsv found"
            else
                printf '  would run: python3 %q -w %q --fdr-dir %q --syn-dir %q\n' \
                    "$MERGER" "$folder" "$fdr" "$SYN_DIR"
            fi
        done
        return
    fi

    if [[ "$AUDIT_ONLY" == "True" ]]; then
        "${cull_cmd[@]}"
        return
    fi

    local log_dir="$DATASET_DIR/${mode}_cull_fdr_merge_logs"
    mkdir -p "$log_dir"
    local stamp main_log summary
    stamp="$(date +%Y%m%d_%H%M%S)"
    main_log="$log_dir/cull_fdr_merge_${stamp}.log"
    summary="$log_dir/postprocess_summary_${stamp}.tsv"
    printf 'partition_folder\tstage\tstatus\tmessage\n' > "$summary"

    "${cull_cmd[@]}" 2>&1 | tee "$main_log"

    local folder name input output plog msg
    for folder in "${folders[@]}"; do
        name="$(basename "$folder")"
        input="$folder/$culled"
        output="$folder/$fdr"
        plog="$log_dir/${name}.log"
        : > "$plog"

        if [[ ! -d "$input" ]] || ! find "$input" -type f -name '*_classical_guided.tsv' -print -quit | grep -q .; then
            msg="missing classical-guided input"
            echo "WARNING [$name]: $msg" | tee -a "$plog" "$main_log"
            printf '%s\tfdr\tskipped\t%s\n' "$name" "$msg" >> "$summary"
            continue
        fi

        if [[ "$FORCE_FDR" != "True" && -d "$output" ]] && \
           find "$output" -type f -name '*_fdrstats.tsv' -print -quit | grep -q .; then
            printf '%s\tfdr\tskipped_existing\tFDR exists\n' "$name" >> "$summary"
        else
            python3 "$FDR_SCRIPT" -i "$input" -o "$output" 2>&1 | tee -a "$plog" "$main_log"
            printf '%s\tfdr\tprocessed\tFDR completed\n' "$name" >> "$summary"
        fi

        if [[ "$SKIP_MERGE" == "True" ]]; then
            printf '%s\tmerge\tskipped\trequested\n' "$name" >> "$summary"
            continue
        fi

        if [[ "$FORCE_MERGE" != "True" && -d "$output" ]] && \
           find "$output" -type f -name '*_withsyn.tsv' -print -quit | grep -q .; then
            printf '%s\tmerge\tskipped_existing\twithsyn exists\n' "$name" >> "$summary"
            continue
        fi

        python3 "$MERGER" -w "$folder" --fdr-dir "$fdr" --syn-dir "$SYN_DIR" 2>&1 | tee -a "$plog" "$main_log"
        printf '%s\tmerge\tprocessed\tSYN merge completed\n' "$name" >> "$summary"
    done

    echo "Summary: $summary"
}

# -----------------------------------------------------------------------------
# Stages 5-7: collect -> compare -> plot
# -----------------------------------------------------------------------------
downstream_mode() {
    local mode="$1"
    local fdr down
    if [[ "$mode" == "pair" ]]; then
        fdr="$PAIR_FDR"
        down="$PAIR_DOWN"
    else
        fdr="$GROUP_FDR"
        down="$GROUP_DOWN"
    fi

    local audit="$DATASET_DIR/$down/classic_vs_${mode}_fdr_audit"
    local delta="$audit/fdr_delta_interval_audit"
    local core="$DATASET_DIR/$down/scan_core_peptide_fdr_summary.tsv"
    local conflict="$DATASET_DIR/$down/${mode}_scan_core_conflicts.tsv"
    local wide="$audit/classic_vs_binned_scan_peptide_fdr_wide.tsv"

    echo
    echo "=== Stages 5-7: collect/compare/plot [$mode] ==="
    echo "Downstream: $DATASET_DIR/$down"

    local log_dir="" log_file=""
    if [[ "$DRY_RUN" != "True" ]]; then
        log_dir="$DATASET_DIR/${mode}_collect_compare_plot_logs"
        mkdir -p "$log_dir"
        log_file="$log_dir/run_$(date +%Y%m%d_%H%M%S).log"
    fi

    local -a cmd
    if [[ "$SKIP_COLLECTOR" != "True" ]]; then
        cmd=(
            python3 "$COLLECTOR"
            --dataset_dir "$DATASET_DIR"
            --partition_type "$mode"
            --fdr_dir_name "$fdr"
            --pattern "$PATTERN"
            --downstream_dir "$down"
        )
        if [[ "$DRY_RUN" == "True" ]]; then
            run_or_preview "${cmd[@]}"
        else
            "${cmd[@]}" 2>&1 | tee -a "$log_file"
        fi
    fi

    if [[ "$DRY_RUN" != "True" ]]; then
        [[ -f "$core" ]] || { echo "ERROR: missing core summary: $core" >&2; exit 1; }
        [[ -f "$conflict" ]] || { echo "ERROR: missing conflict audit: $conflict" >&2; exit 1; }
        (( $(wc -l < "$conflict") <= 1 )) || {
            echo "ERROR: nonempty conflict audit: $conflict" >&2
            exit 1
        }
    fi

    cmd=(
        python3 "$COMPARE"
        --base_withsyn "$REFERENCE"
        --binned_summary "$core"
        --outdir "$audit"
        --thresholds "$THRESHOLDS"
    )
    if [[ "$DRY_RUN" == "True" ]]; then
        run_or_preview "${cmd[@]}"
    else
        "${cmd[@]}" 2>&1 | tee -a "$log_file"
    fi

    if [[ "$SKIP_PLOTS" != "True" ]]; then
        cmd=(python3 "$PLOT" --wide_table "$wide" --outdir "$delta")
        [[ "$INCLUDE_UNSHARED" != "True" ]] && cmd+=(--shared_only)
        if [[ "$DRY_RUN" == "True" ]]; then
            run_or_preview "${cmd[@]}"
        else
            "${cmd[@]}" 2>&1 | tee -a "$log_file"
        fi
    fi
}

# -----------------------------------------------------------------------------
# Stage 8: derive/apply correction formula -> corrected peptide/protein rollup
# -----------------------------------------------------------------------------
correction_mode() {
    local mode="$1"
    local down audit shared condition kind_label

    if [[ "$mode" == "pair" ]]; then
        down="$PAIR_DOWN"
        audit="$DATASET_DIR/$down/classic_vs_pair_fdr_audit"
        kind_label="pairwise"
    else
        down="$GROUP_DOWN"
        audit="$DATASET_DIR/$down/classic_vs_group_fdr_audit"
        kind_label="group"
    fi

    shared="$audit/classic_vs_binned_scan_peptide_fdr_wide.tsv"
    condition="$(basename "$DATASET_DIR")_${kind_label}"

    echo
    echo "=== Stage 8: correction formula/application/rollup [$mode] ==="
    echo "Downstream:  $DATASET_DIR/$down"
    echo "Shared TSV:  $shared"
    echo "Full run:    $FULL_DIR"
    echo "Condition:   $condition"

    if [[ "$DRY_RUN" != "True" ]]; then
        [[ -f "$shared" ]] || {
            echo "ERROR: missing classical-vs-binned wide table required for correction: $shared" >&2
            exit 1
        }
    fi

    local -a cmd=(
        bash "$CORRECTION_RUNNER"
        --project-root "$PROJECT_ROOT"
        --downstream-dir "$DATASET_DIR/$down"
        --full-dir "$FULL_DIR"
        --shared-tsv "$shared"
        --condition-name "$condition"
        --thresholds "$THRESHOLDS"
    )
    [[ "$OVERWRITE_CORRECTION" == "True" ]] && cmd+=(--overwrite)
    [[ "$SKIP_CORRECTION_DERIVE" == "True" ]] && cmd+=(--skip-derive)
    [[ "$SKIP_CORRECTION_APPLY" == "True" ]] && cmd+=(--skip-apply)
    [[ "$DRY_RUN" == "True" ]] && cmd+=(--dry-run)

    if [[ "$DRY_RUN" == "True" ]]; then
        run_or_preview "${cmd[@]}"
    else
        "${cmd[@]}"
    fi
}

# Execute stages in the same order as the current validated shell stack.
if [[ "$SKIP_PREP" != "True" ]]; then
    for mode in "${ACTIVE[@]}"; do
        prepare_mode "$mode"
    done
fi

if [[ "$SKIP_CULL" != "True" ]]; then
    for mode in "${ACTIVE[@]}"; do
        cull_fdr_merge_mode "$mode"
    done
fi

if [[ "$AUDIT_ONLY" == "True" ]]; then
    echo
    echo "Audit-only run complete; downstream stages intentionally skipped."
    exit 0
fi

if [[ "$SKIP_DOWN" != "True" ]]; then
    for mode in "${ACTIVE[@]}"; do
        downstream_mode "$mode"
    done
fi

if [[ "$CORRECTION_ENABLED" == "True" ]]; then
    for mode in "${ACTIVE[@]}"; do
        correction_mode "$mode"
    done
elif [[ "$SKIP_CORRECTION" != "True" && -z "$FULL_DIR" ]]; then
    echo
    echo "Stage 8 skipped: no --full-dir supplied."
fi

echo
echo "Binning downstream workflow completed."
