#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BUILDER="${REPO_ROOT}/scripts/binning/peptide_bin_pair_builder.py"

RUN_NAME=""
CONDITION=""
PEPTIDE_FASTA=""
BINNING_DATA_ROOT=""
NUM_GROUPS="9"
OUTPUT_MODE="both"
SEED=""
DATASET_NAME=""
BUILD_MODE="main"
SIZE_TARGETS="1.2M,600k,300k,150k,75k"
SIZE_TARGETS_SET=0
SIZE_GROUPS=""
SIZE_LABELS=""
DRY_RUN=0

usage() {
cat <<'EOF'
Usage:
  prepare_binning_databases.sh --run-name NAME --peptide-fasta FILE [options]

Prepare peptide database partitions directly in a Binning-compatible data tree.

Required:
  --run-name NAME            Binning run folder name, e.g. Soil_1.
  --peptide-fasta FILE       Input peptide FASTA.

Destination:
  --binning-data-root DIR    Path to the Binning repository data directory.
                             If omitted, tries ../Binning/data relative to
                             FASTA_Processing.
  --condition NAME           Dataset/condition name. If omitted, inferred from
                             --run-name by removing a trailing _<number>.

Build selection:
  --build MODE               main, size-series, or all (default: main).

Main Individual/Pairwise build:
  --dataset-name NAME        Main dataset folder name.
                             Default: <Condition>_Binning.
  --num-groups INT           Number of Individual bins (default: 9).
  --output-mode MODE         groups, pairs, or both (default: both).

Individual-bin size series:
  --size-targets LIST        Comma-separated target peptides per Individual bin.
                             Default: 1.2M,600k,300k,150k,75k.
                             If --size-groups is omitted, each target is
                             converted to the nearest whole number of groups
                             using the FASTA record count.
  --size-groups LIST         Optional comma-separated exact group counts. When
                             supplied, these counts are used directly instead
                             of being calculated from peptide targets.
                             Any number of entries may be supplied.
                             Example:
                               --size-groups 10,20,40
  --size-labels LIST         Optional comma-separated output labels. Labels are
                             written as <Condition>_Groups_<LABEL>.
                             If omitted, labels are derived from --size-targets.
                             Custom labels may contain letters, numbers, ".",
                             "_", and "-".
                             Examples:
                               --size-labels Large,Medium,Small
                               --size-labels 1.2M,600k,300k,150k,75k

                             With --size-groups, label count must match group
                             count. Without --size-groups, label count must
                             match --size-targets.

                             If exact --size-groups are supplied without
                             --size-labels or an explicit --size-targets list,
                             labels default to G<group-count>, e.g. G10,G20,G40.

Shared:
  --seed INT                 Optional reproducible shuffle seed. For size-series
                             builds the same seed is passed to every target.
  --dry-run                  Print planned commands without executing them.
  -h, --help                 Show this help.

Examples:
  # Standard 9-Individual / 36-Pairwise database:
  prepare_binning_databases.sh \
      --run-name Soil_1 \
      --peptide-fasta /path/to/Soil_clustered_peptide.fasta

  # Size series using target peptides/bin (group counts calculated):
  prepare_binning_databases.sh \
      --run-name Soil_1 \
      --peptide-fasta /path/to/Soil_clustered_peptide.fasta \
      --build size-series

  # Size series using exact group counts (paper-style control):
  prepare_binning_databases.sh \
      --run-name Soil_1 \
      --peptide-fasta /path/to/Soil_clustered_peptide.fasta \
      --build size-series \
      --size-groups 11,22,44,89,177

  # Exact group counts with custom dataset labels:
  prepare_binning_databases.sh \
      --run-name Soil_1 \
      --peptide-fasta /path/to/Soil_clustered_peptide.fasta \
      --build size-series \
      --size-groups 10,20,40 \
      --size-labels Small,Medium,Large

  # Build both standard and size-series datasets:
  prepare_binning_databases.sh \
      --run-name Soil_1 \
      --peptide-fasta /path/to/Soil_clustered_peptide.fasta \
      --build all
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-name) RUN_NAME="$2"; shift 2 ;;
        --condition) CONDITION="$2"; shift 2 ;;
        --peptide-fasta) PEPTIDE_FASTA="$2"; shift 2 ;;
        --binning-data-root) BINNING_DATA_ROOT="$2"; shift 2 ;;
        --num-groups) NUM_GROUPS="$2"; shift 2 ;;
        --output-mode) OUTPUT_MODE="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --dataset-name) DATASET_NAME="$2"; shift 2 ;;
        --build) BUILD_MODE="$2"; shift 2 ;;
        --size-targets) SIZE_TARGETS="$2"; SIZE_TARGETS_SET=1; shift 2 ;;
        --size-groups) SIZE_GROUPS="$2"; shift 2 ;;
        --size-labels) SIZE_LABELS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$RUN_NAME" ]] || { echo "ERROR: --run-name is required." >&2; exit 2; }
[[ -n "$PEPTIDE_FASTA" ]] || { echo "ERROR: --peptide-fasta is required." >&2; exit 2; }
[[ -f "$PEPTIDE_FASTA" ]] || { echo "ERROR: peptide FASTA not found: $PEPTIDE_FASTA" >&2; exit 2; }

case "$BUILD_MODE" in
    main|size-series|all) ;;
    *) echo "ERROR: --build must be main, size-series, or all." >&2; exit 2 ;;
esac

case "$OUTPUT_MODE" in
    groups|pairs|both) ;;
    *) echo "ERROR: --output-mode must be groups, pairs, or both." >&2; exit 2 ;;
esac

[[ "$NUM_GROUPS" =~ ^[0-9]+$ ]] && (( NUM_GROUPS >= 1 )) || {
    echo "ERROR: --num-groups must be a positive integer." >&2
    exit 2
}

if [[ -z "$CONDITION" ]]; then
    if [[ "$RUN_NAME" =~ ^(.+)_([0-9]+)$ ]]; then
        CONDITION="${BASH_REMATCH[1]}"
    else
        CONDITION="$RUN_NAME"
    fi
fi

if [[ -z "$DATASET_NAME" ]]; then
    DATASET_NAME="${CONDITION}_Binning"
fi

if [[ -z "$BINNING_DATA_ROOT" ]]; then
    SIBLING_BINNING_ROOT="$(cd "${REPO_ROOT}/.." && pwd)/Binning/data"
    if [[ -d "$SIBLING_BINNING_ROOT" ]]; then
        BINNING_DATA_ROOT="$SIBLING_BINNING_ROOT"
    else
        echo "ERROR: --binning-data-root was not provided and sibling Binning/data was not found:" >&2
        echo "  $SIBLING_BINNING_ROOT" >&2
        echo "Provide --binning-data-root explicitly." >&2
        exit 2
    fi
fi

RUN_ROOT="${BINNING_DATA_ROOT}/${RUN_NAME}"

# Count peptide FASTA records once for size-series group-count calculations.
PEPTIDE_COUNT="$(grep -c '^>' "$PEPTIDE_FASTA" || true)"
[[ "$PEPTIDE_COUNT" =~ ^[0-9]+$ ]] && (( PEPTIDE_COUNT > 0 )) || {
    echo "ERROR: no FASTA records were found in: $PEPTIDE_FASTA" >&2
    exit 2
}

check_output_dir() {
    local out_dir="$1"
    if [[ -e "$out_dir" ]] && [[ -n "$(find "$out_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
        echo "ERROR: output dataset directory already exists and is not empty:" >&2
        echo "  $out_dir" >&2
        echo "Refusing to overwrite an existing Binning dataset." >&2
        exit 2
    fi
}

run_builder() {
    local out_dir="$1"
    local groups="$2"
    local mode="$3"

    check_output_dir "$out_dir"

    local -a cmd=(
        python3 "$BUILDER"
        --input_fasta "$PEPTIDE_FASTA"
        --out_dir "$out_dir"
        --num_groups "$groups"
        --output-mode "$mode"
        --group_subfolders
    )
    if [[ -n "$SEED" ]]; then
        cmd+=(--seed "$SEED")
    fi

    printf '  '
    printf '%q ' "${cmd[@]}"
    printf '\n'

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$(dirname "$out_dir")"
        "${cmd[@]}"
    fi
}

parse_size_target() {
    local raw="$1"
    python3 - "$raw" <<'PY'
import re
import sys

raw = sys.argv[1].strip()
m = re.fullmatch(r"(?i)(\d+(?:\.\d+)?)\s*([kKmM]?)", raw)
if not m:
    raise SystemExit(2)
value = float(m.group(1))
suffix = m.group(2).lower()
if suffix == "k":
    value *= 1_000
elif suffix == "m":
    value *= 1_000_000
target = int(round(value))
if target < 1:
    raise SystemExit(2)
print(target)
PY
}

canonical_size_label() {
    local target="$1"
    python3 - "$target" <<'PY'
import sys
n = int(sys.argv[1])
if n >= 1_000_000 and n % 100_000 == 0:
    s = f"{n / 1_000_000:.1f}".rstrip("0").rstrip(".")
    print(f"{s}M")
elif n >= 1_000 and n % 1_000 == 0:
    print(f"{n // 1_000}k")
else:
    print(str(n))
PY
}

nearest_group_count() {
    local peptide_count="$1"
    local target="$2"
    # Integer half-up rounding of peptide_count / target, minimum one group.
    local groups=$(( (peptide_count + target / 2) / target ))
    if (( groups < 1 )); then
        groups=1
    fi
    echo "$groups"
}

echo "FASTA records: $PEPTIDE_COUNT"
echo "Binning run root: $RUN_ROOT"

if [[ "$BUILD_MODE" == "main" || "$BUILD_MODE" == "all" ]]; then
    MAIN_OUT="${RUN_ROOT}/${DATASET_NAME}"
    echo
    echo "Main Individual/Pairwise build:"
    echo "  Dataset: $DATASET_NAME"
    echo "  Groups:  $NUM_GROUPS"
    echo "  Mode:    $OUTPUT_MODE"
    run_builder "$MAIN_OUT" "$NUM_GROUPS" "$OUTPUT_MODE"
fi

if [[ "$BUILD_MODE" == "size-series" || "$BUILD_MODE" == "all" ]]; then
    echo
    echo "Individual-bin size-series builds:"

    IFS=',' read -r -a TARGET_ITEMS <<< "$SIZE_TARGETS"
    [[ "${#TARGET_ITEMS[@]}" -gt 0 ]] || {
        echo "ERROR: --size-targets did not contain any targets." >&2
        exit 2
    }

    GROUP_ITEMS=()
    if [[ -n "$SIZE_GROUPS" ]]; then
        IFS=',' read -r -a GROUP_ITEMS <<< "$SIZE_GROUPS"
        [[ "${#GROUP_ITEMS[@]}" -gt 0 ]] || {
            echo "ERROR: --size-groups did not contain any group counts." >&2
            exit 2
        }
    fi

    LABEL_ITEMS=()
    if [[ -n "$SIZE_LABELS" ]]; then
        IFS=',' read -r -a LABEL_ITEMS <<< "$SIZE_LABELS"
        [[ "${#LABEL_ITEMS[@]}" -gt 0 ]] || {
            echo "ERROR: --size-labels did not contain any labels." >&2
            exit 2
        }
    fi

    if [[ -n "$SIZE_GROUPS" ]]; then
        SERIES_COUNT="${#GROUP_ITEMS[@]}"
        if [[ -n "$SIZE_LABELS" ]] && [[ "${#LABEL_ITEMS[@]}" -ne "$SERIES_COUNT" ]]; then
            echo "ERROR: --size-labels must contain exactly one label for each --size-groups entry." >&2
            echo "  Groups: ${#GROUP_ITEMS[@]}" >&2
            echo "  Labels: ${#LABEL_ITEMS[@]}" >&2
            exit 2
        fi
        if [[ -z "$SIZE_LABELS" && "$SIZE_TARGETS_SET" -eq 1 ]] && [[ "${#TARGET_ITEMS[@]}" -ne "$SERIES_COUNT" ]]; then
            echo "ERROR: when explicit --size-targets are used as labels with --size-groups, their entry counts must match." >&2
            echo "  Targets: ${#TARGET_ITEMS[@]}" >&2
            echo "  Groups:  ${#GROUP_ITEMS[@]}" >&2
            echo "Alternatively provide --size-labels, or omit --size-targets to use automatic G<count> labels." >&2
            exit 2
        fi
    else
        SERIES_COUNT="${#TARGET_ITEMS[@]}"
        if [[ -n "$SIZE_LABELS" ]] && [[ "${#LABEL_ITEMS[@]}" -ne "$SERIES_COUNT" ]]; then
            echo "ERROR: --size-labels must contain exactly one label for each --size-targets entry." >&2
            echo "  Targets: ${#TARGET_ITEMS[@]}" >&2
            echo "  Labels:  ${#LABEL_ITEMS[@]}" >&2
            exit 2
        fi
    fi

    for (( i=0; i<SERIES_COUNT; i++ )); do
        target=""
        target_desc=""

        if [[ -z "$SIZE_GROUPS" || "$SIZE_TARGETS_SET" -eq 1 ]]; then
            raw_target="${TARGET_ITEMS[$i]//[[:space:]]/}"
            if ! target="$(parse_size_target "$raw_target")"; then
                echo "ERROR: invalid size target: $raw_target" >&2
                exit 2
            fi
        fi

        if [[ -n "$SIZE_LABELS" ]]; then
            label="${LABEL_ITEMS[$i]//[[:space:]]/}"
            [[ -n "$label" ]] || {
                echo "ERROR: size-series labels may not be empty." >&2
                exit 2
            }
            [[ "$label" =~ ^[A-Za-z0-9._-]+$ ]] || {
                echo "ERROR: invalid size-series label: ${LABEL_ITEMS[$i]}" >&2
                echo 'Labels may contain only letters, numbers, ".", "_", and "-".' >&2
                exit 2
            }
        elif [[ -n "$SIZE_GROUPS" && "$SIZE_TARGETS_SET" -eq 0 ]]; then
            raw_groups="${GROUP_ITEMS[$i]//[[:space:]]/}"
            label="G${raw_groups}"
        else
            label="$(canonical_size_label "$target")"
        fi

        if [[ -n "$SIZE_GROUPS" ]]; then
            groups="${GROUP_ITEMS[$i]//[[:space:]]/}"
            [[ "$groups" =~ ^[0-9]+$ ]] && (( groups >= 1 )) || {
                echo "ERROR: invalid exact group count for $label: ${GROUP_ITEMS[$i]}" >&2
                exit 2
            }
            source_desc="exact group count"
            if [[ -n "$target" ]]; then
                target_desc="$target peptides/bin label reference"
            elif [[ -n "$SIZE_LABELS" ]]; then
                target_desc="custom label"
            else
                target_desc="automatic group-count label"
            fi
        else
            groups="$(nearest_group_count "$PEPTIDE_COUNT" "$target")"
            source_desc="calculated from target size"
            target_desc="$target peptides/bin target"
        fi

        actual=$(( (PEPTIDE_COUNT + groups / 2) / groups ))
        out_dir="${RUN_ROOT}/${CONDITION}_Groups_${label}"

        echo
        echo "  Label:       $label"
        echo "  Basis:       $target_desc"
        echo "  Groups:      $groups ($source_desc)"
        echo "  Approx/bin:  $actual"
        echo "  Dataset:     ${CONDITION}_Groups_${label}"
        run_builder "$out_dir" "$groups" "groups"
    done
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo
    echo "Dry run only; no commands executed."
else
    echo
    echo "Binning database preparation complete."
    echo "Run root: $RUN_ROOT"
fi
