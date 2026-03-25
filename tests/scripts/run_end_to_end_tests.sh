set -e

# Repo root: venv, ruff, pytest, and tests/data paths are relative to here.
_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
_REPO_ROOT="$(cd "$_SCRIPT_DIR/../.." && pwd)"
cd "$_REPO_ROOT"

# Assumes uv is available. Use a uv-managed venv with CPython so the build
# (hatchling) does not run under GraalPy, which can break editable installs.
# Install: editable orographer + [project] deps (e.g. pysam) + optional [e2e] (pytest).
UV_VENV=".venv_e2e"
if [ ! -d "$UV_VENV" ] || [ ! -x "$UV_VENV/bin/orographer" ]; then
  rm -rf "$UV_VENV"
  uv python install 3.10
  uv venv "$UV_VENV" --python 3.10
  uv pip install -e ".[e2e]" --python "$UV_VENV/bin/python"
fi
uv pip install -e ".[e2e]" --python "$UV_VENV/bin/python"
OROGRAPHER="$UV_VENV/bin/orographer"
PYTHON_E2E="$UV_VENV/bin/python"
RUFF="$UV_VENV/bin/ruff"

# Lint + format (same rules as pyproject [tool.ruff]); must pass before pytest / plot e2e.
"$RUFF" check .
"$RUFF" format --check .


RED=$'\033[0;31m'
GREEN=$'\033[0;32m'
RESET=$'\033[0m'
final_color=$GREEN

REF_FASTA="$("$_SCRIPT_DIR/manage_refs.sh")"
GTF_FILE="tests/data/inputs/small_gencode.gtf.gz"
BED_SUFFIX="_regions.bed"
TMP_TEST_DIR=".tmp_test"
# Set RENEW=true to refresh expected_outputs instead of comparing (e.g. RENEW=true ./run_end_to_end_tests.sh)
RENEW="${RENEW:-false}"
# RENEW="${RENEW:-true}"

"$UV_VENV/bin/pytest"

for runtype_dir in tests/data/inputs/*/; do
  runtype=$(basename "$runtype_dir")
  # Map input subdirectory name to orographer --region-type
  case "$runtype" in
    paraviewer) region_type="paraphase" ;;
    complex_sv)  region_type="complex_sv" ;;
    *)          echo "Skipping unknown runtype: $runtype"; continue ;;
  esac

  for bam in "${runtype_dir}"*.bam; do
    [ -f "$bam" ] || continue
    bam_stem=$(basename "$bam" .bam)
    # Sample-specific BED: first component of bam stem (e.g. HG002.paraphase.small -> HG002_regions.bed)
    sample_id="${bam_stem%%.*}"
    bed_file="${runtype_dir}${sample_id}${BED_SUFFIX}"
    if [ ! -f "$bed_file" ]; then
      echo "No ${sample_id}${BED_SUFFIX} in ${runtype_dir}; skipping BAM $bam_stem."
      continue
    fi

    while IFS=$'\t' read -r chrom start end name || [ -n "$chrom" ]; do
      [ -z "$chrom" ] && continue
      [[ "$chrom" =~ ^# ]] && continue

      coord="${chrom}:$((start + 1))-${end}"

      rm -rf "$TMP_TEST_DIR"
      mkdir -p "$TMP_TEST_DIR"

      cmd="$OROGRAPHER plot \
        --bam $bam \
        --coord $coord \
        --region-type $region_type \
        --ref $REF_FASTA \
        --gtf $GTF_FILE \
        --outdir $TMP_TEST_DIR"
      if ! $cmd; then
        printf "${RED}FAILED (orographer): %s / %s (%s)${RESET}\n" "$bam_stem" "$name" "$coord"
        final_color="$RED"
      fi

      expected_dir="tests/data/expected_outputs/${runtype}/${bam_stem}/${name}"
      if [ "$RENEW" = "true" ]; then
        mkdir -p "$expected_dir"
        rm -f "${expected_dir}"/*_bokeh.html "${expected_dir}"/*_bokeh.json.gz 2>/dev/null || true
        for f in "${TMP_TEST_DIR}"/*_bokeh.html; do
          [ -f "$f" ] && cp "$f" "$expected_dir/"
        done
        for f in "${TMP_TEST_DIR}"/*_bokeh.json.gz; do
          [ -f "$f" ] && cp "$f" "$expected_dir/"
        done
      else
        for f in "${TMP_TEST_DIR}"/*_bokeh.html; do
          if [ -f "$f" ]; then
            expected="${expected_dir}/$(basename "$f")"
            if [ -f "$expected" ]; then
              if ! "$PYTHON_E2E" tests/scripts/compare_files.py -r "$expected" -t "$f"; then
                printf "${RED}FAILED (compare): %s vs %s${RESET}\n" "$expected" "$f"
                final_color="$RED"
              fi
            else
              printf "${RED}Missing expected: %s${RESET}\n" "$expected"
              final_color="$RED"
            fi
          fi
        done
        for f in "${TMP_TEST_DIR}"/*_bokeh.json.gz; do
          if [ -f "$f" ]; then
            expected="${expected_dir}/$(basename "$f")"
            if [ -f "$expected" ]; then
              if ! "$PYTHON_E2E" tests/scripts/compare_files.py -r "$expected" -t "$f"; then
                printf "${RED}FAILED (compare): %s vs %s${RESET}\n" "$expected" "$f"
                final_color="$RED"
              fi
            else
              printf "${RED}Missing expected: %s${RESET}\n" "$expected"
              final_color="$RED"
            fi
          fi
        done
      fi
    done < "$bed_file"
  done

  # Multi-sample test (paraviewer only): both BAMs, m84106 regions, write to multisample/
  if [ "$runtype" = "paraviewer" ]; then
    multisample_bed="${runtype_dir}m84106_250429_131828_s4_regions.bed"
    if [ -f "$multisample_bed" ]; then
      bams=( "${runtype_dir}"*.bam )
      if [ ${#bams[@]} -eq 2 ]; then
        bam1="${bams[0]}"
        bam2="${bams[1]}"
        sample1=$(basename "$bam1" .bam)
        sample2=$(basename "$bam2" .bam)
        while IFS=$'\t' read -r chrom start end name || [ -n "$chrom" ]; do
          [ -z "$chrom" ] && continue
          [[ "$chrom" =~ ^# ]] && continue
          coord="${chrom}:$((start + 1))-${end}"
          rm -rf "$TMP_TEST_DIR"
          mkdir -p "$TMP_TEST_DIR"
          cmd="$OROGRAPHER plot \
            --bam $bam1 \
            --other-bam $bam2 \
            --sample-label $sample1 \
            --other-sample-label $sample2 \
            --coord $coord \
            --region-type paraphase \
            --ref $REF_FASTA \
            --gtf $GTF_FILE \
            --outdir $TMP_TEST_DIR"
          if ! $cmd; then
            printf "${RED}FAILED (orographer multisample): %s (%s)${RESET}\n" "$name" "$coord"
            final_color="$RED"
          fi
          expected_dir="tests/data/expected_outputs/paraviewer/multisample/${name}"
          if [ "$RENEW" = "true" ]; then
            mkdir -p "$expected_dir"
            rm -f "${expected_dir}"/*_bokeh.html "${expected_dir}"/*_bokeh.json.gz 2>/dev/null || true
            for f in "${TMP_TEST_DIR}"/*_bokeh.html; do
              [ -f "$f" ] && cp "$f" "$expected_dir/"
            done
            for f in "${TMP_TEST_DIR}"/*_bokeh.json.gz; do
              [ -f "$f" ] && cp "$f" "$expected_dir/"
            done
          else
            for f in "${TMP_TEST_DIR}"/*_bokeh.html; do
              if [ -f "$f" ]; then
                expected="${expected_dir}/$(basename "$f")"
                if [ -f "$expected" ]; then
                  if ! "$PYTHON_E2E" tests/scripts/compare_files.py -r "$expected" -t "$f"; then
                    printf "${RED}FAILED (compare multisample): %s vs %s${RESET}\n" "$expected" "$f"
                    final_color="$RED"
                  fi
                else
                  printf "${RED}Missing expected: %s${RESET}\n" "$expected"
                  final_color="$RED"
                fi
              fi
            done
            for f in "${TMP_TEST_DIR}"/*_bokeh.json.gz; do
              if [ -f "$f" ]; then
                expected="${expected_dir}/$(basename "$f")"
                if [ -f "$expected" ]; then
                  if ! "$PYTHON_E2E" tests/scripts/compare_files.py -r "$expected" -t "$f"; then
                    printf "${RED}FAILED (compare multisample): %s vs %s${RESET}\n" "$expected" "$f"
                    final_color="$RED"
                  fi
                else
                  printf "${RED}Missing expected: %s${RESET}\n" "$expected"
                  final_color="$RED"
                fi
              fi
            done
          fi
        done < "$multisample_bed"
      else
        echo "Skipping multisample: expected 2 BAMs in ${runtype_dir}, found ${#bams[@]}"
      fi
    fi
  fi
done

if [ "$RENEW" = "true" ]; then
  printf "\n${final_color}Test data renewed${RESET}\n"
else
  printf "\n${final_color}Tests finished${RESET}\n"
fi
