#!/usr/bin/env bash
set -euo pipefail

# Complex-SV smoke runner.
#
# Fixture data (BAMs, VCF, GTF) must already be present under
# tests/data/inputs/complex_sv. Run from the repo root or let the script
# cd there automatically.
#
# It is meant to verify that commands execute; later we can turn the fixture
# set and command matrix into a real pytest/integration suite.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

INPUT_DIR="tests/data/inputs/complex_sv"
BAM_DIR="$INPUT_DIR"
VCF_DIR="$INPUT_DIR/vcf"
ANNOTATION_DIR="$INPUT_DIR/annotation"
REGIONS_BED="$INPUT_DIR/complex_sv_regions.bed"
OUTPUT_ROOT="tests/data/expected_outputs/complex_sv"

VCF_FILE="$VCF_DIR/sawfish.vcf.gz"
GTF_FILE="$ANNOTATION_DIR/gencode.v47.genes.gtf.gz"

# A full hg38 reference is required for these fixtures. The repo's
# hg38_small.fa.gz only covers the paraviewer smoke data.
REF_FASTA="${REF_FASTA:-${HG38_FASTA:-/home/UNIXHOME/jbelyeu/references/grch38-reference-genome-gencode-v1.fa}}"

# Use the editable/dev command if available in PATH, otherwise fall back to the
# local module invocation. Override with OROGRAPHER=/path/to/orographer.
if [ -n "${OROGRAPHER:-}" ]; then
  read -r -a OROGRAPHER_CMD <<< "$OROGRAPHER"
else
  if command -v orographer >/dev/null 2>&1; then
    OROGRAPHER_CMD=("$(command -v orographer)")
  else
    OROGRAPHER_CMD=("python" "-m" "orographer")
  fi
fi

if [ -z "$REF_FASTA" ] || [ ! -f "$REF_FASTA" ]; then
  cat >&2 <<EOF
Set REF_FASTA or HG38_FASTA to a full hg38 FASTA before running this script.

Example:
  REF_FASTA=/path/to/hg38.fa tests/scripts/smoketest_tmp.sh
EOF
  exit 1
fi

mkdir -p "$OUTPUT_ROOT"

# The VCF may be regular gzip rather than bgzip; tabix requires bgzip.
# Re-compress and index if the index is missing.
if command -v tabix >/dev/null 2>&1 && command -v bgzip >/dev/null 2>&1; then
  if [ ! -f "${VCF_FILE}.tbi" ]; then
    echo "Preparing bgzip-indexed VCF: $VCF_FILE"
    tmp_vcf="$(mktemp --suffix=.vcf.gz)"
    gzip -dc "$VCF_FILE" | bgzip -c > "$tmp_vcf"
    mv "$tmp_vcf" "$VCF_FILE"
    tabix -f -p vcf "$VCF_FILE"
  fi
fi

if [ "${SKIP_GTF:-false}" != "true" ]; then
  if ! command -v tabix >/dev/null 2>&1 || ! command -v bgzip >/dev/null 2>&1; then
    cat >&2 <<EOF
tabix and bgzip are required to index the GTF for orographer.
Run with SKIP_GTF=true to smoke-test without gene annotations.
EOF
    exit 1
  fi

  if [ ! -f "$GTF_FILE" ]; then
    echo "GTF not found at $GTF_FILE; run with SKIP_GTF=true to skip." >&2
    exit 1
  fi

  if [ ! -f "${GTF_FILE}.tbi" ]; then
    echo "Preparing bgzip-indexed GTF: $GTF_FILE"
    if ! tabix -f -p gff "$GTF_FILE"; then
      echo "Direct tabix indexing failed; rebuilding GTF as sorted bgzip."
      tmp_gtf="$(mktemp)"
      gzip -dc "$GTF_FILE" > "$tmp_gtf"
      {
        grep '^#' "$tmp_gtf" || true
        grep -v '^#' "$tmp_gtf" | sort -k1,1 -k4,4n
      } | bgzip -c > "$GTF_FILE"
      rm -f "$tmp_gtf"
      tabix -f -p gff "$GTF_FILE"
    fi
  fi
fi

echo "Writing generated region manifest: $REGIONS_BED"
{
  printf "#chrom\tstart\tend\tname\tbam\truntype\n"
  for bam in "$BAM_DIR"/*.bam; do
    [ -f "$bam" ] || continue
    bam_base="$(basename "$bam")"
    stem="${bam_base%.bam}"
    region_stem="$stem"
    runtype="novcf"
    case "$region_stem" in
      *_vcfonly)
        runtype="vcfonly"
        region_stem="${region_stem%_vcfonly}"
        ;;
      *_vcf)
        runtype="vcf"
        region_stem="${region_stem%_vcf}"
        ;;
    esac

    chrom="$(printf "%s" "$region_stem" | cut -d_ -f1)"
    start="$(printf "%s" "$region_stem" | cut -d_ -f2)"
    end="$(printf "%s" "$region_stem" | cut -d_ -f3)"
    name="$region_stem"

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$chrom" "$start" "$end" "$name" "$bam_base" "$runtype"
  done
} > "$REGIONS_BED"

echo "Running complex-SV smoke commands into $OUTPUT_ROOT"
while IFS=$'\t' read -r chrom start end name bam_base runtype || [ -n "${chrom:-}" ]; do
  [ -z "${chrom:-}" ] && continue
  [[ "$chrom" =~ ^# ]] && continue

  bam="$BAM_DIR/$bam_base"
  # BAM filename coordinates are 0-based (BED convention); orographer --coord is
  # 1-based inclusive, matching the project convention in run_end_to_end_tests.sh.
  coord="${chrom}:$((start + 1))-${end}"
  outdir="$OUTPUT_ROOT/$name/$runtype"
  mkdir -p "$outdir"

  cmd=(
    "${OROGRAPHER_CMD[@]}"
    plot
    --bam "$bam"
    --coord "$coord"
    --region-type complex_sv
    --ref "$REF_FASTA"
    --outdir "$outdir"
  )

  # Match the old script's intent: _vcf had both variant-readnames and VCF,
  # _vcfonly had VCF only. orographer does not use variant-readnames, so both
  # VCF-backed modes simply pass --vcf.
  if [ "$runtype" = "vcf" ] || [ "$runtype" = "vcfonly" ]; then
    cmd+=(--vcf "$VCF_FILE")
  fi

  if [ "${SKIP_GTF:-false}" != "true" ] && [ -f "${GTF_FILE}.tbi" ]; then
    cmd+=(--gtf "$GTF_FILE")
  fi

  echo
  echo "Input: $bam_base"
  echo "Coord: $coord"
  printf 'Command:'
  printf ' %q' "${cmd[@]}"
  printf '\n'

  "${cmd[@]}"
done < "$REGIONS_BED"

echo
echo "Complex-SV smoke commands finished."
