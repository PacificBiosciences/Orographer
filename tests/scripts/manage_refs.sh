#!/usr/bin/env bash
# Resolve a reference FASTA for e2e tests. Prints exactly one path on stdout (absolute).
# Progress and errors go to stderr. Takes no arguments.
set -euo pipefail

if [ "$#" -ne 0 ]; then
  echo "usage: $(basename "$0") (no arguments)" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SMALL="$REPO_ROOT/tests/data/inputs/hg38_small.fa.gz"
URL="https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz"
# Subset of UCSC hg38 kept when building $SMALL (matches typical e2e BAM contigs).
MANAGED_CHROMS="chr2 chr5 chr6 chr7 chrX chr16"

log() {
  printf '%s\n' "$*" >&2
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    log "manage_refs.sh: '$1' is required but not found in PATH"
    exit 1
  }
}

if [ -f "$SMALL" ]; then
  require_cmd samtools
  if [ ! -f "${SMALL}.fai" ] || [ ! -f "${SMALL}.gzi" ]; then
    log "Indexing $SMALL with samtools faidx..."
    samtools faidx "$SMALL"
  fi
  printf '%s\n' "$SMALL"
  exit 0
fi

require_cmd wget
require_cmd samtools
require_cmd bgzip
require_cmd gunzip
require_cmd awk

mkdir -p "$REPO_ROOT/tests/data/inputs"
TMP_DL="$(mktemp "${TMPDIR:-/tmp}/orographer-hg38.fa.gz.XXXXXX")"
cleanup() {
  rm -f "$TMP_DL"
}
trap cleanup EXIT

log "Downloading reference (hg38), subsetting, and writing $SMALL ..."
wget -nv -O "$TMP_DL" "$URL"

SMALL_TMP="${SMALL}.tmp.$$"
log "Filtering FASTA to chromosomes: $MANAGED_CHROMS"
log "Recompressing with bgzip (required for samtools faidx on compressed FASTA)..."
# UCSC headers look like ">chr2" or ">chr2  ..."; keep exact primary contig names only.
export MANAGED_CHROMS
# shellcheck disable=SC2016
gunzip -c "$TMP_DL" | awk '
  BEGIN {
    ok = 0
    n = split(ENVIRON["MANAGED_CHROMS"], _w, /[[:space:]]+/)
    for (i = 1; i <= n; i++) if (_w[i] != "") allow[_w[i]] = 1
  }
  /^>/ {
    name = substr($0, 2)
    sub(/ .*/, "", name)
    ok = (name in allow)
  }
  ok
' | bgzip -c >"$SMALL_TMP"
rm -f "$TMP_DL"
trap - EXIT
mv "$SMALL_TMP" "$SMALL"

log "Indexing $SMALL with samtools faidx..."
samtools faidx "$SMALL"

printf '%s\n' "$SMALL"
