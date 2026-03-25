import gzip
import json
import subprocess
import sys
from pathlib import Path

from tests.helpers.plot_json_schema import (
    assert_genomic_ranges_in_docs_json,
    assert_orig_region_bounds_in_docs_json,
    assert_sample_labels_in_docs_json,
    validate_plot_json_schema,
)


def _bed_first_coord_to_cli_region(bed_path: Path, chrom: str) -> str:
    """
    Convert BED (0-based start, end) into orographer CLI coordinate (1-based inclusive).
    e2e scripts use: coord="${chrom}:$((start + 1))-${end}"
    """
    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] != chrom:
                continue
            bed_chrom, start0, end, _name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            return f"{bed_chrom}:{start0 + 1}-{end}"
    raise AssertionError(f"No BED line for chrom={chrom} in {bed_path}")


def test_multisample_label_propagation(tmp_path: Path):
    root = Path("tests/data/inputs")

    bam_hg002 = root / "paraviewer/HG002.paraphase.small.bam"
    bam_m84106 = root / "paraviewer/m84106_250429_131828_s4.paraphase.small.bam"
    ref_path = root / "hg38_small.fa.gz"
    bed_path = root / "paraviewer/m84106_250429_131828_s4_regions.bed"

    # Pick a region from the shared multisample BED.
    # (Use the first line we can find for a stable chromosome token.)
    coord = _bed_first_coord_to_cli_region(bed_path, chrom="chr5")

    outdir = tmp_path / "out"
    outdir.mkdir(parents=True, exist_ok=True)

    primary_label = "PrimarySample"
    other_label = "OtherSample"

    # Explicitly set which BAM is primary (--bam) vs comparison (--other-bam),
    # so label->track mapping is deterministic.
    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam_hg002),
        "--other-bam",
        str(bam_m84106),
        "--coord",
        coord,
        "--region-type",
        "paraphase",
        "--ref",
        str(ref_path),
        "--outdir",
        str(outdir),
        "--sample-label",
        primary_label,
        "--other-sample-label",
        other_label,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr

    json_files = list(outdir.glob("*_bokeh.json.gz"))
    assert json_files, "Expected at least one *_bokeh.json.gz output"

    json_path = json_files[0]
    with gzip.open(json_path, "rt", encoding="utf-8") as f:
        data = json.load(f)

    validate_plot_json_schema(data)
    _c, rest = coord.split(":")
    s1, e1 = int(rest.split("-")[0]), int(rest.split("-")[1])
    assert_orig_region_bounds_in_docs_json(data["docs_json"], [(s1, e1)])
    assert_genomic_ranges_in_docs_json(data["docs_json"], [(s1, e1)])
    assert_sample_labels_in_docs_json(data["docs_json"], [primary_label, other_label])
