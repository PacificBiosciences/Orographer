import gzip
import json
import subprocess
import sys
from pathlib import Path

from tests.helpers.plot_json_schema import (
    assert_genomic_ranges_in_docs_json,
    assert_orig_region_bounds_in_docs_json,
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
            _bed_chrom, start0, end, _name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            return f"{_bed_chrom}:{start0 + 1}-{end}"
    raise AssertionError(f"No BED line for chrom={chrom} in {bed_path}")


def test_cli_complex_sv_creates_outputs(tmp_path: Path):
    root = Path("tests/data/inputs")

    bam_path = root / "paraviewer/m84106_250429_131828_s4.paraphase.small.bam"
    ref_path = root / "hg38_small.fa.gz"
    bed_path = root / "paraviewer/m84106_250429_131828_s4_regions.bed"

    # Choose a known chromosome present in the fixture BED.
    coord = _bed_first_coord_to_cli_region(bed_path, chrom="chr16")

    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam_path),
        "--coord",
        coord,
        "--region-type",
        "complex_sv",
        "--ref",
        str(ref_path),
        "--outdir",
        str(tmp_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stderr

    html_files = list(tmp_path.glob("*_bokeh.html"))
    json_files = list(tmp_path.glob("*_bokeh.json.gz"))
    assert len(html_files) >= 1
    assert len(json_files) >= 1

    json_path = json_files[0]
    with gzip.open(json_path, "rt", encoding="utf-8") as f:
        data = json.load(f)

    validate_plot_json_schema(data)
    _chrom, rest = coord.split(":")
    start_1, end_1 = int(rest.split("-")[0]), int(rest.split("-")[1])
    assert_orig_region_bounds_in_docs_json(data["docs_json"], [(start_1, end_1)])
    assert_genomic_ranges_in_docs_json(data["docs_json"], [(start_1, end_1)])
