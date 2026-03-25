import subprocess
import sys
from pathlib import Path

from orographer.utils import Region
from tests.helpers.plot_json_schema import (
    assert_genomic_ranges_in_docs_json,
    assert_orig_region_bounds_in_docs_json,
    validate_plot_json_schema,
)


def _bed_first_coord_to_cli_region(bed_path: Path, chrom: str) -> Region:
    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] != chrom:
                continue
            bed_chrom, start0, end, _name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            cli_start1 = start0 + 1
            coordinate_str = f"{bed_chrom}:{cli_start1}-{end}"
            return Region(bed_chrom, cli_start1, end, coordinate_str)
    raise AssertionError(f"No BED line for chrom={chrom} in {bed_path}")


def test_cli_plot_creates_outputs(tmp_path: Path):
    root = Path("tests/data/inputs")
    bam_path = root / "paraviewer/HG002.paraphase.small.bam"
    ref_path = root / "hg38_small.fa.gz"
    bed_path = root / "paraviewer/HG002_regions.bed"

    region = _bed_first_coord_to_cli_region(bed_path, chrom="chr16")

    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam_path),
        "--coord",
        region.coordinate_str,
        "--region-type",
        "paraphase",
        "--ref",
        str(ref_path),
        "--outdir",
        str(tmp_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr

    html_files = list(tmp_path.glob("*_bokeh.html"))
    json_files = list(tmp_path.glob("*_bokeh.json.gz"))

    assert len(html_files) >= 1
    assert len(json_files) >= 1

    import gzip
    import json

    json_path = next(tmp_path.glob("*_bokeh.json.gz"))
    with gzip.open(json_path, "rt", encoding="utf-8") as f:
        data = json.load(f)

    validate_plot_json_schema(data)
    assert_orig_region_bounds_in_docs_json(data["docs_json"], [(region.start, region.end)])
    assert_genomic_ranges_in_docs_json(data["docs_json"], [(region.start, region.end)])


def test_cli_plot_multi_region_filename(tmp_path: Path):
    root = Path("tests/data/inputs")
    bam_path = root / "paraviewer/HG002.paraphase.small.bam"
    ref_path = root / "hg38_small.fa.gz"
    bed_path = root / "paraviewer/HG002_regions.bed"

    region1 = _bed_first_coord_to_cli_region(bed_path, chrom="chr16")
    region2 = _bed_first_coord_to_cli_region(bed_path, chrom="chr2")

    def _sanitize(coord: str) -> str:
        return coord.replace(":", "_").replace("-", "_")

    expected_base = f"{_sanitize(region1.coordinate_str)}_{_sanitize(region2.coordinate_str)}_bokeh"
    expected_html = tmp_path / f"{expected_base}.html"
    expected_json = tmp_path / f"{expected_base}.json.gz"

    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam_path),
        "--coord",
        region1.coordinate_str,
        "--coord",
        region2.coordinate_str,
        "--region-type",
        "paraphase",
        "--ref",
        str(ref_path),
        "--outdir",
        str(tmp_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr

    assert expected_html.exists()
    assert expected_json.exists()

    import gzip
    import json

    with gzip.open(expected_json, "rt", encoding="utf-8") as f:
        data = json.load(f)
    validate_plot_json_schema(data)
    assert_orig_region_bounds_in_docs_json(
        data["docs_json"],
        [(region1.start, region1.end), (region2.start, region2.end)],
    )
    assert_genomic_ranges_in_docs_json(
        data["docs_json"],
        [(region1.start, region1.end), (region2.start, region2.end)],
    )
