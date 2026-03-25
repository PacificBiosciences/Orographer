"""CLI workflow edge cases: multi-coord, prefix+multi-region, deploy failures."""

import subprocess
import sys
from pathlib import Path


def _bed_coord(bed_path: Path, chrom: str) -> str:
    with bed_path.open(encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] != chrom:
                continue
            start0, end = int(parts[1]), int(parts[2])
            return f"{parts[0]}:{start0 + 1}-{end}"
    raise AssertionError(f"No row for {chrom} in {bed_path}")


def test_plot_cli_accepts_multiple_coord_and_succeeds(tmp_path):
    root = Path("tests/data/inputs")
    bam = root / "paraviewer/HG002.paraphase.small.bam"
    ref = root / "hg38_small.fa.gz"
    bed = root / "paraviewer/HG002_regions.bed"
    c1 = _bed_coord(bed, "chr16")
    c2 = _bed_coord(bed, "chr2")
    outdir = tmp_path / "out_multi_coord"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam),
        "--coord",
        c1,
        "--coord",
        c2,
        "--region-type",
        "paraphase",
        "--ref",
        str(ref),
        "--outdir",
        str(outdir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stderr
    html = list(outdir.glob("*_bokeh.html"))
    json_gz = list(outdir.glob("*_bokeh.json.gz"))
    assert len(html) == 1 and len(json_gz) == 1


def test_plot_prefix_multi_region_output_basename(tmp_path):
    root = Path("tests/data/inputs")
    bam = root / "paraviewer/HG002.paraphase.small.bam"
    ref = root / "hg38_small.fa.gz"
    bed = root / "paraviewer/HG002_regions.bed"
    c1 = _bed_coord(bed, "chr16")
    c2 = _bed_coord(bed, "chr2")
    prefix = "pfx_multi"

    def _san(c: str) -> str:
        return c.replace(":", "_").replace("-", "_")

    expected = f"{prefix}_{_san(c1)}_{_san(c2)}_bokeh.html"
    outdir = tmp_path / "out_prefix_multi"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam),
        "--coord",
        c1,
        "--coord",
        c2,
        "--region-type",
        "paraphase",
        "--ref",
        str(ref),
        "--outdir",
        str(outdir),
        "--prefix",
        prefix,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stderr
    assert (outdir / expected).is_file()
    assert (outdir / expected.replace(".html", ".json.gz")).is_file()


def test_deploy_cli_fails_when_outdir_missing(tmp_path):
    missing = tmp_path / "does_not_exist_deploy_dir"
    assert not missing.exists()
    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "deploy",
        "--outdir",
        str(missing),
        "--port",
        "8765",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    assert result.returncode == 1  # run_deploy sys.exit(1) for missing dir
    combined = result.stderr + result.stdout
    assert "Directory not found" in combined or "Error" in combined


def test_plot_second_coord_invalid_fails_first_valid(tmp_path):
    root = Path("tests/data/inputs")
    bam = root / "paraviewer/HG002.paraphase.small.bam"
    ref = root / "hg38_small.fa.gz"
    bed = root / "paraviewer/HG002_regions.bed"
    c1 = _bed_coord(bed, "chr16")
    outdir = tmp_path / "out_bad_second"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam),
        "--coord",
        c1,
        "--coord",
        "chr1-badcoord",
        "--region-type",
        "paraphase",
        "--ref",
        str(ref),
        "--outdir",
        str(outdir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert result.returncode == 1  # ValueError from parse_coordinate
    assert "Invalid coordinate format" in result.stderr
    # No plot artifacts from a failed run (only empty outdir we created).
    assert list(outdir.glob("*_bokeh.html")) == []
