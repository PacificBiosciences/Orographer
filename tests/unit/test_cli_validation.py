import subprocess
import sys
from pathlib import Path


def _run_orographer_plot(
    extra_args: list[str],
    outdir: str,
) -> subprocess.CompletedProcess[str]:
    # Use existing fixture files so we only fail for the specific validation under test.
    root = Path("tests/data/inputs")
    bam = root / "paraviewer/HG002.paraphase.small.bam"
    ref = root / "hg38_small.fa.gz"

    base_args = [
        sys.executable,
        "-m",
        "orographer",
        "plot",
        "--bam",
        str(bam),
        "--coord",
        "chr16:171801-175500",
        "--region-type",
        "paraphase",
        "--ref",
        str(ref),
        "--outdir",
        outdir,
    ]

    cmd = base_args + extra_args
    return subprocess.run(cmd, capture_output=True, text=True)


def test_prefix_validation_rejects_hyphen(tmp_path):
    # --prefix must be alphanumeric with underscores only.
    # argparse parser.error() → exit code 2 (before BAM I/O).
    outdir = str(tmp_path / "out_prefix")
    result = _run_orographer_plot(["--prefix", "bad-prefix"], outdir=outdir)
    assert result.returncode == 2
    assert "--prefix must be alphanumeric" in result.stderr


def test_other_bam_max_two_errors(tmp_path):
    # setup_args should error before BAM content validation.
    root = Path("tests/data/inputs")
    bam_other = root / "paraviewer/HG002.paraphase.small.bam"
    extra = [
        "--other-bam",
        str(bam_other),
        "--other-bam",
        str(bam_other),
        "--other-bam",
        str(bam_other),
    ]
    outdir = str(tmp_path / "out_other_bam")
    result = _run_orographer_plot(extra, outdir=outdir)
    assert result.returncode == 2
    assert "--other-bam may be specified at most two times" in result.stderr


def test_other_vcf_max_two_errors(tmp_path):
    root = Path("tests/data/inputs")
    vcf_other = root / "paraviewer/tmp.yaml"  # any existing file is fine (we won't reach plot)
    extra = [
        "--vcf",
        str(vcf_other),
        "--other-vcf",
        str(vcf_other),
        "--other-vcf",
        str(vcf_other),
        "--other-vcf",
        str(vcf_other),
        "--other-bam",
        str(root / "paraviewer/HG002.paraphase.small.bam"),
        "--other-bam",
        str(root / "paraviewer/HG002.paraphase.small.bam"),
    ]
    outdir = str(tmp_path / "out_other_vcf")
    result = _run_orographer_plot(extra, outdir=outdir)
    assert result.returncode == 2
    assert "--other-vcf may be specified at most two times" in result.stderr


def test_other_sample_label_length_mismatch_errors(tmp_path):
    root = Path("tests/data/inputs")
    bam_other = root / "paraviewer/HG002.paraphase.small.bam"

    # Two other-bams, but only one other-sample-label.
    extra = [
        "--other-bam",
        str(bam_other),
        "--other-bam",
        str(bam_other),
        "--other-sample-label",
        "only_one_label",
    ]
    outdir = str(tmp_path / "out_other_sample_label")
    result = _run_orographer_plot(extra, outdir=outdir)
    assert result.returncode == 2
    assert "--other-sample-label must be specified once per --other-bam" in result.stderr


def test_invalid_coord_format_errors(tmp_path):
    # parse_coordinate raises ValueError after argparse succeeds → uncaught → exit 1.
    outdir = str(tmp_path / "out_invalid_coord")
    result = _run_orographer_plot(["--coord", "chr16-171801-175500"], outdir=outdir)
    assert result.returncode == 1
    assert "Invalid coordinate format" in result.stderr
