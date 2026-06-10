import subprocess
import sys
from pathlib import Path


def _run_orographer_plot(
    extra_args: list[str],
    outdir: str,
    region_type: str = "paraphase",
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
        region_type,
        "--ref",
        str(ref),
        "--outdir",
        outdir,
    ]

    cmd = base_args + extra_args
    return subprocess.run(cmd, capture_output=True, text=True)


def _run_orographer_plot_help() -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "orographer", "plot", "--help"],
        capture_output=True,
        text=True,
    )


def test_prefix_validation_rejects_hyphen(tmp_path: Path) -> None:
    # --prefix must be alphanumeric with underscores only.
    # argparse parser.error() → exit code 2 (before BAM I/O).
    outdir = str(tmp_path / "out_prefix")
    result = _run_orographer_plot(["--prefix", "bad-prefix"], outdir=outdir)
    assert result.returncode == 2
    assert "--prefix must be alphanumeric" in result.stderr


def test_other_bam_max_two_errors(tmp_path: Path) -> None:
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


def test_other_vcf_max_two_errors(tmp_path: Path) -> None:
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


def test_other_sample_label_length_mismatch_errors(tmp_path: Path) -> None:
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


def test_invalid_coord_format_errors(tmp_path: Path) -> None:
    # parse_coordinate raises ValueError after argparse succeeds → uncaught → exit 1.
    outdir = str(tmp_path / "out_invalid_coord")
    result = _run_orographer_plot(["--coord", "chr16-171801-175500"], outdir=outdir)
    assert result.returncode == 1
    assert "Invalid coordinate format" in result.stderr


def test_plot_help_hides_unreleased_region_type_options() -> None:
    result = _run_orographer_plot_help()

    assert result.returncode == 0
    assert "complex_sv" not in result.stdout.lower()
    assert "isoseq" not in result.stdout.lower()
    assert "--region-connection-threshold" not in result.stdout
    assert "--no-expansion" not in result.stdout
    assert "--force-isoseq" not in result.stdout
    assert "--enable-experimental-isoseq" not in result.stdout
    assert "--enable-experimental-region-types" not in result.stdout


def test_complex_sv_region_type_is_hidden_by_default(tmp_path: Path) -> None:
    result = _run_orographer_plot(
        [],
        str(tmp_path / "out_complex_sv"),
        region_type="complex_sv",
    )

    assert result.returncode == 2
    assert "--region-type must be one of: paraphase" in result.stderr


def test_isoseq_region_type_is_hidden_by_default(tmp_path: Path) -> None:
    result = _run_orographer_plot([], str(tmp_path / "out_isoseq"), region_type="isoseq")

    assert result.returncode == 2
    assert "--region-type must be one of: paraphase" in result.stderr


def test_isoseq_requires_gtf_when_experimental_region_types_enabled(tmp_path: Path) -> None:
    result = _run_orographer_plot(
        ["--enable-experimental-region-types"],
        str(tmp_path / "out_isoseq"),
        region_type="isoseq",
    )

    assert result.returncode == 2
    assert "--gtf is required when --region-type isoseq" in result.stderr


def test_force_isoseq_requires_experimental_gate(tmp_path: Path) -> None:
    result = _run_orographer_plot(["--force-isoseq"], str(tmp_path / "out_force_isoseq"))

    assert result.returncode == 2
    assert "--force-isoseq requires --enable-experimental-region-types" in result.stderr


def test_isoseq_hidden_gate_accepts_force_isoseq_arg(tmp_path: Path) -> None:
    result = _run_orographer_plot(
        ["--enable-experimental-region-types", "--force-isoseq"],
        str(tmp_path / "out_isoseq_force"),
        region_type="isoseq",
    )

    assert result.returncode == 2
    assert "--gtf is required when --region-type isoseq" in result.stderr
    assert "unrecognized arguments: --force-isoseq" not in result.stderr
