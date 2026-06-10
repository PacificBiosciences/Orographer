from __future__ import annotations

import logging
import sys
from types import SimpleNamespace

import pytest

import orographer.__main__ as cli
from orographer.utils import OutputConfig


def parse_args(monkeypatch: pytest.MonkeyPatch, argv: list[str]):
    monkeypatch.setattr(sys, "argv", ["orographer", *argv])
    return cli.setup_args()


def test_setup_args_parses_plot_command_with_multi_bam_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    args = parse_args(
        monkeypatch,
        [
            "plot",
            "--enable-experimental-region-types",
            "--bam",
            "primary.bam",
            "--coord",
            "chr1:100-200",
            "--coord",
            "chr2:300-400",
            "--region-type",
            "complex_sv",
            "--ref",
            "ref.fa",
            "--outdir",
            "out",
            "--prefix",
            "sample_1",
            "--gtf",
            "genes.gtf",
            "--vcf",
            "primary.vcf",
            "--other-bam",
            "other1.bam",
            "--other-vcf",
            "other1.vcf",
            "--sample-label",
            "primary",
            "--other-sample-label",
            "other1",
            "--no-dotplot",
            "--region-connection-threshold",
            "3",
            "--no-expansion",
            "--verbose",
        ],
    )

    assert args.command == "plot"
    assert args.coord == ["chr1:100-200", "chr2:300-400"]
    assert args.other_bam == ["other1.bam"]
    assert args.other_vcf == ["other1.vcf"]
    assert args.other_sample_label == ["other1"]
    assert args.no_dotplot is True
    assert args.region_connection_threshold == 3
    assert args.no_expansion is True
    assert args.verbose is True


def test_setup_args_rejects_region_connection_threshold_below_one(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    with pytest.raises(SystemExit) as exc_info:
        parse_args(
            monkeypatch,
            [
                "plot",
                "--enable-experimental-region-types",
                "--bam",
                "input.bam",
                "--coord",
                "chr1:100-200",
                "--region-type",
                "complex_sv",
                "--ref",
                "ref.fa",
                "--outdir",
                "out",
                "--region-connection-threshold",
                "0",
            ],
        )

    assert exc_info.value.code == 2


def test_setup_args_rejects_bad_prefix(monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(SystemExit) as exc_info:
        parse_args(
            monkeypatch,
            [
                "plot",
                "--bam",
                "input.bam",
                "--coord",
                "chr1:100-200",
                "--region-type",
                "paraphase",
                "--ref",
                "ref.fa",
                "--outdir",
                "out",
                "--prefix",
                "bad-prefix",
            ],
        )

    assert exc_info.value.code == 2


def test_setup_args_rejects_other_sample_label_count_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    with pytest.raises(SystemExit) as exc_info:
        parse_args(
            monkeypatch,
            [
                "plot",
                "--bam",
                "input.bam",
                "--coord",
                "chr1:100-200",
                "--region-type",
                "paraphase",
                "--ref",
                "ref.fa",
                "--outdir",
                "out",
                "--other-bam",
                "other1.bam",
                "--other-bam",
                "other2.bam",
                "--other-sample-label",
                "only-one",
            ],
        )

    assert exc_info.value.code == 2


def test_setup_args_parses_deploy_command(monkeypatch: pytest.MonkeyPatch) -> None:
    args = parse_args(monkeypatch, ["deploy", "--outdir", "out", "--port", "8765"])

    assert args.command == "deploy"
    assert args.outdir == "out"
    assert args.port == 8765


def test_setup_logging_sets_requested_levels() -> None:
    logger = cli.setup_logging(verbose=True)

    assert logger.name == cli.__name__
    assert logging.getLogger("numexpr").level == logging.WARNING
    assert logging.getLogger("bokeh.io.state").level == logging.WARNING


def test_run_plot_command_validates_and_dispatches(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict = {}
    validations: list[str] = []

    monkeypatch.setattr(cli, "setup_logging", lambda _verbose: logging.getLogger("test"))
    monkeypatch.setattr(cli, "validate_bam_file", lambda path: validations.append(path))

    def run_orographer(*args, **kwargs) -> None:
        captured["args"] = args
        captured["kwargs"] = kwargs

    monkeypatch.setattr(cli, "orographer", run_orographer)

    args = SimpleNamespace(
        verbose=False,
        bam="primary.bam",
        other_bam=["other.bam"],
        other_vcf=["other.vcf"],
        other_sample_label=["other"],
        coord=["chr1:1_000-2_000"],
        outdir="out",
        prefix="sample",
        region_type="paraphase",
        ref="ref.fa",
        gtf="genes.gtf",
        vcf="primary.vcf",
        sample_label="primary",
        no_dotplot=True,
        region_connection_threshold=4,
        no_expansion=True,
    )

    cli.run_plot_command(args)

    assert validations == ["primary.bam", "other.bam"]
    assert captured["args"] == (
        "paraphase",
        "primary.bam",
        ["chr1:1000-2000"],
        "ref.fa",
        OutputConfig("out", "sample"),
        "genes.gtf",
        "primary.vcf",
        ["other.bam"],
        ["other.vcf"],
        "primary",
        ["other"],
    )
    assert captured["kwargs"] == {
        "show_dotplot": False,
        "region_connection_threshold": 4,
        "no_expansion": True,
        "force_isoseq": False,
    }


def test_run_deploy_command_dispatches(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict = {}

    monkeypatch.setattr(
        cli,
        "run_deploy",
        lambda outdir, port: captured.update({"outdir": outdir, "port": port}),
    )

    cli.run_deploy_command(SimpleNamespace(outdir="out", port=9000))

    assert captured == {"outdir": "out", "port": 9000}


def test_main_dispatches_plot_command(monkeypatch: pytest.MonkeyPatch, capsys) -> None:
    calls: list[object] = []
    args = SimpleNamespace(command="plot")

    monkeypatch.setattr(cli, "setup_args", lambda: args)
    monkeypatch.setattr(cli, "run_plot_command", lambda parsed_args: calls.append(parsed_args))

    cli.main()

    assert calls == [args]
    assert "Orographer v" in capsys.readouterr().err


def test_main_dispatches_deploy_command(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[object] = []
    args = SimpleNamespace(command="deploy")

    monkeypatch.setattr(cli, "setup_args", lambda: args)
    monkeypatch.setattr(cli, "run_deploy_command", lambda parsed_args: calls.append(parsed_args))

    cli.main()

    assert calls == [args]


def test_main_rejects_unknown_command(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli, "setup_args", lambda: SimpleNamespace(command="unknown"))

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 1
