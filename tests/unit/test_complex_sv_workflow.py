from __future__ import annotations

from collections.abc import Callable

import pytest

from orographer import complex_sv
from orographer.utils import Region


class FakeBam:
    def __init__(self) -> None:
        self.closed = False

    def close(self) -> None:
        self.closed = True


def install_common_workflow_mocks(
    monkeypatch: pytest.MonkeyPatch,
    *,
    discovered_reads: set[str],
    final_regions: list[Region] | None = None,
    render_spy: Callable[[str, Region], tuple[dict, dict, dict]] | None = None,
) -> dict:
    captured: dict = {"closed_bams": []}

    def alignment_file(path: str, mode: str) -> FakeBam:
        assert mode == "rb"
        captured["opened_bam"] = path
        bam = FakeBam()
        captured["bam"] = bam
        return bam

    def discover_connected_reads(
        _bam,
        initial_frontier,
        mapq_threshold: int,
        region_connection_threshold: int,
    ):
        captured["initial_frontier"] = initial_frontier
        captured["discovery_mapq"] = mapq_threshold
        captured["region_connection_threshold"] = region_connection_threshold
        return discovered_reads, {"read-a": []}, [("chr1", 1, 10)]

    def build_final_regions(_summaries_by_read, _visited_intervals):
        captured["build_final_regions_called"] = True
        return final_regions or [Region("chr7", 100, 200, "")]

    def discover_reads_in_regions(_bam, regions, mapq_threshold: int):
        captured["explicit_regions"] = regions
        captured["explicit_discovery_mapq"] = mapq_threshold
        return discovered_reads, {"read-a": []}, regions

    def parse_annotation_file(gtf_file: str | None, region: Region) -> list[str]:
        return [f"gene:{gtf_file}:{region.chromosome}"]

    def parse_vcf_file(vcf_file: str | None, region: Region) -> list[str]:
        return [f"vcf:{vcf_file}:{region.chromosome}"]

    def render_region_for_bam(
        bam_path: str,
        region: Region,
        _discovered_reads: set[str],
        _reference_path: str | None,
        _haplotypes_to_track: set[int] | None,
        mapq_threshold: int,
    ) -> tuple[dict, dict, dict]:
        captured.setdefault("render_calls", []).append((bam_path, region, mapq_threshold))
        if render_spy:
            return render_spy(bam_path, region)
        return (
            {f"segments:{bam_path}": []},
            {-1: ([region.start], [1])},
            {0: [{"pos": region.start}]},
        )

    monkeypatch.setattr(complex_sv.pysam, "AlignmentFile", alignment_file)
    monkeypatch.setattr(complex_sv, "discover_connected_reads", discover_connected_reads)
    monkeypatch.setattr(complex_sv, "discover_reads_in_regions", discover_reads_in_regions)
    monkeypatch.setattr(complex_sv, "_build_final_regions", build_final_regions)
    monkeypatch.setattr(complex_sv, "parse_annotation_file", parse_annotation_file)
    monkeypatch.setattr(complex_sv, "parse_vcf_file", parse_vcf_file)
    monkeypatch.setattr(complex_sv, "render_region_for_bam", render_region_for_bam)
    return captured


def test_process_complex_sv_falls_back_to_seed_regions_when_no_reads(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured = install_common_workflow_mocks(monkeypatch, discovered_reads=set())

    result = complex_sv.process_complex_sv(
        coordinate_strs=["chr1:100-200"],
        bam_files=["primary.bam"],
        vcf_files=[None],
        reference_path="ref.fa",
        gtf_file="genes.gtf",
        sample_labels=["primary"],
        region_type="complex_sv",
    )

    assert captured["opened_bam"] == "primary.bam"
    assert captured["bam"].closed
    assert "build_final_regions_called" not in captured
    assert result[0]["region"] == Region("chr1", 100, 200, "chr1:100-200")
    assert result[0]["gene_annotations"] == ["gene:genes.gtf:chr1"]
    assert result[0]["bam_rows"][0]["sample_label"] == "primary"
    assert result[0]["bam_rows"][0]["vcf_variants"] == ["vcf:None:chr1"]


def test_process_complex_sv_no_expansion_uses_explicit_seed_regions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured = install_common_workflow_mocks(
        monkeypatch,
        discovered_reads={"read-a"},
        final_regions=[Region("chr7", 100, 200, "")],
    )

    result = complex_sv.process_complex_sv(
        coordinate_strs=["chr1:100-200", "chr2:300-400"],
        bam_files=["primary.bam"],
        vcf_files=[None],
        reference_path="ref.fa",
        gtf_file="genes.gtf",
        sample_labels=["primary"],
        region_type="complex_sv",
        no_expansion=True,
    )

    assert captured["opened_bam"] == "primary.bam"
    assert captured["bam"].closed
    assert captured["explicit_regions"] == [("chr1", 100, 200), ("chr2", 300, 400)]
    assert captured["explicit_discovery_mapq"] == 0
    assert "initial_frontier" not in captured
    assert "build_final_regions_called" not in captured
    assert [region_data["region"] for region_data in result] == [
        Region("chr1", 100, 200, "chr1:100-200"),
        Region("chr2", 300, 400, "chr2:300-400"),
    ]


def test_process_complex_sv_uses_final_regions_and_orders_bam_rows(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    final_region = Region("chr7", 100, 200, "")
    captured = install_common_workflow_mocks(
        monkeypatch,
        discovered_reads={"read-a"},
        final_regions=[final_region],
    )

    result = complex_sv.process_complex_sv(
        coordinate_strs=["chr1:100-200", "chr1:180-250"],
        bam_files=["primary.bam", "parent1.bam", "parent2.bam"],
        vcf_files=["primary.vcf", "parent1.vcf", "parent2.vcf"],
        reference_path="ref.fa",
        gtf_file="genes.gtf",
        sample_labels=["primary", "parent1", "parent2"],
        region_type="complex_sv",
        region_connection_threshold=3,
    )

    assert captured["build_final_regions_called"]
    assert captured["discovery_mapq"] == 0
    assert captured["region_connection_threshold"] == 3
    assert len(captured["initial_frontier"]) == 1
    assert [call[0] for call in captured["render_calls"]] == [
        "primary.bam",
        "parent1.bam",
        "parent2.bam",
    ]
    bam_rows = result[0]["bam_rows"]
    assert [row["sample_label"] for row in bam_rows] == ["parent1", "parent2", "primary"]
    assert [row["sample_index"] for row in bam_rows] == [1, 2, 0]
    assert bam_rows[0]["vcf_variants"] == ["vcf:parent1.vcf:chr7"]


def test_process_complex_sv_truncates_final_regions_when_total_span_is_too_large(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    final_regions = [
        Region("chr1", 1, 100, ""),
        Region("chr2", 1, 100, ""),
    ]
    captured = install_common_workflow_mocks(
        monkeypatch,
        discovered_reads={"read-a"},
        final_regions=final_regions,
    )

    def truncate_regions(regions: list[Region], max_bp: int) -> list[Region]:
        captured["truncated_regions"] = regions
        captured["truncate_max_bp"] = max_bp
        return [regions[0]]

    monkeypatch.setattr(complex_sv, "MAX_TOTAL_RENDER_BP", 150)
    monkeypatch.setattr(complex_sv, "_truncate_regions", truncate_regions)

    result = complex_sv.process_complex_sv(
        coordinate_strs=["chr1:1-10"],
        bam_files=["primary.bam"],
        vcf_files=[],
        reference_path=None,
        gtf_file=None,
        sample_labels=[],
        region_type="complex_sv",
    )

    assert captured["truncated_regions"] == final_regions
    assert captured["truncate_max_bp"] == 150
    assert [region_data["region"] for region_data in result] == [final_regions[0]]
    assert result[0]["bam_rows"][0]["sample_label"] is None


def test_process_complex_sv_wraps_primary_bam_open_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def alignment_file(_path: str, _mode: str):
        raise OSError("nope")

    monkeypatch.setattr(complex_sv.pysam, "AlignmentFile", alignment_file)

    with pytest.raises(RuntimeError, match=r"Cannot open primary BAM: missing\.bam"):
        complex_sv.process_complex_sv(
            coordinate_strs=["chr1:1-10"],
            bam_files=["missing.bam"],
            vcf_files=[],
            reference_path=None,
            gtf_file=None,
            sample_labels=[],
            region_type="complex_sv",
        )
