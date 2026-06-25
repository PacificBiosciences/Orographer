from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pytest

from orographer.dotplot import MAX_DOTPLOT_REGION_BP
from orographer.utils import (
    COMPLEX_SV_REGION_TYPE,
    PARAPHASE_REGION_TYPE,
    OutputConfig,
    Region,
)

workflow = importlib.import_module("orographer.orographer")


def test_attach_dotplot_images_adds_combined_payload_when_reference_sequence_is_available(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first_region = Region("chr1", 100, 120, "chr1:100-120")
    second_region = Region("chr2", 200, 230, "chr2:200-230")
    region_data = [{"region": first_region}, {"region": second_region}]
    expected = np.eye(2, dtype=np.float32)
    fetched_regions: list[Region] = []

    def fetch_reference(_ref: str, region: Region) -> str:
        fetched_regions.append(region)
        return "AAAA" if region == first_region else "CCCC"

    captured_sequence = ""

    def compute_identity(sequence: str) -> np.ndarray:
        nonlocal captured_sequence
        captured_sequence = sequence
        return expected

    monkeypatch.setattr(workflow, "fetch_reference_sequence", fetch_reference)
    monkeypatch.setattr(workflow, "compute_self_identity", compute_identity)

    workflow._attach_dotplot_images(region_data, "ref.fa")

    payload = region_data[0]["dotplot_payload"]
    assert payload["matrix"] is expected
    # total_span uses actual sequence lengths (4+4=8), not nominal region spans (21+31=52).
    assert payload["total_span"] == 8
    assert captured_sequence == "AAAACCCC"
    assert fetched_regions == [first_region, second_region]
    assert region_data[1]["dotplot_payload"] is None
    assert payload["blocks"] == [
        {
            "label": "chr1:0.0-0.0 Mb",
            "coordinate_label": "chr1:100-120",
            "chromosome": "chr1",
            "start": 100,
            "end": 120,
            "offset_start": 0,
            "offset_end": 4,   # actual sequence length, not nominal span 21
            "span": 21,
            "sequence_length": 4,
        },
        {
            "label": "chr2:0.0-0.0 Mb",
            "coordinate_label": "chr2:200-230",
            "chromosome": "chr2",
            "start": 200,
            "end": 230,
            "offset_start": 4,   # immediately after first actual sequence
            "offset_end": 8,     # actual cumulative length, not nominal 52
            "span": 31,
            "sequence_length": 4,
        },
    ]


def test_attach_dotplot_images_block_offsets_use_actual_sequence_lengths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Block offsets must reflect the actual fetched sequence length, not the nominal span.

    When a region extends past the chromosome end, pysam returns fewer bases than
    region.end - region.start + 1.  The matrix is built from those shorter sequences,
    so the extent and pink-line coordinates must use the actual lengths; otherwise the
    pink boundary line is drawn far to the right of the actual chr14/chr20 junction.
    """
    # chr14 region nominally 1,740,001 bp but only 1,493,719 bases returned (past chrom end).
    # chr20 region nominally and actually 1,300,001 bp.
    first_region = Region("chr14", 105_550_000, 107_290_000, "chr14:105550000-107290000")
    second_region = Region("chr20", 39_500_000, 40_800_000, "chr20:39500000-40800000")
    region_data = [{"region": first_region}, {"region": second_region}]

    chr14_actual = "A" * 1_493_719   # truncated at chromosome end
    chr20_actual = "C" * 1_300_001

    def fetch_reference(_ref: str, region: Region) -> str:
        return chr14_actual if region is first_region else chr20_actual

    monkeypatch.setattr(workflow, "fetch_reference_sequence", fetch_reference)
    monkeypatch.setattr(workflow, "compute_self_identity", lambda _s: np.eye(2, dtype=np.float32))

    workflow._attach_dotplot_images(region_data, "ref.fa")

    payload = region_data[0]["dotplot_payload"]
    blocks = payload["blocks"]

    # Block offsets must be built from actual sequence lengths, not nominal spans.
    assert blocks[0]["offset_start"] == 0
    assert blocks[0]["offset_end"] == 1_493_719, (
        f"chr14 block offset_end should be actual seq length 1,493,719, "
        f"got {blocks[0]['offset_end']:,}"
    )
    assert blocks[1]["offset_start"] == 1_493_719
    assert blocks[1]["offset_end"] == 1_493_719 + 1_300_001

    # total_span in the payload must also reflect actual sequence total.
    assert payload["total_span"] == 1_493_719 + 1_300_001


def test_attach_dotplot_images_includes_individual_payloads_for_multi_region(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first_region = Region("chr1", 100, 120, "chr1:100-120")
    second_region = Region("chr2", 200, 230, "chr2:200-230")
    region_data = [{"region": first_region}, {"region": second_region}]
    matrices_by_seq: dict[str, np.ndarray] = {}

    def fetch_reference(_ref: str, region: Region) -> str:
        return "AAAA" if region == first_region else "CCCC"

    def compute_identity(sequence: str) -> np.ndarray:
        m = np.eye(2, dtype=np.float32) * len(matrices_by_seq)
        matrices_by_seq[sequence] = m
        return m

    monkeypatch.setattr(workflow, "fetch_reference_sequence", fetch_reference)
    monkeypatch.setattr(workflow, "compute_self_identity", compute_identity)

    workflow._attach_dotplot_images(region_data, "ref.fa")

    payload = region_data[0]["dotplot_payload"]
    assert "individual_payloads" in payload
    indivs = payload["individual_payloads"]
    assert len(indivs) == 2
    assert indivs[0]["region"] is first_region
    assert indivs[1]["region"] is second_region
    assert indivs[0]["matrix"] is matrices_by_seq["AAAA"]
    assert indivs[1]["matrix"] is matrices_by_seq["CCCC"]
    assert indivs[0]["label"] == "chr1:100-120"
    assert indivs[1]["label"] == "chr2:200-230"


def test_attach_dotplot_images_skips_individual_payloads_for_single_region(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    region = Region("chr1", 100, 120, "chr1:100-120")
    region_data = [{"region": region}]

    monkeypatch.setattr(workflow, "fetch_reference_sequence", lambda _r, _reg: "AAAA")
    monkeypatch.setattr(workflow, "compute_self_identity", lambda _s: np.eye(2, dtype=np.float32))

    workflow._attach_dotplot_images(region_data, "ref.fa")

    payload = region_data[0]["dotplot_payload"]
    assert payload["individual_payloads"] == []


def test_attach_dotplot_images_adds_repeat_density_to_each_region(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first_region = Region("chr1", 100, 200, "chr1:100-200")
    second_region = Region("chr2", 300, 400, "chr2:300-400")
    region_data = [{"region": first_region}, {"region": second_region}]

    monkeypatch.setattr(workflow, "fetch_reference_sequence", lambda _r, _reg: "ACGT" * 25)

    workflow._attach_dotplot_images(region_data, "ref.fa")

    assert "repeat_density" in region_data[0]
    assert "repeat_density" in region_data[1]
    assert isinstance(region_data[0]["repeat_density"], np.ndarray)
    assert isinstance(region_data[1]["repeat_density"], np.ndarray)
    assert region_data[0]["repeat_density"].shape == (512,)
    assert region_data[1]["repeat_density"].shape == (512,)


def test_attach_dotplot_images_single_region_label_is_coordinate_string(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    region = Region("chr1", 100, 120, "chr1:100-120")
    region_data = [{"region": region}]

    monkeypatch.setattr(workflow, "fetch_reference_sequence", lambda _r, _reg: "AAAA")
    monkeypatch.setattr(workflow, "compute_self_identity", lambda _s: np.eye(2, dtype=np.float32))

    workflow._attach_dotplot_images(region_data, "ref.fa")

    payload = region_data[0]["dotplot_payload"]
    assert "Regions are concatenated" not in payload["label"]
    assert payload["label"] == "chr1:100-120"


def test_attach_dotplot_images_adds_repeat_density_for_single_region(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    region = Region("chr1", 100, 200, "chr1:100-200")
    region_data = [{"region": region}]

    monkeypatch.setattr(workflow, "fetch_reference_sequence", lambda _r, _reg: "ACGT" * 25)

    workflow._attach_dotplot_images(region_data, "ref.fa")

    assert "repeat_density" in region_data[0]
    assert isinstance(region_data[0]["repeat_density"], np.ndarray)
    assert region_data[0]["repeat_density"].shape == (512,)


def test_attach_dotplot_images_skips_oversized_and_missing_reference(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    oversized_region = Region("chr1", 1, MAX_DOTPLOT_REGION_BP + 2, "chr1:1-large")
    missing_region = Region("chr2", 10, 20, "chr2:10-20")
    oversized_region_data = [{"region": oversized_region}]
    missing_region_data = [{"region": missing_region}]
    fetch_calls: list[Region] = []

    def fetch_reference(_reference_path: str, region: Region) -> None:
        fetch_calls.append(region)
        return None

    monkeypatch.setattr(workflow, "fetch_reference_sequence", fetch_reference)
    monkeypatch.setattr(
        workflow,
        "compute_self_identity",
        lambda _seq: pytest.fail("dotplot should not be computed"),
    )

    workflow._attach_dotplot_images(oversized_region_data, "ref.fa")

    assert oversized_region_data[0]["dotplot_payload"] is None
    assert fetch_calls == []

    workflow._attach_dotplot_images(missing_region_data, "ref.fa")

    assert missing_region_data[0]["dotplot_payload"] is None
    assert fetch_calls == [missing_region]


def test_orographer_paraphase_orders_other_bams_above_primary(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: dict = {}

    def process_paraphase(paths, region: Region) -> dict:
        return {
            "segments_by_read": {},
            "gene_annotations": [f"genes-{paths.bam_file}"],
            "vcf_variants": [],
            "coverage_tracks": {},
            "bam_path": paths.bam_file,
        }

    def plot_reads_bokeh(region_data_list: list[dict], output_config: OutputConfig) -> str:
        captured["region_data_list"] = region_data_list
        captured["output_config"] = output_config
        return "plot.html"

    monkeypatch.setattr(workflow, "validate_bam_file", lambda _path: None)
    monkeypatch.setattr(workflow, "validate_output_dir", lambda _path: None)
    monkeypatch.setattr(workflow, "process_paraphase", process_paraphase)
    monkeypatch.setattr(workflow, "_attach_dotplot_images", lambda _regions, _ref: None)
    monkeypatch.setattr(workflow, "plot_reads_bokeh", plot_reads_bokeh)

    workflow.orographer(
        region_type=PARAPHASE_REGION_TYPE,
        bam_file="primary.bam",
        coordinate_strs="chr1:10-20",
        reference_path="ref.fa",
        output_config=OutputConfig(str(tmp_path), "sample"),
        gtf_file="genes.gtf",
        vcf_file="primary.vcf",
        other_bam_files=["parent1.bam", "parent2.bam"],
        other_vcf_files=["parent1.vcf", "parent2.vcf"],
        sample_label="proband",
        other_sample_labels=["parent1", "parent2"],
    )

    region_data = captured["region_data_list"][0]
    bam_rows = region_data["bam_rows"]
    assert [row["bam_path"] for row in bam_rows] == [
        "parent1.bam",
        "parent2.bam",
        "primary.bam",
    ]
    assert [row["sample_label"] for row in bam_rows] == ["parent1", "parent2", "proband"]
    assert region_data["gene_annotations"] == ["genes-primary.bam"]


def test_orographer_can_skip_dotplot_attachment(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    dotplot_called = False

    def mark_dotplot_called(_regions: list[dict], _reference_path: str) -> None:
        nonlocal dotplot_called
        dotplot_called = True

    monkeypatch.setattr(workflow, "validate_bam_file", lambda _path: None)
    monkeypatch.setattr(workflow, "validate_output_dir", lambda _path: None)
    monkeypatch.setattr(
        workflow,
        "process_paraphase",
        lambda _paths, _region: {
            "segments_by_read": {},
            "gene_annotations": [],
            "vcf_variants": [],
            "coverage_tracks": {},
        },
    )
    monkeypatch.setattr(workflow, "_attach_dotplot_images", mark_dotplot_called)
    monkeypatch.setattr(workflow, "plot_reads_bokeh", lambda _regions, _config: "plot.html")

    workflow.orographer(
        region_type=PARAPHASE_REGION_TYPE,
        bam_file="primary.bam",
        coordinate_strs=["chr1:10-20"],
        reference_path="ref.fa",
        output_config=OutputConfig(str(tmp_path), None),
        gtf_file="genes.gtf",
        vcf_file=None,
        other_bam_files=[],
        other_vcf_files=[],
        sample_label=None,
        other_sample_labels=[],
        show_dotplot=False,
    )

    assert not dotplot_called


def test_orographer_complex_sv_uses_seed_coordinates_for_output_filename(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: dict = {}
    seed_coordinates = ["chr1:10-20", "chr2:30-40"]

    def process_complex_sv(**kwargs) -> list[dict]:
        captured["process_complex_sv_kwargs"] = kwargs
        return [{"region": Region("chr1", 10, 20, "chr1:10-20")}]

    def plot_reads_bokeh(region_data_list: list[dict], output_config: OutputConfig) -> str:
        captured["region_data_list"] = region_data_list
        captured["output_config"] = output_config
        return "plot.html"

    monkeypatch.setattr(workflow, "validate_bam_file", lambda _path: None)
    monkeypatch.setattr(workflow, "validate_output_dir", lambda _path: None)
    monkeypatch.setattr(workflow, "process_complex_sv", process_complex_sv)
    monkeypatch.setattr(workflow, "_attach_dotplot_images", lambda _regions, _ref: None)
    monkeypatch.setattr(workflow, "plot_reads_bokeh", plot_reads_bokeh)

    workflow.orographer(
        region_type=COMPLEX_SV_REGION_TYPE,
        bam_file="primary.bam",
        coordinate_strs=seed_coordinates,
        reference_path="ref.fa",
        output_config=OutputConfig(str(tmp_path), "sv"),
        gtf_file="genes.gtf",
        vcf_file="primary.vcf",
        other_bam_files=["other.bam"],
        other_vcf_files=[],
        sample_label="primary",
        other_sample_labels=["other"],
        region_connection_threshold=5,
        no_expansion=True,
    )

    output_config = captured["output_config"]
    process_kwargs = captured["process_complex_sv_kwargs"]
    assert output_config.filename_regions == seed_coordinates
    assert process_kwargs["coordinate_strs"] == seed_coordinates
    assert process_kwargs["bam_files"] == ["primary.bam", "other.bam"]
    assert process_kwargs["vcf_files"] == ["primary.vcf", None]
    assert process_kwargs["sample_labels"] == ["primary", "other"]
    assert process_kwargs["region_connection_threshold"] == 5
    assert process_kwargs["no_expansion"] is True
