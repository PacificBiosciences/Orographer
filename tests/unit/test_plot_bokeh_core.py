from __future__ import annotations

import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from bokeh.models import ColumnDataSource, CustomJS, Div, GlyphRenderer, HoverTool, Range1d, TapTool
from bokeh.plotting import figure

from orographer.plot_bokeh import plot_bokeh
from orographer.plot_bokeh.callbacks import compact_bokeh_json, expand_compact_bokeh_json
from orographer.plot_bokeh.segment_processing import chevron_tip_fraction
from orographer.utils import (
    COMPLEX_SV_REGION_TYPE,
    ISOSEQ_REGION_TYPE,
    PARAPHASE_REGION_TYPE,
    OutputConfig,
    Region,
)


def segment(
    *,
    pos: int = 100,
    end: int = 200,
    is_fwd_strand: bool = True,
    haplotype_tag: int | None = 1,
    readname: str = "readA",
    chrom: str = "chr1",
    color_tag: str | None = None,
    alignment_order: int = 1,
    fwd_read_start: int = 0,
    fwd_read_end: int = 100,
    phaseset_tag: int | str | None = None,
    inclusion_reason: str | None = None,
    mismatches: list[tuple[int, str]] | None = None,
    insertions: list[tuple[int, str]] | None = None,
    deletions: list[tuple[int, int]] | None = None,
    aligned_blocks: list[tuple[int, int]] | None = None,
):
    return SimpleNamespace(
        pos=pos,
        end=end,
        is_fwd_strand=is_fwd_strand,
        haplotype_tag=haplotype_tag,
        readname=readname,
        chrom=chrom,
        color_tag=color_tag,
        alignment_order=alignment_order,
        fwd_read_start=fwd_read_start,
        fwd_read_end=fwd_read_end,
        phaseset_tag=phaseset_tag,
        inclusion_reason=inclusion_reason,
        mismatches=mismatches or [],
        insertions=insertions or [],
        deletions=deletions or [],
        aligned_blocks=aligned_blocks or [(pos, end)],
    )


def test_compact_bokeh_json_round_trips_common_bokeh_columns() -> None:
    payload = {
        "read_name": ["readA"] * 20 + ["readB"] * 20,
        "layout_read_name": ["readA", "readB"] * 20,
        "x": list(range(100, 260, 4)),
        "y": [0] * 12 + [1] * 14 + [2] * 14,
        "xs": [list(range(100, 260, 4))] * 8,
    }

    compacted = compact_bokeh_json(payload)

    assert expand_compact_bokeh_json(compacted) == payload
    assert "__orog_dict__" in json.dumps(compacted)
    assert "__orog_range__" in json.dumps(compacted)
    assert "__orog_repeat__" in json.dumps(compacted)


def test_get_segment_color_prefers_valid_yc_tag() -> None:
    seg = segment(color_tag="1,128,255", haplotype_tag=None)

    assert plot_bokeh.get_segment_color(seg, "complex_sv") == "#0180ff"


def test_get_segment_color_handles_region_and_haplotype_defaults() -> None:
    invalid_yc = segment(color_tag="bad", haplotype_tag=1, is_fwd_strand=False)
    unassigned = segment(haplotype_tag=None)

    assert plot_bokeh.get_segment_color(invalid_yc, "paraphase") == "grey"
    assert plot_bokeh.get_segment_color(invalid_yc, "complex_sv") == "blue"
    assert plot_bokeh.get_segment_color(unassigned, "complex_sv") == "lightgrey"


def test_isoseq_lazy_chunks_write_static_read_page_json(tmp_path: Path) -> None:
    groups = [
        {
            "group_id": "T1",
            "read_names": ["readA"],
            "assigned_read_count": 1,
        },
        {
            "group_id": "UNASSIGNED",
            "read_names": ["readB"],
            "assigned_read_count": 1,
        },
    ]
    segments_by_read = {
        "readA": [
            segment(
                pos=100,
                end=200,
                readname="readA",
                fwd_read_start=0,
                fwd_read_end=100,
            )
        ],
        "readB": [
            segment(
                pos=120,
                end=180,
                readname="readB",
                fwd_read_start=0,
                fwd_read_end=60,
            )
        ],
    }
    layout = {
        "read_metadata": {
            "readA": {
                "gene_id": "G1",
                "gene_name": "Gene1",
                "transcript_id": "T1",
                "group_id": "T1",
            },
            "readB": {
                "gene_id": "",
                "gene_name": "Unassigned",
                "transcript_id": "UNASSIGNED",
                "group_id": "UNASSIGNED",
            },
        },
        "selected_read_y_start": 1.0,
    }

    total_coverage_data = plot_bokeh._prepare_isoseq_lazy_chunks(
        groups,
        segments_by_read,
        layout,
        90,
        210,
        0,
        0,
        "sample1",
        "plot-class",
        "plot-id",
        chunk_dir=str(tmp_path),
        chunk_url_prefix="chunks",
    )

    assert groups[0]["chunk_url"] == "chunks/r0_row0_reads_manifest.json.gz"
    assert groups[0]["coverage_url"] == "chunks/g0_T1_coverage.json.gz"
    assert groups[1]["chunk_url"] == "chunks/r0_row0_reads_manifest.json.gz"
    assert groups[1]["coverage_url"] == "chunks/g1_UNASSIGNED_coverage.json.gz"
    with gzip.open(tmp_path / "r0_row0_reads_manifest.json.gz", "rt", encoding="utf-8") as handle:
        manifest = json.load(handle)
    with gzip.open(tmp_path / "g0_T1_coverage.json.gz", "rt", encoding="utf-8") as handle:
        coverage_payload = json.load(handle)
    with gzip.open(
        tmp_path / "g1_UNASSIGNED_coverage.json.gz",
        "rt",
        encoding="utf-8",
    ) as handle:
        unassigned_coverage = json.load(handle)
    with gzip.open(tmp_path / "r0_row0_reads_shard0.msgpack.gz", "rb") as handle:
        shard_prefix = handle.read(1)
    assert "bam_path" not in manifest
    assert "coverage_data" not in manifest
    assert manifest["schema"] == "isoseq_read_manifest_v1"
    assert manifest["page_size"] == 100
    assert manifest["shard_size"] == plot_bokeh.ISOSEQ_READ_SHARD_SIZE
    assert manifest["groups"]["g0_T1"]["read_ids"] == [0]
    assert manifest["groups"]["g1_UNASSIGNED"]["read_ids"] == [1]
    assert manifest["shards"] == ["r0_row0_reads_shard0.msgpack.gz"]
    assert shard_prefix
    assert coverage_payload["coverage_data"]["x"]
    assert unassigned_coverage["coverage_data"]["x"]
    assert max(total_coverage_data["y"]) == 2


def test_isoseq_lazy_chunks_accept_custom_read_page_size(tmp_path: Path) -> None:
    groups = [
        {
            "group_id": "T1",
            "read_names": ["readA"],
            "assigned_read_count": 1,
        }
    ]
    segments_by_read = {
        "readA": [
            segment(
                pos=100,
                end=200,
                readname="readA",
                fwd_read_start=0,
                fwd_read_end=100,
            )
        ],
    }
    layout = {
        "read_metadata": {
            "readA": {
                "gene_id": "G1",
                "gene_name": "Gene1",
                "transcript_id": "T1",
                "group_id": "T1",
            },
        },
        "selected_read_y_start": 1.0,
    }

    plot_bokeh._prepare_isoseq_lazy_chunks(
        groups,
        segments_by_read,
        layout,
        90,
        210,
        0,
        0,
        "sample1",
        "plot-class",
        "plot-id",
        chunk_dir=str(tmp_path),
        chunk_url_prefix="chunks",
        read_page_size=50,
    )

    with gzip.open(tmp_path / "r0_row0_reads_manifest.json.gz", "rt", encoding="utf-8") as handle:
        manifest = json.load(handle)
    assert manifest["page_size"] == 50


def test_isoseq_coverage_uses_aligned_blocks_not_introns() -> None:
    segments_by_read = {
        "readA": [
            segment(
                pos=100,
                end=300,
                readname="readA",
                aligned_blocks=[(100, 150), (250, 300)],
            )
        ],
    }
    coverage_blocks = plot_bokeh._isoseq_coverage_block_cache(
        segments_by_read,
        90,
        310,
    )
    coverage = plot_bokeh._isoseq_coverage_for_cached_reads(
        ["readA"],
        coverage_blocks,
        90,
        310,
    )

    def depth_at(position: int) -> int:
        depth = 0
        for x_value, y_value in zip(coverage["x"], coverage["y"], strict=True):
            if x_value > position:
                break
            depth = y_value
        return depth

    assert depth_at(125) == 1
    assert depth_at(100) == 0
    assert depth_at(101) == 1
    assert depth_at(150) == 1
    assert depth_at(151) == 0
    assert depth_at(200) == 0
    assert depth_at(250) == 0
    assert depth_at(251) == 1
    assert depth_at(300) == 1
    assert depth_at(301) == 0
    assert depth_at(275) == 1


def test_isoseq_exon_labels_use_local_display_rank() -> None:
    transcript = SimpleNamespace(
        transcript_id="T1",
        transcript_name="Tx1",
        gene_id="G1",
        gene_name="Gene1",
        chrom="chr1",
        start=100,
        end=200,
        strand="+",
        exons=[(100, 200, 5)],
    )
    groups = [
        {
            "group_id": "T1",
            "transcript": transcript,
            "read_names": [],
            "assigned_read_count": 0,
            "gene_transcript_count": 1,
            "is_unassigned": False,
        }
    ]
    layout = plot_bokeh._build_isoseq_layout(groups, {})
    plot_figure = figure(x_range=(50, 250), y_range=(1, 0))
    tap_tool = TapTool()
    plot_figure.add_tools(tap_tool)

    result = plot_bokeh.add_isoseq_transcripts_to_plot(
        plot_figure,
        tap_tool,
        groups,
        layout,
        50,
        250,
    )

    transcript_source = result[0]
    feature_label_source = result[4]
    assert transcript_source.data["left"] == [100]
    assert transcript_source.data["right"] == [201]
    assert feature_label_source is not None
    assert feature_label_source.data["feature_type"] == ["exon"]
    assert feature_label_source.data["text"] == ["1"]
    assert feature_label_source.data["feature_number"] == [1]
    assert feature_label_source.data["x"] == pytest.approx([150.5])


def test_isoseq_single_exon_transcript_gets_direction_marker() -> None:
    transcript = SimpleNamespace(
        transcript_id="T1",
        transcript_name="Tx1",
        gene_id="G1",
        gene_name="Gene1",
        chrom="chr1",
        start=100,
        end=200,
        strand="-",
        exons=[(100, 200, 1)],
    )
    groups = [
        {
            "group_id": "T1",
            "transcript": transcript,
            "read_names": [],
            "assigned_read_count": 0,
            "gene_transcript_count": 1,
            "is_unassigned": False,
        }
    ]
    layout = plot_bokeh._build_isoseq_layout(groups, {})
    plot_figure = figure(x_range=(150, 250), y_range=(1, 0))
    tap_tool = TapTool()
    plot_figure.add_tools(tap_tool)

    result = plot_bokeh.add_isoseq_transcripts_to_plot(
        plot_figure,
        tap_tool,
        groups,
        layout,
        150,
        250,
    )

    intron_source = result[2]
    intron_arrow_source = result[3]
    assert intron_source is None
    assert intron_arrow_source is not None
    assert intron_arrow_source.data["x"] == pytest.approx([175.5])
    assert intron_arrow_source.data["angle"] == [math.pi / 2]
    assert intron_arrow_source.data["fill_color"] == ["#FFFFFF"]
    assert intron_arrow_source.data["transcript_id"] == ["T1"]
    arrow_renderers = [
        renderer
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer)
        and renderer.data_source is intron_arrow_source
    ]
    assert arrow_renderers
    assert arrow_renderers[0].level == "overlay"


def test_isoseq_reverse_strand_exon_labels_use_transcript_order() -> None:
    transcript = SimpleNamespace(
        transcript_id="T1",
        transcript_name="Tx1",
        gene_id="G1",
        gene_name="Gene1",
        chrom="chr1",
        start=100,
        end=400,
        strand="-",
        exons=[(100, 150, 1), (200, 250, 2), (300, 400, 3)],
    )
    groups = [
        {
            "group_id": "T1",
            "transcript": transcript,
            "read_names": [],
            "assigned_read_count": 0,
            "gene_transcript_count": 1,
            "is_unassigned": False,
        }
    ]
    layout = plot_bokeh._build_isoseq_layout(groups, {})
    plot_figure = figure(x_range=(50, 450), y_range=(1, 0))
    tap_tool = TapTool()
    plot_figure.add_tools(tap_tool)

    result = plot_bokeh.add_isoseq_transcripts_to_plot(
        plot_figure,
        tap_tool,
        groups,
        layout,
        50,
        450,
    )

    transcript_source = result[0]
    intron_source = result[2]
    feature_label_source = result[4]
    assert transcript_source.data["left"] == [100, 200, 300]
    assert transcript_source.data["right"] == [151, 251, 401]
    assert intron_source is not None
    assert intron_source.data["xs"] == [[151, 200], [251, 300]]
    assert feature_label_source is not None
    exon_rows = [
        index
        for index, feature_type in enumerate(feature_label_source.data["feature_type"])
        if feature_type == "exon"
    ]
    intron_rows = [
        index
        for index, feature_type in enumerate(feature_label_source.data["feature_type"])
        if feature_type == "intron"
    ]
    assert [feature_label_source.data["text"][index] for index in exon_rows] == ["3", "2", "1"]
    assert [feature_label_source.data["feature_number"][index] for index in exon_rows] == [3, 2, 1]
    assert [feature_label_source.data["text"][index] for index in intron_rows] == ["2", "1"]
    assert [feature_label_source.data["feature_number"][index] for index in intron_rows] == [2, 1]
    assert [feature_label_source.data["feature_length"][index] for index in intron_rows] == [
        "49 bp",
        "49 bp",
    ]


def test_isoseq_controls_place_dotplot_after_gene_filter() -> None:
    transcript_source = ColumnDataSource(
        {
            "transcript_id": ["T1"],
            "transcript_name": ["Tx1"],
            "source_kind": ["transcript"],
            "gene_id": ["G1"],
            "gene_name": ["Gene1"],
        }
    )
    dotplot_thumbnail = Div(text="Show ref identity", width=156)

    controls = plot_bokeh.create_isoseq_controls(
        [transcript_source],
        [],
        dotplot_thumbnail=dotplot_thumbnail,
    )

    assert controls is not None
    assert controls.children[2] is dotplot_thumbnail
    assert controls.children[3].label == "Hide unselected isoforms"


def test_isoseq_controls_include_hide_comparison_when_comparison_sources_exist() -> None:
    transcript_source = ColumnDataSource(
        {
            "transcript_id": ["T1"],
            "transcript_name": ["Tx1"],
            "source_kind": ["transcript"],
            "gene_id": ["G1"],
            "gene_name": ["Gene1"],
            "annotation_id": ["comparison"],
            "annotation_label": ["Comparison GTF"],
        }
    )
    comparison_div = Div(text="Comparison GTF")

    controls = plot_bokeh.create_isoseq_controls(
        [transcript_source],
        [],
        comparison_components=[comparison_div],
    )

    assert controls is not None
    assert controls.children[-1].label == "Hide comparison GTF"


def test_format_alignment_coordinates_includes_exact_comma_formatted_span() -> None:
    assert (
        plot_bokeh.format_alignment_coordinates("chr1", 174_937_531, 174_938_005)
        == "chr1:174937531-174938005 (474 bp)"
    )
    assert plot_bokeh.format_alignment_coordinates("chr1", 100, 2_600) == "chr1:100-2600 (2,500 bp)"
    assert (
        plot_bokeh.format_alignment_coordinates("chr1", 100, 1_500_100)
        == "chr1:100-1500100 (1,500,000 bp)"
    )


def test_get_vcf_variant_color_by_variant_type() -> None:
    assert plot_bokeh.get_base_color("a") == "#00CC00"
    assert plot_bokeh.get_base_color("n") == "#888888"
    assert (
        plot_bokeh.get_vcf_variant_color(SimpleNamespace(variant_type="SNP", alt_base="T"))
        == "#CC0000"
    )
    assert (
        plot_bokeh.get_vcf_variant_color(SimpleNamespace(variant_type="INSERTION", alt_base=""))
        == "#333333"
    )
    assert (
        plot_bokeh.get_vcf_variant_color(SimpleNamespace(variant_type="DELETION", alt_base=""))
        == "#E8E8E8"
    )


def test_process_segments_clips_arrows_and_extracts_variants() -> None:
    segs = {
        "readA": [
            segment(
                pos=90,
                end=210,
                mismatches=[(119, "T")],
                insertions=[(129, "A"), (139, "A" * 10)],
                deletions=[(179, 181)],
                inclusion_reason="10 bp INS at chr1:130",
            )
        ],
        "readB": [segment(pos=150, end=250, is_fwd_strand=False, readname="readB")],
    }

    arrow_data, clickable_data, variant_data, connector_data, same_region_data = (
        plot_bokeh.process_segments(
            segs,
            ["readA", "readB"],
            {"readA": 1.0, "readB": 2.0},
            {"readA": 0.2, "readB": 0.2},
            100,
            200,
            "complex_sv",
            sample_label="sample",
        )
    )

    assert arrow_data["x0"] == [100, 200]
    assert arrow_data["x1"] == [200, 150]
    assert arrow_data["read_name"] == ["readA", "readB"]
    assert arrow_data["source_kind"] == ["arrow", "arrow"]
    assert arrow_data["chevron_tip_fraction"] == [
        0.65,
        plot_bokeh.PLOT_CONFIG["read_chevron_tip_fraction"],
    ]
    assert clickable_data["customdata"][0]["sample_label"] == "sample"
    assert clickable_data["customdata"][0]["inclusion_reason"] == "10 bp INS at chr1:130"
    assert clickable_data["customdata"][1]["inclusion_reason"] == ""
    assert clickable_data["customdata"][0]["coordinates"] == "chr1:90-210 (120 bp)"
    assert variant_data["mismatch"]["x"] == [120]
    assert variant_data["insertion"]["x"] == [130, 140]
    assert variant_data["insertion"]["is_1bp"] == [True, False]
    assert variant_data["deletion"]["x0"] == [180]
    assert variant_data["deletion"]["x1"] == [180]
    assert variant_data["deletion"]["is_1bp"] == [False]
    assert connector_data["stub_x0"] == []
    assert same_region_data["xs"] == []


def test_chevron_tip_fraction_avoids_deletion_under_glyph_back_edge() -> None:
    assert chevron_tip_fraction(
        100,
        200,
        [(166, 167)],
        100,
        10000,
    ) == 0.65


def test_process_segments_numbers_isoseq_blocks_by_read_orientation() -> None:
    blocks = [(100, 150), (200, 250), (300, 350)]
    reverse_segment = segment(
        pos=100,
        end=350,
        readname="readA",
        is_fwd_strand=False,
        aligned_blocks=blocks,
    )

    arrow_data, clickable_data, _variant_data, _connector_data, _same_region_data = (
        plot_bokeh.process_segments(
            {"readA": [reverse_segment]},
            ["readA"],
            {"readA": 1.0},
            {"readA": 0.2},
            90,
            360,
            plot_bokeh.ISOSEQ_REGION_TYPE,
            sample_label="sample",
        )
    )

    assert arrow_data["alignment_order"] == [3, 2, 1]
    assert arrow_data["x0"] == [151, 251, 351]
    assert arrow_data["x1"] == [101, 201, 301]
    assert arrow_data["segment_id"] == [
        "r0:row0:readA:1:chr1:100-350:block3",
        "r0:row0:readA:1:chr1:100-350:block2",
        "r0:row0:readA:1:chr1:100-350:block1",
    ]
    assert [item["alignment_number"] for item in clickable_data["customdata"]] == [3, 2, 1]
    assert clickable_data["customdata"][0]["all_alignment_numbers"] == [1, 2, 3]
    assert clickable_data["customdata"][0]["all_alignment_coordinates"] == [
        "chr1:301-351 (50 bp)",
        "chr1:201-251 (50 bp)",
        "chr1:101-151 (50 bp)",
    ]


def test_build_read_connection_index_links_by_read_order_not_reference_order() -> None:
    first = segment(
        pos=500,
        end=560,
        alignment_order=1,
        fwd_read_start=0,
        fwd_read_end=60,
    )
    second = segment(
        pos=100,
        end=160,
        alignment_order=2,
        fwd_read_start=61,
        fwd_read_end=120,
    )
    region_data = [
        {
            "region": Region("chr2", 490, 570, "chr2:490-570"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [first]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        },
        {
            "region": Region("chr1", 90, 170, "chr1:90-170"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [second]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        },
    ]

    index = plot_bokeh.build_read_connection_index(region_data)
    first_key = plot_bokeh.connection_key(0, 0, "readA", first)
    second_key = plot_bokeh.connection_key(1, 0, "readA", second)

    assert index[first_key]["next"]["target_region_number"] == 2
    assert index[first_key]["next"]["target_direction"] == "right"
    assert index[second_key]["previous"]["target_region_number"] == 1
    assert index[second_key]["previous"]["target_direction"] == "left"
    assert (
        index[first_key]["next"]["connection_id"] == index[second_key]["previous"]["connection_id"]
    )


def test_assign_region_agnostic_alignment_orders_uses_source_read_order() -> None:
    third_piece = segment(
        chrom="chr1",
        pos=120,
        end=180,
        alignment_order=1,
        fwd_read_start=200,
        fwd_read_end=260,
    )
    first_piece = segment(
        chrom="chr2",
        pos=320,
        end=380,
        alignment_order=1,
        fwd_read_start=0,
        fwd_read_end=60,
    )
    second_piece = segment(
        chrom="chr3",
        pos=520,
        end=580,
        alignment_order=1,
        fwd_read_start=100,
        fwd_read_end=160,
    )
    region_data = [
        {
            "region": Region("chr1", 100, 200, "chr1:100-200"),
            "bam_rows": [{"segments_by_read": {"readA_HP1": [third_piece]}}],
        },
        {
            "region": Region("chr2", 300, 400, "chr2:300-400"),
            "bam_rows": [{"segments_by_read": {"readA": [first_piece]}}],
        },
        {
            "region": Region("chr3", 500, 600, "chr3:500-600"),
            "bam_rows": [{"segments_by_read": {"readA": [second_piece]}}],
        },
    ]

    plot_bokeh.assign_region_agnostic_alignment_orders(region_data)

    assert first_piece.alignment_order == 1
    assert second_piece.alignment_order == 2
    assert third_piece.alignment_order == 3


def test_build_read_connection_index_includes_same_region_links() -> None:
    first = segment(pos=100, end=130, alignment_order=1, fwd_read_start=0, fwd_read_end=30)
    second = segment(pos=150, end=180, alignment_order=2, fwd_read_start=31, fwd_read_end=60)
    third = segment(pos=300, end=330, alignment_order=3, fwd_read_start=61, fwd_read_end=90)
    region_data = [
        {
            "region": Region("chr1", 90, 190, "chr1:90-190"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [first, second]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        },
        {
            "region": Region("chr2", 290, 340, "chr2:290-340"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [third]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        },
    ]

    index = plot_bokeh.build_read_connection_index(region_data)
    first_key = plot_bokeh.connection_key(0, 0, "readA", first)
    second_key = plot_bokeh.connection_key(0, 0, "readA", second)

    assert index[first_key]["next"]["connection_scope"] == "same_region"
    assert index[first_key]["next"]["target_direction"] == "same"
    assert index[first_key]["next"]["target_alignment_order"] == 2
    assert index[second_key]["previous"]["connection_scope"] == "same_region"
    assert (
        index[first_key]["next"]["connection_id"] == index[second_key]["previous"]["connection_id"]
    )
    assert index[second_key]["next"]["target_alignment_order"] == 3


def test_build_read_connection_index_uses_original_read_identity_across_display_rows() -> None:
    phased_segment = segment(
        pos=100,
        end=160,
        haplotype_tag=2,
        readname="readA",
        alignment_order=1,
        fwd_read_start=0,
        fwd_read_end=60,
    )
    unassigned_segment = segment(
        pos=300,
        end=360,
        haplotype_tag=None,
        readname="readA",
        alignment_order=2,
        fwd_read_start=61,
        fwd_read_end=121,
    )
    region_data = [
        {
            "region": Region("chr1", 90, 170, "chr1:90-170"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA_HP2": [phased_segment]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        },
        {
            "region": Region("chr2", 290, 370, "chr2:290-370"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [unassigned_segment]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        },
    ]

    index = plot_bokeh.build_read_connection_index(region_data)
    phased_key = plot_bokeh.connection_key(0, 0, "readA_HP2", phased_segment)
    unassigned_key = plot_bokeh.connection_key(1, 0, "readA", unassigned_segment)

    assert index[phased_key]["next"]["target_region_number"] == 2
    assert index[phased_key]["next"]["source_haplotype"] == "HP2"
    assert index[phased_key]["next"]["target_haplotype"] == "Unassigned"
    assert index[phased_key]["next"]["haplotype_transition"] is True
    assert (
        index[phased_key]["next"]["connection_id"]
        == index[unassigned_key]["previous"]["connection_id"]
    )


def test_build_read_connection_index_keeps_bam_rows_separate() -> None:
    top_row = segment(pos=100, end=150, alignment_order=1, readname="readA")
    bottom_row = segment(pos=300, end=350, alignment_order=2, readname="readA")
    region_data = [
        {
            "region": Region("chr1", 90, 160, "chr1:90-160"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [top_row]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        },
        {
            "region": Region("chr2", 290, 360, "chr2:290-360"),
            "bam_rows": [
                {"segments_by_read": {}, "region_type": COMPLEX_SV_REGION_TYPE},
                {
                    "segments_by_read": {"readA": [bottom_row]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                },
            ],
        },
    ]

    assert plot_bokeh.build_read_connection_index(region_data) == {}


def test_process_segments_emits_same_region_connector_line() -> None:
    first = segment(pos=100, end=140, alignment_order=1, fwd_read_start=0, fwd_read_end=40)
    second = segment(pos=180, end=220, alignment_order=2, fwd_read_start=41, fwd_read_end=80)
    region_data = [
        {
            "region": Region("chr1", 90, 230, "chr1:90-230"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [first, second]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        }
    ]
    connection_index = plot_bokeh.build_read_connection_index(region_data)

    _arrow_data, _clickable_data, _variant_data, connector_data, same_region_data = (
        plot_bokeh.process_segments(
            {"readA": [first, second]},
            ["readA"],
            {"readA": 1.0},
            {"readA": 0.3},
            90,
            230,
            COMPLEX_SV_REGION_TYPE,
            region_idx=0,
            row_index=0,
            connection_index=connection_index,
        )
    )

    assert connector_data["connection_role"] == []
    assert same_region_data["read_name"] == ["readA"]
    assert same_region_data["source_kind"] == ["connector"]
    assert same_region_data["xs"] == [[140, 140, 180, 180]]
    assert same_region_data["ys"][0][0] == pytest.approx(0.95)
    assert same_region_data["ys"][0][1] == pytest.approx(1.15)
    assert same_region_data["ys"][0][2] == pytest.approx(1.15)
    assert same_region_data["ys"][0][3] == pytest.approx(1.05)


def test_process_segments_puts_reverse_strand_previous_connector_on_read_start() -> None:
    reverse_segment = segment(
        pos=300,
        end=380,
        is_fwd_strand=False,
        alignment_order=2,
        fwd_read_start=101,
        fwd_read_end=180,
    )
    key = plot_bokeh.connection_key(1, 0, "readA", reverse_segment)
    connection_index = {
        key: {
            "previous": {
                "connection_id": "connection-1",
                "target_region_idx": 0,
                "target_region_number": 1,
                "target_alignment_order": 1,
                "target_direction": "left",
                "source_haplotype": "HP1",
                "target_haplotype": "HP1",
                "haplotype_transition": False,
            }
        }
    }

    _arrow_data, _clickable_data, _variant_data, connector_data, _same_region_data = (
        plot_bokeh.process_segments(
            {"readA": [reverse_segment]},
            ["readA"],
            {"readA": 1.0},
            {"readA": 0.2},
            300,
            400,
            COMPLEX_SV_REGION_TYPE,
            region_idx=1,
            row_index=0,
            connection_index=connection_index,
        )
    )

    assert connector_data["stub_x0"] == [380]
    assert connector_data["connection_role"] == ["previous"]
    assert connector_data["target_region"] == [1]
    assert connector_data["connection_id"] == ["connection-1"]


def test_process_segments_offsets_endpoint_rings_for_doubly_connected_segment() -> None:
    middle_segment = segment(
        pos=200,
        end=280,
        alignment_order=2,
        fwd_read_start=101,
        fwd_read_end=180,
    )
    key = plot_bokeh.connection_key(1, 0, "readA", middle_segment)
    connection_index = {
        key: {
            "previous": {
                "connection_id": "connection-previous",
                "target_region_idx": 0,
                "target_region_number": 1,
                "target_alignment_order": 1,
                "target_direction": "left",
                "source_haplotype": "HP1",
                "target_haplotype": "HP1",
                "haplotype_transition": False,
            },
            "next": {
                "connection_id": "connection-next",
                "target_region_idx": 2,
                "target_region_number": 3,
                "target_alignment_order": 3,
                "target_direction": "right",
                "source_haplotype": "HP1",
                "target_haplotype": "HP1",
                "haplotype_transition": False,
            },
        }
    }

    arrow_data, _clickable_data, _variant_data, connector_data, _same_region_data = (
        plot_bokeh.process_segments(
            {"readA": [middle_segment]},
            ["readA"],
            {"readA": 1.0},
            {"readA": 0.2},
            200,
            300,
            COMPLEX_SV_REGION_TYPE,
            region_idx=1,
            row_index=0,
            connection_index=connection_index,
        )
    )

    assert connector_data["connection_role"] == ["previous", "next"]
    assert connector_data["y"] == pytest.approx([0.97, 1.03])
    assert arrow_data["y"] == [1.0]
    assert arrow_data["y0"] == pytest.approx([0.97])
    assert arrow_data["y1"] == pytest.approx([1.03])


def test_same_region_connectors_attach_to_flat_double_connected_alignment() -> None:
    first = segment(pos=100, end=140, alignment_order=1, fwd_read_start=0, fwd_read_end=40)
    middle = segment(pos=180, end=220, alignment_order=2, fwd_read_start=41, fwd_read_end=80)
    last = segment(pos=260, end=300, alignment_order=3, fwd_read_start=81, fwd_read_end=120)
    region_data = [
        {
            "region": Region("chr1", 90, 310, "chr1:90-310"),
            "bam_rows": [
                {
                    "segments_by_read": {"readA": [first, middle, last]},
                    "region_type": COMPLEX_SV_REGION_TYPE,
                }
            ],
        }
    ]
    connection_index = plot_bokeh.build_read_connection_index(region_data)

    arrow_data, _clickable_data, _variant_data, _connector_data, same_region_data = (
        plot_bokeh.process_segments(
            {"readA": [first, middle, last]},
            ["readA"],
            {"readA": 1.0},
            {"readA": 0.4},
            90,
            310,
            COMPLEX_SV_REGION_TYPE,
            region_idx=0,
            row_index=0,
            connection_index=connection_index,
        )
    )

    assert arrow_data["y"] == pytest.approx([0.9, 1.0, 1.1])
    assert arrow_data["y0"][1] == pytest.approx(0.97)
    assert arrow_data["y1"][1] == pytest.approx(1.03)
    assert same_region_data["ys"][0][-1] == pytest.approx(1.0)
    assert same_region_data["ys"][1][0] == pytest.approx(1.0)


def test_connector_offset_min_heights_by_read_flags_doubly_connected_reads() -> None:
    middle_segment = segment(pos=200, end=280, alignment_order=2)
    key = plot_bokeh.connection_key(1, 0, "readA", middle_segment)
    connection_index = {
        key: {
            "previous": {"connection_id": "previous"},
            "next": {"connection_id": "next"},
        }
    }

    min_heights = plot_bokeh.connector_offset_min_heights_by_read(
        {"readA": [middle_segment]},
        Region("chr1", 190, 290, "chr1:190-290"),
        1,
        0,
        connection_index,
    )

    assert min_heights == {"readA": 0.18}


def test_add_read_connection_markers_to_plot_returns_selectable_source() -> None:
    fig = figure()

    result = plot_bokeh.add_read_connection_markers_to_plot(
        fig,
        {
            "stub_x0": [150],
            "stub_x1": [160],
            "y": [1],
            "marker_x": [150],
            "marker_angle": [0],
            "read_name": ["readA"],
            "source_kind": ["connector"],
            "tooltip": ["next: alignment 2 in region 2"],
            "target_region": [2],
            "target_alignment": [2],
            "connection_role": ["next"],
            "connection_id": ["connection-1"],
            "source_haplotype": ["HP1"],
            "target_haplotype": ["Unassigned"],
            "haplotype_transition": [True],
            "plot_dom_class": ["orographer-alignment-plot-r0-row0"],
            "plot_model_id": ["p1000"],
        },
    )

    assert result is not None
    source, stub_renderer, marker_renderer = result
    assert source.data["read_name"] == ["readA"]
    assert source.data["source_kind"] == ["connector"]
    assert stub_renderer is None
    assert marker_renderer.name == "orographer_read_connector_marker"


def test_add_same_region_read_connection_lines_to_plot_returns_selectable_source() -> None:
    fig = figure()

    result = plot_bokeh.add_same_region_read_connection_lines_to_plot(
        fig,
        {
            "xs": [[100, 100, 200, 200]],
            "ys": [[1.0, 1.2, 1.2, 1.1]],
            "read_name": ["readA"],
            "source_kind": ["connector"],
            "connection_id": ["connection-1"],
        },
    )

    assert result is not None
    source, renderer = result
    assert source.data["read_name"] == ["readA"]
    assert source.data["source_kind"] == ["connector"]
    assert renderer.glyph.line_dash == "dashed"
    assert renderer.glyph.line_color == plot_bokeh.PLOT_CONFIG["connector_line_color"]
    assert renderer.glyph.line_width == plot_bokeh.PLOT_CONFIG["connector_line_width"]
    hit_renderer = fig.renderers[-1]
    assert hit_renderer.name == "orographer_same_region_connector_hit_area"
    assert hit_renderer.glyph.line_alpha == 0.001
    assert hit_renderer.glyph.line_width == plot_bokeh.PLOT_CONFIG["connector_hit_line_width"]
    assert (
        renderer.selection_glyph.line_width
        == plot_bokeh.PLOT_CONFIG["connector_selection_line_width"]
    )
    assert (
        renderer.selection_glyph.line_color
        == plot_bokeh.PLOT_CONFIG["connector_selection_line_color"]
    )
    assert renderer.selection_glyph.line_alpha == "connector_selected_alpha"
    assert renderer.glyph.line_alpha == "connector_line_alpha"
    assert source.data["connector_line_alpha"] == [plot_bokeh.PLOT_CONFIG["connector_line_alpha"]]


def test_phase_block_boundaries_by_haplotype_uses_phased_visible_segments() -> None:
    segments_by_read = {
        "readA": [
            segment(pos=95, end=140, haplotype_tag=1, phaseset_tag=100),
            segment(pos=150, end=180, haplotype_tag=1, phaseset_tag="100"),
            segment(pos=220, end=260, haplotype_tag=1, phaseset_tag=220),
        ],
        "readB": [
            segment(pos=150, end=200, haplotype_tag=2, phaseset_tag=150),
            segment(pos=150, end=200, haplotype_tag=0, phaseset_tag=150),
            segment(pos=300, end=320, haplotype_tag=2, phaseset_tag=300),
        ],
        "readC": [segment(pos=120, end=160, haplotype_tag=1, phaseset_tag=None)],
    }

    boundaries = plot_bokeh.phase_block_boundaries_by_haplotype(
        segments_by_read,
        90,
        250,
    )

    assert boundaries == {1: [100, 180, 220], 2: [150, 200]}


def test_add_phase_block_boundaries_to_plot_clips_to_haplotype_lanes() -> None:
    fig = figure()

    renderer = plot_bokeh.add_phase_block_boundaries_to_plot(
        fig,
        {1: [100, 180], 2: [150], 3: [120]},
        {1: (0.5, 1.5), 2: (2.0, 3.0)},
        90,
        200,
    )

    assert renderer is not None
    assert renderer.glyph.line_color == "color"
    assert renderer.glyph.line_width == plot_bokeh.PLOT_CONFIG["phase_block_line_width"]
    assert renderer.glyph.line_alpha == plot_bokeh.PLOT_CONFIG["phase_block_line_alpha"]
    assert renderer.glyph.line_dash == plot_bokeh.PLOT_CONFIG["phase_block_line_dash"]
    assert renderer.data_source.data["xs"] == [[100, 100], [180, 180], [150, 150]]
    assert renderer.data_source.data["ys"] == [[0.5, 1.75], [0.5, 1.75], [1.75, 3.0]]
    assert renderer.data_source.data["color"] == [
        plot_bokeh.PLOT_CONFIG["sample_label_color"],
        plot_bokeh.PLOT_CONFIG["sample_label_color"],
        plot_bokeh.PLOT_CONFIG["sample_label_color"],
    ]


def test_add_same_region_connector_lines_offsets_repeated_routes() -> None:
    line_data = plot_bokeh.empty_same_region_connector_data()
    visible_endpoints = {
        "source-a": {
            "next_x": 140,
            "next_y": 1.0,
            "previous_x": 100,
            "previous_y": 1.0,
            "read_y_min": 0.9,
            "read_y_max": 1.3,
            "read_name": "readA",
        },
        "target-a": {
            "next_x": 220,
            "next_y": 1.1,
            "previous_x": 180,
            "previous_y": 1.1,
            "read_y_min": 0.9,
            "read_y_max": 1.3,
            "read_name": "readA",
        },
        "source-b": {
            "next_x": 142,
            "next_y": 1.4,
            "previous_x": 100,
            "previous_y": 1.4,
            "read_y_min": 1.3,
            "read_y_max": 1.7,
            "read_name": "readB",
        },
        "target-b": {
            "next_x": 221,
            "next_y": 1.5,
            "previous_x": 181,
            "previous_y": 1.5,
            "read_y_min": 1.3,
            "read_y_max": 1.7,
            "read_name": "readB",
        },
    }
    connection_index = {
        "source-a": {
            "next": {
                "connection_scope": "same_region",
                "target_key": "target-a",
                "connection_id": "connection-a",
            },
        },
        "source-b": {
            "next": {
                "connection_scope": "same_region",
                "target_key": "target-b",
                "connection_id": "connection-b",
            },
        },
    }

    plot_bokeh.add_same_region_connector_lines(
        line_data,
        visible_endpoints,
        connection_index,
    )

    assert line_data["read_name"] == ["readA", "readB"]
    assert line_data["ys"][0][1] == pytest.approx(0.9)
    assert line_data["ys"][1][1] == pytest.approx(1.3)


def test_same_region_connector_router_separates_overlapping_intervals() -> None:
    requests = [
        {
            "index": 0,
            "source_x": 100,
            "source_y": 1.0,
            "target_x": 220,
            "target_y": 1.05,
            "x_min": 100,
            "x_max": 220,
            "span": 120,
            "lane_min": 0.9,
            "lane_max": 1.3,
            "outer_lane_min": 0.9,
            "outer_lane_max": 1.3,
        },
        {
            "index": 1,
            "source_x": 120,
            "source_y": 1.0,
            "target_x": 240,
            "target_y": 1.05,
            "x_min": 120,
            "x_max": 240,
            "span": 120,
            "lane_min": 0.9,
            "lane_max": 1.3,
            "outer_lane_min": 0.9,
            "outer_lane_max": 1.3,
        },
    ]

    lanes = plot_bokeh.route_same_region_connectors(requests)

    assert abs(lanes[0] - lanes[1]) >= plot_bokeh.connector_lane_spacing()


def test_same_region_connector_candidates_stay_inside_shared_read_lane() -> None:
    request = {
        "source_y": 1.0,
        "target_y": 1.05,
        "lane_min": 0.9,
        "lane_max": 1.3,
        "outer_lane_min": 0.9,
        "outer_lane_max": 1.3,
    }

    candidates = plot_bokeh.connector_candidate_lanes(request)

    assert candidates
    assert all(0.9 <= candidate <= 1.3 for candidate in candidates)


def test_same_region_connector_router_stays_inside_shared_read_lane() -> None:
    requests = [
        {
            "index": 0,
            "source_x": 100,
            "source_y": 1.0,
            "target_x": 220,
            "target_y": 1.05,
            "x_min": 100,
            "x_max": 220,
            "span": 120,
            "lane_min": 0.9,
            "lane_max": 1.3,
            "outer_lane_min": 0.9,
            "outer_lane_max": 1.3,
        },
        {
            "index": 1,
            "source_x": 120,
            "source_y": 1.0,
            "target_x": 240,
            "target_y": 1.05,
            "x_min": 120,
            "x_max": 240,
            "span": 120,
            "lane_min": 0.9,
            "lane_max": 1.3,
            "outer_lane_min": 0.9,
            "outer_lane_max": 1.3,
        },
    ]

    lanes = plot_bokeh.route_same_region_connectors(requests)

    assert all(0.9 <= lane <= 1.3 for lane in lanes)


def test_same_region_connector_candidates_stay_inside_outer_read_lanes() -> None:
    request = {
        "source_y": 1.0,
        "target_y": 2.0,
        "lane_min": 1.8,
        "lane_max": 1.2,
        "outer_lane_min": 0.9,
        "outer_lane_max": 2.1,
    }

    candidates = plot_bokeh.connector_candidate_lanes(request)

    assert candidates
    assert all(0.9 <= candidate <= 2.1 for candidate in candidates)


def test_single_region_callbacks_select_connector_sources_when_present() -> None:
    fig = figure()
    arrow_source = ColumnDataSource(
        data={
            "x0": [100],
            "x1": [200],
            "y": [1],
            "read_name": ["readA"],
        }
    )
    connector_source = ColumnDataSource(
        data={
            "xs": [[100, 100, 200, 200]],
            "ys": [[1, 1.1, 1.1, 1]],
            "read_name": ["readA"],
            "source_kind": ["connector"],
            "connection_id": ["connection-1"],
        }
    )

    plot_bokeh.add_multi_region_callbacks(
        [fig],
        [arrow_source],
        [None],
        [arrow_source, connector_source],
    )

    callback = fig.js_event_callbacks["tap"][0]
    assert callback.args["source"] is arrow_source
    assert callback.args["all_sources"] == [arrow_source, connector_source]


def test_add_arrows_to_plot_returns_source_with_interior_chevrons() -> None:
    fig = figure()
    source, renderer = plot_bokeh.add_arrows_to_plot(
        fig,
        {
            "x0": [100, 200],
            "x1": [200, 100],
            "y": [1, 2],
            "y0": [0.97, 1.97],
            "y1": [1.03, 2.03],
            "color": ["red", "blue"],
            "read_name": ["fwd", "rev"],
        },
    )

    assert isinstance(source, ColumnDataSource)
    assert renderer is not None
    assert renderer.glyph.y0 == "y"
    assert renderer.glyph.y1 == "y"
    assert renderer.nonselection_glyph.line_alpha == "arrow_nonselected_alpha"
    assert source.data["arrow_nonselected_alpha"] == [
        plot_bokeh.PLOT_CONFIG["arrow_nonselection_line_alpha"],
        plot_bokeh.PLOT_CONFIG["arrow_nonselection_line_alpha"],
    ]
    assert source.data["y0"] == [0.97, 1.97]
    assert source.data["y1"] == [1.03, 2.03]
    assert source.data["angle"][0] < 0
    assert source.data["angle"][1] > 0
    chevron_xs = source.data["chevron_xs"]
    chevron_ys = source.data["chevron_ys"]

    assert chevron_xs[0] == pytest.approx([172, 180, 172])
    assert chevron_ys[0] == pytest.approx([1.035, 1, 0.965])
    assert min(value for value in source.data["chevron_xs"][0] if not math.isnan(value)) >= 100
    assert max(value for value in source.data["chevron_xs"][0] if not math.isnan(value)) <= 200
    assert min(value for value in source.data["chevron_xs"][1] if not math.isnan(value)) >= 100
    assert max(value for value in source.data["chevron_xs"][1] if not math.isnan(value)) <= 200
    assert fig.renderers[-1].glyph.xs == "chevron_xs"
    assert fig.renderers[-1].glyph.ys == "chevron_ys"


def test_add_arrows_to_plot_uses_per_row_chevron_tip_fraction() -> None:
    fig = figure()
    source, _renderer = plot_bokeh.add_arrows_to_plot(
        fig,
        {
            "x0": [100, 200],
            "x1": [200, 100],
            "y": [1, 2],
            "color": ["red", "blue"],
            "read_name": ["fwd", "rev"],
            "chevron_tip_fraction": [0.65, 0.65],
        },
    )

    assert source.data["chevron_xs"][0] == pytest.approx([157, 165, 157])
    assert source.data["chevron_ys"][0] == pytest.approx([1.035, 1, 0.965])
    assert source.data["chevron_xs"][1] == pytest.approx([143, 135, 143])
    assert source.data["chevron_ys"][1] == pytest.approx([2.035, 2, 1.965])


def test_add_variants_to_plot_splits_one_bp_renderers() -> None:
    fig = figure()
    renderers = plot_bokeh.add_variants_to_plot(
        fig,
        {
            "mismatch": {"x": [110], "y": [1], "alt": ["A"], "color": ["#00CC00"]},
            "insertion": {
                "x": [120, 130],
                "y": [1, 1],
                "size": [2, 6],
                "count": [1, 10],
                "is_1bp": [True, False],
            },
            "deletion": {"x0": [140], "x1": [140], "y": [1], "is_1bp": [True]},
        },
    )

    assert len(renderers["marker"]) == 3
    assert len(renderers["text"]) == 3
    assert len(renderers["one_bp_markers"]) == 1
    assert len(renderers["one_bp_texts"]) == 1
    assert len(renderers["one_bp_segments"]) == 1


def test_setup_variant_lod_callback_attaches_range_callbacks() -> None:
    fig = figure()
    marker = fig.scatter([1], [1])
    text = fig.text([1], [1], text=["A"])
    checkbox = SimpleNamespace(active=False)

    plot_bokeh.setup_variant_lod_callback(
        fig,
        {"marker": [marker], "text": [text]},
        one_bp_renderers=[marker],
        hide_1bp_checkbox=checkbox,
    )

    callbacks = fig.x_range.js_property_callbacks
    assert "change:start" in callbacks
    assert "change:end" in callbacks
    callback = callbacks["change:start"][0]
    assert isinstance(callback, CustomJS)
    assert callback.args["one_bp_renderers"] == [marker]
    assert callback.args["hide_1bp_checkbox"] is checkbox


def test_add_vcf_variants_filters_by_visible_region() -> None:
    fig = figure()
    variants = [
        SimpleNamespace(
            chrom="chr1",
            pos=120,
            ref="A",
            alt="T",
            alt_base="T",
            variant_type="SNP",
            haplotypes=["1"],
        ),
        SimpleNamespace(
            chrom="chr1",
            pos=500,
            ref="A",
            alt="AC",
            alt_base="",
            variant_type="INSERTION",
            haplotypes=[],
        ),
    ]

    plot_bokeh.add_vcf_variants(fig, variants, 100, 200, sample_label="sample")

    source = fig.renderers[0].data_source
    assert source.data["x"] == [120]
    assert source.data["variant_type"] == ["SNP"]
    assert source.data["sample_label"] == ["sample"]
    assert fig.toolbar.active_tap is not None


def test_add_clickable_labels_wires_tap_and_reset_callbacks() -> None:
    fig = figure()
    tap_tool = TapTool()
    fig.add_tools(tap_tool)
    arrow_source = ColumnDataSource({"x0": [1], "x1": [2], "y": [1]})
    arrow_renderer = fig.segment([1], [1], [2], [1])

    result = plot_bokeh.add_clickable_labels(
        fig,
        tap_tool,
        {
            "x": [1],
            "y": [2],
            "customdata": [
                {
                    "read_name": "readA",
                    "alignment_number": 3,
                    "strand": "Forward (+)",
                    "coordinates": "chr1:1-2",
                    "haplotype": "HP:1",
                    "sample_label": "sample",
                }
            ],
        },
        arrow_source=arrow_source,
        arrow_renderer=arrow_renderer,
    )

    assert result is not None
    source, renderers = result
    assert renderers is not None
    assert source.data["label_alpha"] == [0.8]
    assert tap_tool.renderers
    assert "tap" in fig.js_event_callbacks
    assert "reset" in fig.js_event_callbacks
    assert not fig.select(type=HoverTool)


def test_alignment_label_selection_callbacks_watch_selection_sources() -> None:
    selection_source = ColumnDataSource({"read_name": ["readA"]})
    label_source = ColumnDataSource({"read_name": ["readA"], "label_alpha": [0.8]})

    plot_bokeh.add_alignment_label_selection_callbacks(
        [selection_source],
        [label_source],
    )

    callbacks = selection_source.selected.js_property_callbacks["change:indices"]
    assert len(callbacks) == 1
    assert callbacks[0].args["selection_sources"] == [selection_source]
    assert callbacks[0].args["alignment_label_sources"] == [label_source]


def test_get_expected_haplotypes_only_for_complex_sv_coverage() -> None:
    assert plot_bokeh._get_expected_haplotypes(
        {"coverage_tracks": {-1: ([], []), 0: ([], []), 2: ([], [])}, "region_type": "complex_sv"}
    ) == {0, 2}
    assert (
        plot_bokeh._get_expected_haplotypes(
            {"coverage_tracks": {-1: ([], []), 1: ([], [])}, "region_type": "paraphase"}
        )
        is None
    )


def test_paraphase_global_controls_hide_alignment_numbers_by_default() -> None:
    region_state = plot_bokeh.RegionBuildState()
    region_state.region_type = PARAPHASE_REGION_TYPE

    control_renderers = plot_bokeh._global_control_renderers([(0, {}, region_state)])

    assert control_renderers["default_hide_alignment_numbers"] is True


def test_non_paraphase_global_controls_show_alignment_numbers_by_default() -> None:
    region_state = plot_bokeh.RegionBuildState()
    region_state.region_type = COMPLEX_SV_REGION_TYPE

    control_renderers = plot_bokeh._global_control_renderers([(0, {}, region_state)])

    assert control_renderers["default_hide_alignment_numbers"] is False


def test_complex_sv_rows_use_region_local_y_layouts() -> None:
    region_one_state = plot_bokeh.RegionBuildState()
    region_two_state = plot_bokeh.RegionBuildState()
    region_one = Region("chr1", 100, 300, "chr1:100-300")
    region_two = Region("chr2", 100, 300, "chr2:100-300")
    read_a_first_region = segment(readname="readA", haplotype_tag=2)
    read_b_first_region = segment(readname="readB", haplotype_tag=1, pos=150, end=220)
    read_a_second_region = segment(readname="readA", haplotype_tag=2)

    plot_bokeh._build_bam_row_track(
        {
            "segments_by_read": {
                "readA": [read_a_first_region],
                "readB": [read_b_first_region],
            },
            "region_type": COMPLEX_SV_REGION_TYPE,
            "vcf_variants": [],
        },
        region_one,
        region_one.start,
        region_one.end,
        region_one_state,
        region_idx=0,
        row_index=0,
    )
    plot_bokeh._build_bam_row_track(
        {
            "segments_by_read": {"readA": [read_a_second_region]},
            "region_type": COMPLEX_SV_REGION_TYPE,
            "vcf_variants": [],
        },
        region_two,
        region_two.start,
        region_two.end,
        region_two_state,
        region_idx=1,
        row_index=0,
    )

    region_one_arrow_data = region_one_state.arrow_sources[0].data
    region_two_arrow_data = region_two_state.arrow_sources[0].data
    region_one_read_a_index = region_one_arrow_data["read_name"].index("readA")
    region_two_read_a_index = region_two_arrow_data["read_name"].index("readA")

    assert (
        region_one_arrow_data["y"][region_one_read_a_index]
        != region_two_arrow_data["y"][region_two_read_a_index]
    )
    assert region_two_arrow_data["y"][region_two_read_a_index] == 0.1


def test_add_multi_region_callbacks_single_and_multi_region() -> None:
    single_fig = figure()
    source_a = ColumnDataSource({"x0": [1], "x1": [2], "y": [1]})
    plot_bokeh.add_multi_region_callbacks([single_fig], [source_a], [object()])
    assert "tap" in single_fig.js_event_callbacks
    assert "reset" in single_fig.js_event_callbacks

    first_fig = figure()
    second_fig = figure()
    source_b = ColumnDataSource({"x0": [2], "x1": [3], "y": [1]})
    bounds = {
        "x_range": Range1d(0, 10),
        "x_start": 0,
        "x_end": 10,
        "y_ranges": [Range1d(2, 0)],
        "y_bounds": [[2, 0]],
        "all_sources": [source_a, source_b],
    }
    plot_bokeh.add_multi_region_callbacks(
        [first_fig, second_fig],
        [source_a, source_b],
        [object(), object()],
        all_region_reset_bounds=[([first_fig, second_fig], [source_a], [object()], bounds)],
    )

    assert "tap" in first_fig.js_event_callbacks
    assert "tap" in second_fig.js_event_callbacks
    assert "reset" in first_fig.js_event_callbacks
    assert "reset" in second_fig.js_event_callbacks


def test_create_read_search_button_wires_callback() -> None:
    source = ColumnDataSource({"x": [1]})
    controls = plot_bokeh.create_read_search_button([source])
    read_button, clear_button = controls.children

    assert controls.sizing_mode == "fixed"
    assert controls.width == 212
    assert controls.height == 22
    assert read_button.height == 22
    assert clear_button.height == 22
    assert clear_button.disabled
    assert read_button.label == "Select reads"
    assert clear_button.label == "Clear selected"
    assert "button_click" in read_button.js_event_callbacks
    assert "button_click" in clear_button.js_event_callbacks
    assert clear_button.styles["width"] == "90px"
    assert read_button.styles["width"] == "90px"
    assert "change:indices" in source.selected.js_property_callbacks


def test_plot_reads_bokeh_returns_none_for_empty_input() -> None:
    assert plot_bokeh.plot_reads_bokeh([], OutputConfig("out", None)) is None


def test_plot_reads_bokeh_saves_no_data_layout(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    captured: dict = {}
    region = Region("chr1", 100, 200, "chr1:100-200")

    def save_plot(layout, output_file: str, prefix: str | None) -> None:
        captured["layout"] = layout
        captured["output_file"] = output_file
        captured["prefix"] = prefix

    monkeypatch.setattr(plot_bokeh, "save_plot_with_modal", save_plot)

    output = plot_bokeh.plot_reads_bokeh(
        [
            {
                "region": region,
                "gene_annotations": [],
                "bam_rows": [
                    {
                        "segments_by_read": {},
                        "region_type": "paraphase",
                    }
                ],
            }
        ],
        OutputConfig(str(tmp_path), "sample"),
    )

    assert output == str(tmp_path / "sample_chr1_100_200_bokeh.html")
    assert captured["output_file"] == output
    assert captured["prefix"] == "sample"


def test_plot_reads_bokeh_preserves_isoseq_track_height_for_dual_gtf(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: dict = {}
    region = Region("chr1", 100, 250, "chr1:100-250")

    def save_plot(layout, output_file: str, prefix: str | None) -> None:
        captured["layout"] = layout
        captured["output_file"] = output_file
        captured["prefix"] = prefix

    original_prepare_chunks = plot_bokeh.isoseq_tracks.prepare_lazy_chunks
    read_page_sizes = []

    def prepare_chunks_with_capture(*args, **kwargs):
        read_page_sizes.append(kwargs.get("read_page_size"))
        return original_prepare_chunks(*args, **kwargs)

    def group(transcript_id: str, annotation_id: str, annotation_label: str) -> dict:
        transcript = SimpleNamespace(
            transcript_id=transcript_id,
            transcript_name=transcript_id,
            gene_id="G1",
            gene_name="Gene1",
            chrom="chr1",
            start=100,
            end=200,
            strand="+",
            exons=[(100, 200, 1)],
        )
        return {
            "group_id": transcript_id,
            "transcript": transcript,
            "read_names": [],
            "assigned_read_count": 0,
            "gene_transcript_count": 1,
            "is_unassigned": False,
            "annotation_id": annotation_id,
            "annotation_label": annotation_label,
        }

    monkeypatch.setattr(plot_bokeh, "save_plot_with_modal", save_plot)
    monkeypatch.setattr(
        plot_bokeh.isoseq_tracks,
        "prepare_lazy_chunks",
        prepare_chunks_with_capture,
    )

    output = plot_bokeh.plot_reads_bokeh(
        [
            {
                "region": region,
                "gene_annotations": [],
                "isoseq": True,
                "bam_rows": [
                    {
                        "segments_by_read": {},
                        "region_type": ISOSEQ_REGION_TYPE,
                        "reference_path": None,
                        "vcf_variants": [],
                        "sample_label": "sample",
                        "isoseq_groups": [group("TX1", "primary", "Primary GTF")],
                        "isoseq_annotation_tracks": [
                            {
                                "annotation_id": "primary",
                                "annotation_label": "Primary GTF",
                                "isoseq_groups": [group("TX1", "primary", "Primary GTF")],
                            },
                            {
                                "annotation_id": "comparison",
                                "annotation_label": "Comparison GTF",
                                "isoseq_groups": [
                                    group("TX2", "comparison", "Comparison GTF")
                                ],
                            },
                        ],
                    }
                ],
            }
        ],
        OutputConfig(str(tmp_path), "sample"),
    )

    final_layout = captured["layout"]
    region_row = final_layout.children[1]

    assert output == str(tmp_path / "sample_chr1_100_250_bokeh.html")
    assert final_layout.sizing_mode == "stretch_both"
    assert region_row.sizing_mode == "stretch_both"
    assert region_row.children[0].sizing_mode == "stretch_both"
    isoseq_figures = []
    components = [final_layout]
    while components:
        component = components.pop()
        css_classes = getattr(component, "css_classes", [])
        if any(css_class.startswith("orographer-isoseq-plot") for css_class in css_classes):
            isoseq_figures.append(component)
        components.extend(getattr(component, "children", []))

    assert len(isoseq_figures) == 2
    assert [figure.height for figure in isoseq_figures] == [480, 480]
    assert [figure.sizing_mode for figure in isoseq_figures] == ["stretch_both", "stretch_both"]
    assert read_page_sizes == [50, 50]
