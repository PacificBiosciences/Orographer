from __future__ import annotations

from dataclasses import dataclass

from orographer.plot_bokeh.data import (
    READ_FILTER_LAYOUT_MODES,
    ReadFilterLayoutRequest,
    _has_nearby_segments,
    _has_overlapping_segments,
    build_read_filter_layouts,
    calculate_read_positions,
    get_read_haplotype,
    sort_read_names,
)


@dataclass(frozen=True)
class SegmentStub:
    chrom: str
    pos: int
    end: int
    haplotype_tag: int | str | None = None


def test_get_read_haplotype_normalizes_unassigned_values() -> None:
    assert get_read_haplotype([]) == 0
    assert get_read_haplotype([SegmentStub("chr1", 10, 20, None)]) == 0
    assert get_read_haplotype([SegmentStub("chr1", 10, 20, 0)]) == 0
    assert get_read_haplotype([SegmentStub("chr1", 10, 20, "Unassigned")]) == 0
    assert get_read_haplotype([SegmentStub("chr1", 10, 20, "2")]) == 2


def test_sort_read_names_groups_phased_reads_before_unassigned() -> None:
    segments_by_read = {
        "hp2_low": [SegmentStub("chr1", 10, 50, 2)],
        "hp1": [SegmentStub("chr1", 20, 60, 1)],
        "unassigned": [SegmentStub("chr1", 30, 70, None)],
        "hp2_high": [SegmentStub("chr1", 40, 80, 2)],
    }

    read_names, haplotype_groups, haplotype_order = sort_read_names(
        segments_by_read,
        expected_haplotypes={1, 2, 3},
    )

    assert haplotype_order == [1, 2, 3, 0]
    assert haplotype_groups[2] == ["hp2_high", "hp2_low"]
    assert haplotype_groups[3] == []
    assert read_names == ["hp1", "hp2_high", "hp2_low", "unassigned"]


def test_has_overlapping_segments_only_compares_within_chromosome() -> None:
    same_chrom_overlap = [
        SegmentStub("chr1", 10, 30),
        SegmentStub("chr1", 25, 40),
    ]
    different_chrom_overlap = [
        SegmentStub("chr1", 10, 30),
        SegmentStub("chr2", 25, 40),
    ]

    assert _has_overlapping_segments(same_chrom_overlap)
    assert not _has_overlapping_segments(different_chrom_overlap)


def test_has_nearby_segments_accepts_close_non_overlapping_same_chromosome_segments() -> None:
    same_chrom_nearby = [
        SegmentStub("chr1", 10, 30),
        SegmentStub("chr1", 40, 60),
    ]
    same_chrom_far = [
        SegmentStub("chr1", 10, 30),
        SegmentStub("chr1", 41, 60),
    ]

    assert _has_nearby_segments(same_chrom_nearby)
    assert not _has_nearby_segments(same_chrom_far)
    assert not _has_overlapping_segments(same_chrom_nearby)


def test_calculate_read_positions_expands_overlapping_and_nearby_split_reads() -> None:
    segments_by_read = {
        "split_overlap": [
            SegmentStub("chr1", 10, 40, 1),
            SegmentStub("chr1", 30, 60, 1),
        ],
        "split_nearby": [
            SegmentStub("chr1", 100, 130, 1),
            SegmentStub("chr1", 140, 170, 1),
        ],
        "split_far": [
            SegmentStub("chr1", 200, 230, 1),
            SegmentStub("chr1", 241, 270, 1),
        ],
    }
    haplotype_groups = {1: ["split_overlap", "split_nearby", "split_far"], 2: []}

    read_to_y, read_to_y_bottom, heights, total_height, boundaries = calculate_read_positions(
        read_names=["split_overlap", "split_nearby", "split_far"],
        segments_by_read=segments_by_read,
        haplotype_groups=haplotype_groups,
        haplotype_order=[1, 2],
        min_read_height=0.2,
        multi_alignment_row_spacing=0.25,
    )

    assert heights["split_overlap"] == 0.75
    assert heights["split_nearby"] == 0.75
    assert heights["split_far"] == 0.2
    assert read_to_y["split_overlap"] == 0.375
    assert read_to_y_bottom["split_nearby"] == 0.75
    assert read_to_y_bottom["split_far"] == 1.5
    assert boundaries[1] == (0.0, 1.7)
    assert boundaries[2] == (2.05, 2.11)
    assert total_height == 2.11


def test_calculate_read_positions_compacts_adjacent_empty_haplotype_groups() -> None:
    _read_to_y, _read_to_y_bottom, _heights, total_height, boundaries = calculate_read_positions(
        read_names=[],
        segments_by_read={},
        haplotype_groups={1: [], 2: [], 0: []},
        haplotype_order=[1, 2, 0],
    )

    assert boundaries[1] == (0.0, 0.06)
    assert boundaries[2] == (0.41, 0.47)
    assert boundaries[0] == (0.82, 0.8799999999999999)
    assert total_height == 0.8799999999999999


def test_calculate_read_positions_respects_per_read_minimum_heights() -> None:
    segments_by_read = {
        "connected": [SegmentStub("chr1", 10, 40, 1)],
        "plain": [SegmentStub("chr1", 50, 80, 1)],
    }
    haplotype_groups = {1: ["connected", "plain"]}

    read_to_y, read_to_y_bottom, heights, total_height, boundaries = calculate_read_positions(
        read_names=["connected", "plain"],
        segments_by_read=segments_by_read,
        haplotype_groups=haplotype_groups,
        haplotype_order=[1],
        min_read_height=0.2,
        min_read_heights_by_read={"connected": 0.44},
    )

    assert heights["connected"] == 0.44
    assert heights["plain"] == 0.2
    assert read_to_y["connected"] == 0.22
    assert read_to_y_bottom["plain"] == 0.44
    assert boundaries[1] == (0.0, 0.64)
    assert total_height == 0.64


def test_build_read_filter_layouts_compacts_by_evidence_mode() -> None:
    segments_by_read = {
        "split_only": [SegmentStub("chr1", 10, 40, 1)],
        "plain": [SegmentStub("chr1", 50, 80, 1)],
        "multiregion": [SegmentStub("chr1", 90, 120, 2)],
    }
    request = ReadFilterLayoutRequest(
        read_names=["split_only", "plain", "multiregion"],
        segments_by_read=segments_by_read,
        haplotype_groups={1: ["split_only", "plain"], 2: ["multiregion"]},
        haplotype_order=[1, 2],
        read_filter_flags_by_read={
            "split_only": (True, False),
            "plain": (False, False),
            "multiregion": (True, True),
        },
    )

    layouts = build_read_filter_layouts(request)

    assert tuple(layouts) == READ_FILTER_LAYOUT_MODES
    assert layouts["all"]["read_names"] == ["split_only", "plain", "multiregion"]
    assert layouts["split"]["read_names"] == ["split_only", "multiregion"]
    assert layouts["multiregion"]["read_names"] == ["multiregion"]
    assert layouts["split_multiregion"]["read_names"] == ["multiregion"]
    assert layouts["split"]["alignments_height"] < layouts["all"]["alignments_height"]
    assert layouts["multiregion"]["read_to_y"] == {"multiregion": 0.51}
