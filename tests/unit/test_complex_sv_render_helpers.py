from __future__ import annotations

from types import SimpleNamespace

from orographer.complex_sv import (
    add_record_coverage,
    build_coverage_tracks,
    collect_insertion_summary,
    create_coverage_state,
    group_segments_by_haplotype,
)
from orographer.utils import Region


class CoverageRecord:
    def __init__(
        self,
        blocks: list[tuple[int, int]],
        hp_tag: int | str | None = None,
    ) -> None:
        self._blocks = blocks
        self._hp_tag = hp_tag

    def get_blocks(self) -> list[tuple[int, int]]:
        return self._blocks

    def get_tag(self, tag: str) -> int | str:
        if tag == "HP" and self._hp_tag is not None:
            return self._hp_tag
        raise KeyError(tag)


def segment(
    *,
    haplotype_tag: int | None,
    fwd_read_start: int = 0,
    insertions: list[tuple[int, str]] | None = None,
) -> SimpleNamespace:
    return SimpleNamespace(
        haplotype_tag=haplotype_tag,
        fwd_read_start=fwd_read_start,
        alignment_order=0,
        insertions=insertions or [],
    )


def test_add_record_coverage_populates_total_and_haplotype_tracks() -> None:
    region = Region("chr1", 101, 130, "chr1:101-130")
    coverage_state = create_coverage_state(region, haplotypes_to_track={0, 1, 2}, bin_size=10)
    record = CoverageRecord(blocks=[(100, 120)], hp_tag=1)

    add_record_coverage(record, coverage_state, min_cigar_gap=20)
    tracks = build_coverage_tracks(coverage_state)

    assert tracks[-1] == ([101, 111, 121], [1, 1, 0])
    assert tracks[1] == ([101, 111, 121], [1, 1, 0])
    assert tracks[0] == ([101, 111, 121], [0, 0, 0])
    assert tracks[2] == ([101, 111, 121], [0, 0, 0])


def test_add_record_coverage_tracks_unassigned_hp_when_requested() -> None:
    region = Region("chr1", 101, 120, "chr1:101-120")
    coverage_state = create_coverage_state(region, haplotypes_to_track={0}, bin_size=10)
    record = CoverageRecord(blocks=[(100, 110)])

    add_record_coverage(record, coverage_state, min_cigar_gap=20)
    tracks = build_coverage_tracks(coverage_state)

    assert tracks[-1] == ([101, 111], [1, 0])
    assert tracks[0] == ([101, 111], [1, 0])


def test_add_record_coverage_tracks_observed_string_haplotypes() -> None:
    region = Region("chr1", 101, 120, "chr1:101-120")
    coverage_state = create_coverage_state(
        region,
        haplotypes_to_track={0, 1, 2},
        bin_size=10,
        include_observed_haplotypes=True,
    )
    record = CoverageRecord(blocks=[(100, 110)], hp_tag="custom_hap_a")

    add_record_coverage(record, coverage_state, min_cigar_gap=20)
    tracks = build_coverage_tracks(coverage_state)

    assert tracks[-1] == ([101, 111], [1, 0])
    assert tracks["custom_hap_a"] == ([101, 111], [1, 0])
    assert tracks[1] == ([101, 111], [0, 0])


def test_group_segments_by_haplotype_orders_segments_and_uses_hp_keys() -> None:
    late = segment(haplotype_tag=2, fwd_read_start=100)
    early = segment(haplotype_tag=2, fwd_read_start=10)
    unassigned = segment(haplotype_tag=None, fwd_read_start=50)

    grouped = group_segments_by_haplotype(
        {
            "readA": [late, early],
            "readB": [unassigned],
        }
    )

    assert grouped["readA_HP2"] == [early, late]
    assert grouped["readB"] == [unassigned]
    assert [early.alignment_order, late.alignment_order] == [1, 2]
    assert unassigned.alignment_order == 1


def test_collect_insertion_summary_filters_deduplicates_and_truncates_names() -> None:
    region = Region("chr1", 100, 200, "chr1:100-200")
    segments_by_read = {
        "r1": [
            segment(
                haplotype_tag=1,
                insertions=[(109, "A" * 30), (109, "C" * 45), (150, "G" * 20)],
            )
        ],
        "r2": [segment(haplotype_tag=1, insertions=[(109, "A" * 40)])],
        "r3": [segment(haplotype_tag=1, insertions=[(109, "A" * 50)])],
        "r4": [segment(haplotype_tag=1, insertions=[(109, "A" * 60)])],
        "r5": [segment(haplotype_tag=1, insertions=[(109, "A" * 70)])],
        "r6": [segment(haplotype_tag=1, insertions=[(109, "A" * 80)])],
        "outside": [segment(haplotype_tag=1, insertions=[(250, "A" * 90)])],
        "unassigned": [segment(haplotype_tag=None, insertions=[(119, "T" * 25)])],
    }

    summary = collect_insertion_summary(segments_by_read, region)

    assert summary[1] == [
        {
            "pos": 110,
            "count": 6,
            "median_size": 55,
            "read_names": ["r1", "r2", "r3", "r4", "r5"],
            "total_count": 6,
            "chrom": "chr1",
        }
    ]
    assert summary[0] == [
        {
            "pos": 120,
            "count": 1,
            "median_size": 25,
            "read_names": ["unassigned"],
            "total_count": 1,
            "chrom": "chr1",
        }
    ]
