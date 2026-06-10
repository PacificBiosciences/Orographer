"""Unit tests for BFS state helpers in complex_sv (no real BAM needed)."""

import pytest

from orographer.bam_parser import AlignmentSummary
from orographer.complex_sv import (
    INCLUSION_REASON_CHIMERIC,
    INCLUSION_REASON_LARGE_DELETION,
    INCLUSION_REASON_LARGE_INSERTION,
    INCLUSION_REASON_LARGE_SOFTCLIP,
    _add_visited_span,
    _build_final_regions,
    _cigar_ref_len,
    _truncate_regions,
    discover_connected_reads,
    inclusion_reason_for_record,
)
from orographer.utils import Region


def summary(
    read_name: str,
    *,
    chromosome: str = "chr1",
    start: int = 1000,
    end: int = 1100,
    has_sa: bool = False,
    sa_tag: str | None = None,
    max_indel: int = 0,
    max_softclip: int = 0,
) -> AlignmentSummary:
    return AlignmentSummary(
        read_name=read_name,
        chromosome=chromosome,
        start=start,
        end=end,
        is_supplementary=False,
        has_sa=has_sa,
        sa_tag=sa_tag,
        max_indel=max_indel,
        max_softclip=max_softclip,
    )


def sa_segment(
    *,
    rname: str = "chr2",
    pos: int = 2000,
    cigar: str = "100M",
    mapq: int = 60,
):
    return type(
        "SplitSegmentStub",
        (),
        {
            "rname": rname,
            "pos": pos,
            "cigar": cigar,
            "mapq": mapq,
        },
    )()


class RecordStub:
    def __init__(
        self,
        *,
        is_supplementary: bool = False,
        has_sa: bool = False,
        cigartuples: list[tuple[int, int]] | None = None,
        reference_start: int = 100,
        reference_name: str = "chr1",
    ) -> None:
        self.is_supplementary = is_supplementary
        self._has_sa = has_sa
        self.cigartuples = cigartuples or []
        self.reference_start = reference_start
        self.reference_name = reference_name

    def has_tag(self, tag: str) -> bool:
        return tag == "SA" and self._has_sa


def test_inclusion_reason_for_record_matches_complex_sv_discovery_triggers() -> None:
    assert inclusion_reason_for_record(RecordStub(has_sa=True)) == INCLUSION_REASON_CHIMERIC
    assert (
        inclusion_reason_for_record(RecordStub(is_supplementary=True)) == INCLUSION_REASON_CHIMERIC
    )
    assert (
        inclusion_reason_for_record(RecordStub(cigartuples=[(0, 30), (1, 21)]))
        == f"21 bp {INCLUSION_REASON_LARGE_INSERTION} at chr1:130"
    )
    assert (
        inclusion_reason_for_record(
            RecordStub(
                cigartuples=[(1, 1568)],
                reference_start=187_478_993,
            )
        )
        == f"1,568 bp {INCLUSION_REASON_LARGE_INSERTION} at chr1:187478994"
    )
    assert (
        inclusion_reason_for_record(RecordStub(cigartuples=[(0, 30), (2, 21)]))
        == f"21 bp {INCLUSION_REASON_LARGE_DELETION} at chr1:131"
    )
    assert (
        inclusion_reason_for_record(
            RecordStub(
                cigartuples=[(2, 1568)],
                reference_start=187_478_993,
            )
        )
        == f"1,568 bp {INCLUSION_REASON_LARGE_DELETION} at chr1:187478994"
    )
    assert (
        inclusion_reason_for_record(RecordStub(cigartuples=[(4, 21)]))
        == INCLUSION_REASON_LARGE_SOFTCLIP
    )
    assert inclusion_reason_for_record(RecordStub(cigartuples=[(0, 100)])) == ""


class TestAddVisitedSpan:
    def test_first_span_creates_list(self):
        visited: dict = {}
        _add_visited_span(visited, "chr1", 100, 200)
        assert visited == {"chr1": [(100, 200)]}

    def test_non_overlapping_append(self):
        visited = {"chr1": [(100, 200)]}
        _add_visited_span(visited, "chr1", 300, 400)
        assert visited["chr1"] == [(100, 200), (300, 400)]

    def test_overlapping_merges(self):
        visited = {"chr1": [(100, 200)]}
        _add_visited_span(visited, "chr1", 150, 300)
        assert visited["chr1"] == [(100, 300)]

    def test_adjacent_merges(self):
        visited = {"chr1": [(100, 200)]}
        _add_visited_span(visited, "chr1", 200, 300)
        assert visited["chr1"] == [(100, 300)]

    def test_contained_span(self):
        visited = {"chr1": [(50, 400)]}
        _add_visited_span(visited, "chr1", 100, 200)
        assert visited["chr1"] == [(50, 400)]

    def test_multiple_merges(self):
        visited = {"chr1": [(100, 200), (250, 350), (400, 500)]}
        _add_visited_span(visited, "chr1", 150, 450)
        assert visited["chr1"] == [(100, 500)]

    def test_different_chroms_independent(self):
        visited = {"chr1": [(100, 200)]}
        _add_visited_span(visited, "chr2", 100, 200)
        assert "chr1" in visited
        assert "chr2" in visited
        assert visited["chr1"] == [(100, 200)]
        assert visited["chr2"] == [(100, 200)]


class TestCigarRefLen:
    def test_match(self):
        assert _cigar_ref_len("100M") == 100

    def test_deletion_counts(self):
        assert _cigar_ref_len("50M50D50M") == 150

    def test_insertion_does_not_count(self):
        assert _cigar_ref_len("50M10I50M") == 100

    def test_soft_clip_does_not_count(self):
        assert _cigar_ref_len("10S80M10S") == 80

    def test_complex_cigar(self):
        assert _cigar_ref_len("30M5I20M3D10M") == 63

    def test_equal_and_diff(self):
        assert _cigar_ref_len("50=50X") == 100


class TestBuildFinalRegions:
    def test_empty_inputs(self):
        result = _build_final_regions({}, [])
        assert result == []

    def test_single_visited_interval(self):
        # visited_intervals are 1-based inclusive
        result = _build_final_regions({}, [("chr1", 1000, 2000)])
        assert len(result) >= 1
        r = result[0]
        # Expanded by 40% again, so result is wider than input
        assert r.chromosome == "chr1"
        assert r.start <= 1000
        assert r.end >= 2000

    def test_empty_coordinate_str(self):
        result = _build_final_regions({}, [("chr1", 1000, 2000)])
        for r in result:
            assert r.coordinate_str == ""

    def test_summary_extents_included(self):
        summaries = {
            "r1": [
                AlignmentSummary(
                    read_name="r1",
                    chromosome="chr1",
                    start=5000,  # 0-based
                    end=5100,  # 0-based exclusive
                    is_supplementary=False,
                    has_sa=False,
                    sa_tag=None,
                )
            ]
        }
        result = _build_final_regions(summaries, [])
        # Summary at 0-based [5000, 5100) = 1-based [5001, 5100]
        # After expand_interval the region should contain this range
        assert len(result) >= 1
        r = result[0]
        assert r.start <= 5001
        assert r.end >= 5100

    def test_multiple_chroms_sorted(self):
        result = _build_final_regions(
            {},
            [("chr12", 1000, 2000), ("chr8", 1000, 2000), ("chr1", 1000, 2000)],
        )
        chroms = [r.chromosome for r in result]
        assert chroms == ["chr1", "chr8", "chr12"]

    def test_nearby_intervals_merge(self):
        # Two intervals close enough that after 10% expansion they merge.
        # [1000, 2000]: width=1001, margin=101 → expanded [899, 2101]
        # [2001, 3000]: width=1000, margin=100 → expanded [1901, 3100]
        # 1901 < 2101 → they overlap → merge into one region
        result = _build_final_regions(
            {},
            [("chr1", 1000, 2000), ("chr1", 2001, 3000)],
        )
        assert len(result) == 1

    def test_distant_intervals_stay_separate(self):
        # Intervals farther apart than the final-region merge distance stay separate.
        result = _build_final_regions(
            {},
            [("chr1", 1000, 2000), ("chr1", 300_000, 400_000)],
        )
        assert len(result) == 2

    def test_same_chromosome_final_regions_merge_within_200kb(self):
        result = _build_final_regions(
            {},
            [("chr1", 1000, 2000), ("chr1", 100_000, 120_000)],
        )

        assert len(result) == 1
        assert result[0].chromosome == "chr1"

    def test_different_chromosome_final_regions_do_not_merge_within_200kb(self):
        result = _build_final_regions(
            {},
            [("chr1", 1000, 2000), ("chr2", 100_000, 120_000)],
        )

        assert len(result) == 2


class TestTruncateRegions:
    def test_keeps_regions_until_cumulative_limit(self):
        regions = [
            Region("chr1", 1, 100, "chr1:1-100"),
            Region("chr2", 1, 100, "chr2:1-100"),
            Region("chr3", 1, 100, "chr3:1-100"),
        ]

        assert _truncate_regions(regions, max_bp=200) == regions[:2]

    def test_continues_after_dropping_oversized_region(self):
        regions = [
            Region("chr1", 1, 100, "chr1:1-100"),
            Region("chr2", 1, 500, "chr2:1-500"),
            Region("chr3", 1, 50, "chr3:1-50"),
        ]

        assert _truncate_regions(regions, max_bp=160) == [regions[0], regions[2]]


class TestDiscoverConnectedReads:
    def test_discovers_reads_from_sa_indel_and_softclip_summaries(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ):
        fetch_regions: list[Region] = []

        def fetch_alignment_summaries(_bam_file, region: Region, _mapq_threshold: int):
            fetch_regions.append(region)
            if region.chromosome == "chr1":
                return [
                    summary("sa-read", has_sa=True, sa_tag="sa-tag"),
                    summary("indel-read", max_indel=21),
                    summary("softclip-read", max_softclip=21),
                    summary("boring-read"),
                ]
            return []

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            fetch_alignment_summaries,
        )
        monkeypatch.setattr(
            "orographer.complex_sv.parse_sa_aux_val",
            lambda _tag: [sa_segment(rname="chr2", pos=1999, cigar="50M", mapq=60)],
        )

        discovered, summaries_by_read, visited = discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
            region_connection_threshold=1,
        )

        assert discovered == {"sa-read", "indel-read", "softclip-read"}
        assert set(summaries_by_read) == {"sa-read", "indel-read", "softclip-read"}
        assert [region.chromosome for region in fetch_regions] == ["chr1", "chr2"]
        assert visited[0][0] == "chr1"
        assert visited[1][0] == "chr2"

    def test_ignores_low_mapq_sa_segments(self, monkeypatch: pytest.MonkeyPatch):
        fetch_regions: list[Region] = []

        def fetch_alignment_summaries(_bam_file, region: Region, _mapq_threshold: int):
            fetch_regions.append(region)
            return [summary("sa-read", has_sa=True, sa_tag="sa-tag")]

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            fetch_alignment_summaries,
        )
        monkeypatch.setattr(
            "orographer.complex_sv.parse_sa_aux_val",
            lambda _tag: [sa_segment(mapq=10)],
        )

        discovered, summaries_by_read, visited = discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
        )

        assert discovered == {"sa-read"}
        assert set(summaries_by_read) == {"sa-read"}
        assert len(fetch_regions) == 1
        assert len(visited) == 1

    def test_default_threshold_does_not_follow_single_read_sa_target(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ):
        fetch_regions: list[Region] = []

        def fetch_alignment_summaries(_bam_file, region: Region, _mapq_threshold: int):
            fetch_regions.append(region)
            if region.chromosome == "chr1":
                return [summary("read-a", has_sa=True, sa_tag="tag-a")]
            return []

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            fetch_alignment_summaries,
        )
        monkeypatch.setattr(
            "orographer.complex_sv.parse_sa_aux_val",
            lambda _tag: [sa_segment(rname="chr2", pos=1999, cigar="50M", mapq=60)],
        )

        discovered, summaries_by_read, visited = discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
        )

        assert discovered == {"read-a"}
        assert set(summaries_by_read) == {"read-a"}
        assert [region.chromosome for region in fetch_regions] == ["chr1"]
        assert len(visited) == 1

    def test_threshold_follows_merged_target_locus_with_two_unique_reads(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ):
        fetch_regions: list[Region] = []

        def fetch_alignment_summaries(_bam_file, region: Region, _mapq_threshold: int):
            fetch_regions.append(region)
            if region.chromosome == "chr1":
                return [
                    summary("read-a", has_sa=True, sa_tag="tag-a"),
                    summary("read-b", has_sa=True, sa_tag="tag-b"),
                ]
            return []

        def parse_sa(tag: str):
            if tag == "tag-a":
                return [sa_segment(rname="chr2", pos=1999, cigar="100M", mapq=60)]
            return [sa_segment(rname="chr2", pos=2049, cigar="100M", mapq=60)]

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            fetch_alignment_summaries,
        )
        monkeypatch.setattr("orographer.complex_sv.parse_sa_aux_val", parse_sa)

        _, _, visited = discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
        )

        assert [region.chromosome for region in fetch_regions] == ["chr1", "chr2"]
        assert visited[1][0] == "chr2"

    def test_multiple_sa_entries_from_one_read_count_once(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ):
        fetch_regions: list[Region] = []

        def fetch_alignment_summaries(_bam_file, region: Region, _mapq_threshold: int):
            fetch_regions.append(region)
            if region.chromosome == "chr1":
                return [summary("read-a", has_sa=True, sa_tag="tag-a")]
            return []

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            fetch_alignment_summaries,
        )
        monkeypatch.setattr(
            "orographer.complex_sv.parse_sa_aux_val",
            lambda _tag: [
                sa_segment(rname="chr2", pos=1999, cigar="100M", mapq=60),
                sa_segment(rname="chr2", pos=2049, cigar="100M", mapq=60),
            ],
        )

        discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
        )

        assert [region.chromosome for region in fetch_regions] == ["chr1"]

    def test_low_mapq_sa_segment_does_not_contribute_threshold_support(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ):
        fetch_regions: list[Region] = []

        def fetch_alignment_summaries(_bam_file, region: Region, _mapq_threshold: int):
            fetch_regions.append(region)
            if region.chromosome == "chr1":
                return [
                    summary("read-a", has_sa=True, sa_tag="tag-a"),
                    summary("read-b", has_sa=True, sa_tag="tag-b"),
                ]
            return []

        def parse_sa(tag: str):
            if tag == "tag-a":
                return [sa_segment(rname="chr2", pos=1999, cigar="100M", mapq=60)]
            return [sa_segment(rname="chr2", pos=2049, cigar="100M", mapq=10)]

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            fetch_alignment_summaries,
        )
        monkeypatch.setattr("orographer.complex_sv.parse_sa_aux_val", parse_sa)

        discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
        )

        assert [region.chromosome for region in fetch_regions] == ["chr1"]

    def test_parses_each_sa_tag_once(self, monkeypatch: pytest.MonkeyPatch):
        parse_calls: list[str] = []

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            lambda _bam, _region, _mapq: [
                summary("read-a", has_sa=True, sa_tag="same-tag"),
                summary("read-b", has_sa=True, sa_tag="same-tag"),
            ],
        )

        def parse_sa(tag: str):
            parse_calls.append(tag)
            return []

        monkeypatch.setattr("orographer.complex_sv.parse_sa_aux_val", parse_sa)

        discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
        )

        assert parse_calls == ["same-tag"]

    def test_truncates_when_interval_limit_is_reached(self, monkeypatch: pytest.MonkeyPatch):
        monkeypatch.setattr("orographer.complex_sv.MAX_DISCOVERY_INTERVALS", 2)
        fetched: list[Region] = []

        def fetch_alignment_summaries(_bam_file, region: Region, _mapq_threshold: int):
            fetched.append(region)
            return []

        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            fetch_alignment_summaries,
        )

        _, _, visited = discover_connected_reads(
            bam_file=object(),
            initial_frontier=[
                ("chr1", 1000, 1100),
                ("chr2", 2000, 2100),
                ("chr3", 3000, 3100),
            ],
            mapq_threshold=20,
        )

        assert len(fetched) == 2
        assert len(visited) == 2

    def test_stops_when_connected_read_limit_is_reached(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ):
        monkeypatch.setattr("orographer.complex_sv.MAX_CONNECTED_READS", 1)
        monkeypatch.setattr(
            "orographer.complex_sv.fetch_alignment_summaries",
            lambda _bam, _region, _mapq: [
                summary("read-a", max_indel=21),
                summary("read-b", max_softclip=21),
            ],
        )

        discovered, _, _ = discover_connected_reads(
            bam_file=object(),
            initial_frontier=[("chr1", 1000, 1100)],
            mapq_threshold=20,
        )

        assert discovered == {"read-a", "read-b"}
