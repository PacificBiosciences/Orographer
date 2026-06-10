"""Unit tests for complex_sv interval free functions."""

import math

from orographer.complex_sv import (
    expand_interval,
    merge_sorted_intervals,
    subtract_visited_spans,
)


class TestExpandInterval:
    def test_basic_expansion(self):
        chrom, start, end = expand_interval("chr1", 1000, 2000)
        margin = math.ceil(1001 * 0.10)  # 101
        assert chrom == "chr1"
        assert start == 1000 - margin
        assert end == 2000 + margin

    def test_clamps_start_to_one(self):
        _, start, _ = expand_interval("chr1", 50, 100)
        assert start >= 1

    def test_custom_fraction(self):
        _, start, end = expand_interval("chr1", 1000, 2000, fraction=0.10)
        margin = math.ceil(1001 * 0.10)
        assert start == 1000 - margin
        assert end == 2000 + margin

    def test_end_always_expands(self):
        _, _, end = expand_interval("chr1", 1000, 2000)
        assert end > 2000

    def test_start_always_contracts(self):
        _, start, _ = expand_interval("chr1", 1000, 2000)
        assert start < 1000


class TestMergeSortedIntervals:
    def test_empty(self):
        assert merge_sorted_intervals([]) == []

    def test_single(self):
        result = merge_sorted_intervals([("chr1", 100, 200)])
        assert result == [("chr1", 100, 200)]

    def test_non_overlapping(self):
        intervals = [("chr1", 100, 200), ("chr1", 300, 400)]
        result = merge_sorted_intervals(intervals)
        assert result == [("chr1", 100, 200), ("chr1", 300, 400)]

    def test_overlapping(self):
        intervals = [("chr1", 100, 250), ("chr1", 200, 400)]
        result = merge_sorted_intervals(intervals)
        assert result == [("chr1", 100, 400)]

    def test_adjacent(self):
        intervals = [("chr1", 100, 200), ("chr1", 200, 300)]
        result = merge_sorted_intervals(intervals)
        assert result == [("chr1", 100, 300)]

    def test_merges_same_chromosome_intervals_within_configured_gap(self):
        intervals = [("chr1", 100, 200), ("chr1", 350, 400)]
        result = merge_sorted_intervals(intervals, max_gap_bp=200)
        assert result == [("chr1", 100, 400)]

    def test_does_not_merge_different_chromosomes_within_configured_gap(self):
        intervals = [("chr1", 100, 200), ("chr2", 250, 300)]
        result = merge_sorted_intervals(intervals, max_gap_bp=200)
        assert result == [("chr1", 100, 200), ("chr2", 250, 300)]

    def test_multiple_chroms(self):
        intervals = [
            ("chr2", 100, 200),
            ("chr1", 300, 400),
            ("chr1", 350, 500),
        ]
        result = merge_sorted_intervals(intervals)
        assert ("chr1", 300, 500) in result
        assert ("chr2", 100, 200) in result

    def test_unsorted_input(self):
        intervals = [("chr1", 300, 400), ("chr1", 100, 200)]
        result = merge_sorted_intervals(intervals)
        assert result == [("chr1", 100, 200), ("chr1", 300, 400)]

    def test_contained_interval(self):
        intervals = [("chr1", 100, 500), ("chr1", 200, 300)]
        result = merge_sorted_intervals(intervals)
        assert result == [("chr1", 100, 500)]

    def test_chrom_order(self):
        intervals = [("chr2", 100, 200), ("chr1", 100, 200)]
        result = merge_sorted_intervals(intervals)
        assert result[0][0] == "chr1"
        assert result[1][0] == "chr2"


class TestSubtractVisitedSpans:
    def test_empty_candidates(self):
        assert subtract_visited_spans([], {}) == []

    def test_no_visited(self):
        candidates = [("chr1", 100, 300)]
        result = subtract_visited_spans(candidates, {})
        assert result == [("chr1", 100, 300)]

    def test_fully_visited(self):
        candidates = [("chr1", 100, 300)]
        visited = {"chr1": [(50, 350)]}
        result = subtract_visited_spans(candidates, visited)
        assert result == []

    def test_partial_left(self):
        candidates = [("chr1", 100, 300)]
        visited = {"chr1": [(50, 200)]}
        result = subtract_visited_spans(candidates, visited)
        assert len(result) == 1
        chrom, start, end = result[0]
        assert chrom == "chr1"
        assert start == 200
        assert end == 300

    def test_partial_right(self):
        candidates = [("chr1", 100, 300)]
        visited = {"chr1": [(200, 350)]}
        result = subtract_visited_spans(candidates, visited)
        assert result == [("chr1", 100, 200)]

    def test_visited_in_middle(self):
        candidates = [("chr1", 100, 400)]
        visited = {"chr1": [(200, 300)]}
        result = subtract_visited_spans(candidates, visited)
        assert ("chr1", 100, 200) in result
        assert ("chr1", 300, 400) in result

    def test_different_chrom_not_subtracted(self):
        candidates = [("chr1", 100, 300)]
        visited = {"chr2": [(50, 350)]}
        result = subtract_visited_spans(candidates, visited)
        assert result == [("chr1", 100, 300)]

    def test_multiple_candidates(self):
        candidates = [("chr1", 100, 300), ("chr1", 500, 700)]
        visited = {"chr1": [(50, 200)]}
        result = subtract_visited_spans(candidates, visited)
        assert ("chr1", 200, 300) in result
        assert ("chr1", 500, 700) in result
