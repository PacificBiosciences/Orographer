"""Unit tests for dotplot rendering helpers in plot_bokeh.figures."""

import pytest

from orographer.plot_bokeh.figures import (
    _build_block_ticks,
    _format_genomic_pos,
    _nice_tick_interval,
)


class TestNiceTickInterval:
    def test_returns_positive_integer(self):
        assert _nice_tick_interval(1_000) > 0
        assert isinstance(_nice_tick_interval(1_000), int)

    def test_small_span_still_returns_positive(self):
        assert _nice_tick_interval(5) >= 1

    def test_136kb_gives_reasonable_interval(self):
        interval = _nice_tick_interval(136_153)
        # Should produce roughly 4 ticks → interval in range 25000-50000
        assert 10_000 <= interval <= 100_000

    def test_1mb_gives_reasonable_interval(self):
        interval = _nice_tick_interval(1_200_000)
        assert 100_000 <= interval <= 500_000

    def test_tick_count_is_near_target(self):
        for span in [5_000, 33_670, 136_153, 1_200_000]:
            interval = _nice_tick_interval(span, target_count=4)
            approx_ticks = span / interval
            assert 2 <= approx_ticks <= 8, f"span={span} interval={interval} ticks={approx_ticks:.1f}"


class TestFormatGenomicPos:
    def test_megabase_range(self):
        label = _format_genomic_pos(55_600_000)
        assert "Mb" in label
        assert "55" in label

    def test_kilobase_range(self):
        label = _format_genomic_pos(79_747)
        assert "kb" in label

    def test_small_bp(self):
        label = _format_genomic_pos(500)
        assert "bp" in label

    def test_one_mb_boundary(self):
        label = _format_genomic_pos(1_000_000)
        assert "Mb" in label


class TestBuildBlockTicks:
    def _make_block(self, chrom, start, end, offset_start):
        span = end - start
        return {
            "chromosome": chrom,
            "start": start,
            "end": end,
            "offset_start": offset_start,
            "offset_end": offset_start + span,
            "label": f"{chrom}:{start}-{end}",
        }

    def test_returns_two_parallel_lists(self):
        block = self._make_block("chr11", 55_569_788, 55_705_940, 0)
        positions, labels = _build_block_ticks([block])
        assert len(positions) == len(labels)

    def test_positions_are_within_block_offsets(self):
        block = self._make_block("chr11", 55_569_788, 55_705_940, 0)
        positions, _ = _build_block_ticks([block])
        for pos in positions:
            assert block["offset_start"] <= pos < block["offset_end"]

    def test_two_blocks_produce_non_overlapping_positions(self):
        b1 = self._make_block("chr11", 55_569_788, 55_705_940, 0)
        b2 = self._make_block("chr16", 78_500_000, 79_700_000, b1["offset_end"])
        positions, labels = _build_block_ticks([b1, b2])
        assert len(positions) >= 2
        # All b2 positions must be > all b1 positions
        b1_end = b1["offset_end"]
        b1_positions = [p for p in positions if p < b1_end]
        b2_positions = [p for p in positions if p >= b1_end]
        assert len(b1_positions) >= 1
        assert len(b2_positions) >= 1

    def test_labels_contain_chromosome_name(self):
        block = self._make_block("chr7", 40_820_003, 40_856_731, 0)
        _, labels = _build_block_ticks([block])
        for label in labels:
            assert "chr7" in label

    def test_positions_are_floats(self):
        block = self._make_block("chr11", 55_569_788, 55_705_940, 0)
        positions, _ = _build_block_ticks([block])
        for pos in positions:
            assert isinstance(pos, float)
