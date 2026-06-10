"""Coverage data builders shared by region-specific plot tracks."""

from collections.abc import Iterable, Mapping
from typing import Any

from .sources import empty_isoseq_coverage_data


def coverage_blocks_for_read(
    read_name: str,
    segments_by_read: Mapping[str, list[Any]],
    coordinate_start: int,
    coordinate_end: int,
) -> Iterable[tuple[int, int]]:
    """Yield clipped half-open coverage blocks for one read."""
    for segment in segments_by_read.get(read_name, []):
        blocks = getattr(segment, "aligned_blocks", None) or [(segment.pos, segment.end)]
        for block_start, block_end in blocks:
            start = max(coordinate_start, block_start + 1)
            end = min(coordinate_end + 1, block_end + 1)
            if start >= end:
                continue
            yield start, end


def coverage_block_cache(
    segments_by_read: Mapping[str, list[Any]],
    coordinate_start: int,
    coordinate_end: int,
) -> dict[str, list[tuple[int, int]]]:
    """Cache per-read coverage blocks once for coverage assembly."""
    return {
        read_name: list(
            coverage_blocks_for_read(
                read_name,
                segments_by_read,
                coordinate_start,
                coordinate_end,
            )
        )
        for read_name in segments_by_read
    }


def coverage_for_cached_reads(
    read_names: Iterable[str],
    coverage_blocks: Mapping[str, list[tuple[int, int]]],
    coordinate_start: int,
    coordinate_end: int,
) -> dict[str, list]:
    """Return step coverage data for cached read names."""
    blocks = (block for read_name in read_names for block in coverage_blocks.get(read_name, []))
    return coverage_for_blocks(blocks, coordinate_start, coordinate_end)


def coverage_for_blocks(
    blocks: Iterable[tuple[int, int]],
    coordinate_start: int,
    coordinate_end: int,
) -> dict[str, list]:
    """Return step coverage source data from half-open coverage blocks."""
    events: dict[int, int] = {}
    for start, end in blocks:
        events[start] = events.get(start, 0) + 1
        events[end] = events.get(end, 0) - 1
    if not events:
        return empty_isoseq_coverage_data()

    x_vals = [coordinate_start]
    y_vals = [0]
    depth = 0
    for pos in sorted(events):
        if pos < coordinate_start or pos > coordinate_end:
            continue
        x_vals.append(pos)
        y_vals.append(depth)
        depth += events[pos]
        x_vals.append(pos)
        y_vals.append(depth)
    x_vals.append(coordinate_end)
    y_vals.append(depth)
    return {
        "x": x_vals,
        "y": y_vals,
        "coverage_alpha": [0.72] * len(x_vals),
    }
