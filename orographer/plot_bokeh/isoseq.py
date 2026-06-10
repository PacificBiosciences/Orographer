"""IsoSeq-specific plot data helpers shared by Bokeh rendering paths."""

import hashlib
import math
import re
from typing import Any

from .segment_processing import display_block_coordinates
from .sources import add_read_chevron_data, empty_arrow_source_data
from .utils import PLOT_CONFIG

ISOSEQ_READ_PAGE_SIZE = 100
ISOSEQ_PAGE_TEMPLATE_TOKEN = "__PAGE__"
ISOSEQ_READ_SHARD_SIZE = 2000


def stable_gene_color(gene_id: str) -> str:
    """Return a stable color for transcript glyphs from the gene of origin."""
    palette = [
        "#3366CC",
        "#109618",
        "#990099",
        "#FF9900",
        "#0099C6",
        "#DD4477",
        "#66AA00",
        "#B82E2E",
        "#316395",
        "#994499",
    ]
    digest = hashlib.sha1(str(gene_id or "Unknown").encode("utf-8")).hexdigest()
    return palette[int(digest[:8], 16) % len(palette)]


def format_span(start: int, end: int) -> str:
    """Return a comma-formatted inclusive span label."""
    return f"{max(0, end - start + 1):,} bp"


def intron_direction_positions(
    intron_start: int,
    intron_end: int,
    region_span: int,
) -> list[float]:
    """Return capped, evenly spaced x positions for transcript strand markers."""
    intron_length = max(0, intron_end - intron_start)
    if intron_length == 0:
        return []
    marker_count = max(1, min(6, int(intron_length / max(60, region_span / 30)) + 1))
    return [
        intron_start + (intron_length * (marker_idx + 1) / (marker_count + 1))
        for marker_idx in range(marker_count)
    ]


def safe_chunk_token(value: Any) -> str:
    """Return a filesystem-safe token for lazy-loaded chunk files."""
    token = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "unknown")).strip("._")
    return token[:80] or "unknown"


def natural_sort_key(value: Any) -> list:
    """Return a case-insensitive sort key that treats digit runs numerically."""
    key = []
    for part in re.split(r"(\d+)", str(value)):
        if not part:
            continue
        key.append(int(part) if part.isdigit() else part.casefold())
    return key


def block_number_by_index(segment_blocks: list[tuple[int, int]], is_fwd_strand: bool) -> list[int]:
    """Return block labels in read-orientation order."""
    block_count = len(segment_blocks)
    return [
        block_index + 1 if is_fwd_strand else block_count - block_index
        for block_index in range(block_count)
    ]


def read_shard_record(read_name: str, segments: list[Any], read_metadata: dict) -> list:
    """Return the compact read record stored in IsoSeq MessagePack shards."""
    segment = segments[0]
    segment_blocks = getattr(segment, "aligned_blocks", None) or [(segment.pos, segment.end)]
    block_numbers = block_number_by_index(segment_blocks, segment.is_fwd_strand)
    blocks = [
        [
            int(display_start),
            int(display_end),
            int(block_numbers[block_index]),
        ]
        for block_index, (start, end) in enumerate(segment_blocks)
        for display_start, display_end in [display_block_coordinates(start, end)]
    ]
    haplotype = getattr(segment, "haplotype_tag", None) or 0
    return [
        read_name,
        read_metadata.get("gene_id", ""),
        read_metadata.get("gene_name", ""),
        read_metadata.get("transcript_id", ""),
        read_metadata.get("group_id", ""),
        int(haplotype),
        bool(segment.is_fwd_strand),
        int(getattr(segment, "fwd_read_start", 0)),
        int(getattr(segment, "fwd_read_end", 0)),
        blocks,
    ]


def build_layout(groups: list[dict], _segments_by_read: dict) -> dict:
    """Return row positions for IsoSeq transcript/read groups."""
    transcript_height = 0.36
    read_height = 0.24
    group_gap = 0.18
    selected_read_limit = ISOSEQ_READ_PAGE_SIZE
    y = 0.0
    transcript_y = {}
    read_to_y = {}
    read_to_y_bottom = {}
    read_heights = {}
    read_names = []
    read_metadata = {}

    for group in groups:
        if y > 0:
            y += group_gap
        transcript_y[group["group_id"]] = y + transcript_height / 2
        y += transcript_height
        transcript = group.get("transcript")
        for read_name in group.get("read_names", []):
            read_metadata[read_name] = {
                "gene_id": transcript.gene_id if transcript else "",
                "gene_name": transcript.gene_name if transcript else "Unassigned",
                "transcript_id": transcript.transcript_id if transcript else "UNASSIGNED",
                "group_id": group["group_id"],
            }

    selected_read_y_start = y + 1.0
    selected_view_height = selected_read_y_start + selected_read_limit * read_height

    return {
        "transcript_y": transcript_y,
        "read_names": read_names,
        "read_to_y": read_to_y,
        "read_to_y_bottom": read_to_y_bottom,
        "read_heights": read_heights,
        "height": max(selected_view_height, 1.0),
        "transcript_area_height": max(y, 1.0),
        "selected_read_y_start": selected_read_y_start,
        "read_metadata": read_metadata,
    }


def add_metadata_to_sources(
    arrow_data: dict,
    clickable_data: dict,
    read_metadata: dict,
) -> None:
    """Attach transcript/gene metadata to read glyph and label sources."""
    for key in ("gene_id", "gene_name", "transcript_id", "group_id"):
        arrow_data[key] = []
    for read_name in arrow_data.get("read_name", []):
        meta = read_metadata.get(read_name, {})
        arrow_data["gene_id"].append(meta.get("gene_id", ""))
        arrow_data["gene_name"].append(meta.get("gene_name", ""))
        arrow_data["transcript_id"].append(meta.get("transcript_id", ""))
        arrow_data["group_id"].append(meta.get("group_id", ""))

    for info in clickable_data.get("customdata", []):
        meta = read_metadata.get(info.get("read_name", ""), {})
        info["gene_id"] = meta.get("gene_id", "")
        info["gene_name"] = meta.get("gene_name", "")
        info["transcript_id"] = meta.get("transcript_id", "")


def normalise_arrow_source_data(arrow_data: dict) -> dict:
    """Return read glyph data with every browser source column populated."""
    arrow_source_data = dict(empty_arrow_source_data())
    arrow_source_data.update(arrow_data)
    row_count = len(arrow_source_data["x0"])
    if len(arrow_source_data.get("y0", [])) != row_count:
        arrow_source_data["y0"] = list(arrow_source_data["y"])
    if len(arrow_source_data.get("y1", [])) != row_count:
        arrow_source_data["y1"] = list(arrow_source_data["y"])
    defaults = {
        "read_filter_alpha": 1.0,
        "has_split_alignment": False,
        "has_multiregion_connection": False,
        "arrow_line_alpha": PLOT_CONFIG["arrow_line_alpha"],
        "arrowhead_alpha": PLOT_CONFIG["arrowhead_fill_alpha"],
        "arrow_selected_alpha": 1.0,
        "arrow_nonselected_alpha": PLOT_CONFIG["arrow_nonselection_line_alpha"],
    }
    for key, default_value in defaults.items():
        if len(arrow_source_data.get(key, [])) != row_count:
            arrow_source_data[key] = [default_value] * row_count
    arrow_source_data["angle"] = [
        -math.pi / 2 if x1 > x0 else math.pi / 2
        for x0, x1 in zip(arrow_source_data["x0"], arrow_source_data["x1"], strict=True)
    ]
    add_read_chevron_data(arrow_source_data)
    return arrow_source_data
