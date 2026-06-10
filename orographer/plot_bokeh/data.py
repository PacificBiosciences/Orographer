"""Data and layout helpers for Bokeh plotting (read ordering, positions, filenames)."""

import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

READ_FILTER_LAYOUT_MODES = ("all", "split", "multiregion", "split_multiregion")


def add_read_filter_visibility_columns(source_data: dict) -> None:
    """Add one precomputed visibility column per read-filter layout mode."""
    row_count = 0
    for values in source_data.values():
        if isinstance(values, list):
            row_count = len(values)
            break
    if row_count == 0:
        return

    split_values = source_data.get("has_split_alignment", [True] * row_count)
    multiregion_values = source_data.get("has_multiregion_connection", [True] * row_count)
    source_data["read_filter_visible_all"] = [1.0] * row_count
    source_data["read_filter_visible_split"] = [
        1.0 if bool(value) else 0.0 for value in split_values
    ]
    source_data["read_filter_visible_multiregion"] = [
        1.0 if bool(value) else 0.0 for value in multiregion_values
    ]
    source_data["read_filter_visible_split_multiregion"] = [
        1.0 if bool(split_value) and bool(multiregion_value) else 0.0
        for split_value, multiregion_value in zip(
            split_values,
            multiregion_values,
            strict=True,
        )
    ]


@dataclass(frozen=True)
class ReadFilterLayoutRequest:
    """Inputs needed to precompute read layouts for evidence-filter modes."""

    read_names: Sequence[str]
    segments_by_read: Mapping[str, Sequence[object]]
    haplotype_groups: Mapping[object, Sequence[str]]
    haplotype_order: Sequence[object]
    read_filter_flags_by_read: Mapping[str, tuple[bool, bool]]
    min_read_heights_by_read: Mapping[str, float] | None = None


def _read_visible_for_filter_mode(has_split: bool, has_multiregion: bool, mode: str) -> bool:
    """Return whether a read should be included in a compacted layout mode."""
    if mode == "split":
        return has_split
    if mode == "multiregion":
        return has_multiregion
    if mode == "split_multiregion":
        return has_split and has_multiregion
    return True


def _filtered_haplotype_groups(
    haplotype_groups: Mapping[object, Sequence[str]],
    selected_reads: set[str],
) -> dict[object, list[str]]:
    """Return haplotype groups with non-selected reads removed."""
    return {
        haplotype: [read_name for read_name in group_reads if read_name in selected_reads]
        for haplotype, group_reads in haplotype_groups.items()
    }


def build_read_filter_layouts(
    request: ReadFilterLayoutRequest,
) -> dict[str, dict[str, object]]:
    """Precompute y layouts for all read evidence filter modes."""
    min_read_heights = request.min_read_heights_by_read or {}
    layouts: dict[str, dict[str, object]] = {}
    for mode in READ_FILTER_LAYOUT_MODES:
        selected_reads = {
            read_name
            for read_name in request.read_names
            if _read_visible_for_filter_mode(
                *request.read_filter_flags_by_read.get(read_name, (False, False)),
                mode,
            )
        }
        mode_read_names = [
            read_name for read_name in request.read_names if read_name in selected_reads
        ]
        mode_segments_by_read = {
            read_name: request.segments_by_read[read_name] for read_name in mode_read_names
        }
        mode_haplotype_groups = _filtered_haplotype_groups(
            request.haplotype_groups,
            selected_reads,
        )
        mode_min_heights = {
            read_name: height
            for read_name, height in min_read_heights.items()
            if read_name in selected_reads
        }
        read_to_y, read_to_y_bottom, read_heights, alignments_height, group_boundaries = (
            calculate_read_positions(
                mode_read_names,
                mode_segments_by_read,
                mode_haplotype_groups,
                request.haplotype_order,
                min_read_heights_by_read=mode_min_heights,
            )
        )
        layouts[mode] = {
            "read_names": mode_read_names,
            "read_to_y": read_to_y,
            "read_to_y_bottom": read_to_y_bottom,
            "read_heights": read_heights,
            "alignments_height": alignments_height,
            "group_boundaries": group_boundaries,
            "haplotype_order": list(request.haplotype_order),
        }
    return layouts


def generate_output_filename(
    segments_by_read, coordinate_start, coordinate_end, output_dir, prefix
):
    """Generate output filename from coordinates."""
    first_segment = next(iter(segments_by_read.values()))[0]
    chrom = first_segment.chrom

    # Build filename with optional prefix
    if prefix:
        output_file = f"{prefix}_{chrom}_{coordinate_start}_{coordinate_end}_bokeh.html"
    else:
        output_file = f"{chrom}_{coordinate_start}_{coordinate_end}_bokeh.html"

    # Ensure output_dir exists and join with filename
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, output_file)

    return output_file


def _has_nearby_segments(segments, max_gap_bp=10) -> bool:
    """Return True if any same-chromosome segments overlap or are very close."""
    by_chrom: dict = {}
    for seg in segments:
        by_chrom.setdefault(seg.chrom, []).append(seg)
    for chrom_segs in by_chrom.values():
        if len(chrom_segs) < 2:
            continue
        chrom_segs = sorted(chrom_segs, key=lambda s: s.pos)
        for i in range(len(chrom_segs) - 1):
            if chrom_segs[i].end + max_gap_bp >= chrom_segs[i + 1].pos:
                return True
    return False


def _has_overlapping_segments(segments) -> bool:
    """Return True if any two segments of a read overlap on the same chromosome."""
    return _has_nearby_segments(segments, max_gap_bp=0)


def get_read_haplotype(segments):
    """Get the haplotype tag from the first segment of a read.
    Treats None, 0, and `Unassigned` as unassigned (group 0)."""
    if not segments:
        return 0
    hp = segments[0].haplotype_tag
    if hp is None or hp == 0 or str(hp).lower() == "unassigned":
        return 0
    try:
        return int(hp)
    except (TypeError, ValueError):
        pass
    return hp


def sort_read_names(segments_by_read, expected_haplotypes=None):
    """
    Sort read names by haplotype group, then by start position.
    Dynamically detects all haplotypes present in the data.

    expected_haplotypes: optional set of haplotype ints that should appear in
    haplotype_order even if no reads belong to them (e.g. from coverage_tracks).

    Returns:
        tuple: (sorted_read_names, haplotype_groups, haplotype_order) where:
               - haplotype_groups: haplotype -> list of read names (empty list for absent HPs)
               - haplotype_order: display order (phased first, then unassigned)
    """
    # First pass: discover all unique haplotypes
    haplotype_groups = {}

    for read_name in segments_by_read:
        haplotype = get_read_haplotype(segments_by_read[read_name])
        if haplotype not in haplotype_groups:
            haplotype_groups[haplotype] = []
        haplotype_groups[haplotype].append(read_name)

    # Include expected haplotypes that have no reads (empty placeholders)
    if expected_haplotypes:
        for hp in expected_haplotypes:
            if hp not in haplotype_groups:
                haplotype_groups[hp] = []

    # Sort each group by start position (descending; lower positions at top)
    for hp in haplotype_groups:
        haplotype_groups[hp] = sorted(
            haplotype_groups[hp],
            key=lambda read_name: (
                segments_by_read[read_name][0].pos if segments_by_read[read_name] else float("inf")
            ),
            reverse=True,
        )

    # Display order: phased (1,2,3,...) first, then unassigned (0) at bottom
    phased_haplotypes = sorted([hp for hp in haplotype_groups if hp != 0], key=str)
    haplotype_order = phased_haplotypes + ([0] if 0 in haplotype_groups else [])

    # Combine in order
    sorted_read_names = []
    for hp in haplotype_order:
        sorted_read_names.extend(haplotype_groups[hp])

    return sorted_read_names, haplotype_groups, haplotype_order


def calculate_read_positions(
    read_names,
    segments_by_read,
    haplotype_groups,
    haplotype_order=None,
    min_read_height=0.2,
    group_spacing=2.0,
    empty_group_height=0.06,
    empty_group_spacing=0.35,
    multi_alignment_row_spacing=0.20,
    min_read_heights_by_read=None,
):
    """
    Calculate y positions and heights for each read based on number of alignments.
    Also tracks haplotype group boundaries for labeling.

    haplotype_order: if provided, groups are laid out in that order and any haplotype
    with no reads gets a small placeholder slot so it still receives a visible label.

    multi_alignment_row_spacing: vertical gap (in data units) between each sub-row when
    a read has more than one alignment segment. The rendering formula distributes N
    segments into (N+1) equal gaps, so read_height = (N+1) * multi_alignment_row_spacing.
    Single-alignment reads use min_read_height unchanged unless
    min_read_heights_by_read asks for more room.

    Returns:
        tuple: (read_to_y, read_to_y_bottom, read_heights, total_height,
               group_boundaries)
               group_boundaries: haplotype -> (y_start, y_end)
    """
    read_heights = {}
    min_read_heights_by_read = min_read_heights_by_read or {}
    cumulative_y = 0.0
    read_to_y = {}
    read_to_y_bottom = {}
    group_boundaries = {}

    order = (
        haplotype_order
        if haplotype_order
        else list(
            sorted([hp for hp in haplotype_groups if hp != 0])
            + ([0] if 0 in haplotype_groups else [])
        )
    )

    for i, hp in enumerate(order):
        group_reads = haplotype_groups.get(hp, [])
        previous_hp = order[i - 1] if i > 0 else None
        previous_group_empty = (
            previous_hp is not None and len(haplotype_groups.get(previous_hp, [])) == 0
        )

        if i > 0:
            if previous_group_empty or not group_reads:
                cumulative_y += empty_group_spacing
            else:
                cumulative_y += group_spacing

        if not group_reads:
            group_boundaries[hp] = (cumulative_y, cumulative_y + empty_group_height)
            cumulative_y += empty_group_height
        else:
            y_start = cumulative_y
            for read_name in group_reads:
                num_alignments = len(segments_by_read[read_name])
                if num_alignments > 1 and _has_nearby_segments(segments_by_read[read_name]):
                    read_height = (num_alignments + 1) * multi_alignment_row_spacing
                else:
                    read_height = min_read_height
                read_height = max(read_height, min_read_heights_by_read.get(read_name, 0))
                read_heights[read_name] = read_height
                read_to_y[read_name] = cumulative_y + read_height / 2
                read_to_y_bottom[read_name] = cumulative_y
                cumulative_y += read_height
            group_boundaries[hp] = (y_start, cumulative_y)

    total_height = cumulative_y
    return read_to_y, read_to_y_bottom, read_heights, total_height, group_boundaries


def generate_multi_region_filename(region_data_list, output_dir, prefix, filename_regions=None):
    """
    Generate output filename for multiple regions, including all region coordinates.

    Args:
        region_data_list: List of region data dictionaries
        output_dir: Output directory path
        prefix: Optional prefix string
        filename_regions: Optional list of coordinate strings to use instead of region_data_list

    Returns:
        str: Full path to output file
    """
    # Extract region coordinates and sanitize for filename
    region_parts = []
    if filename_regions is not None:
        for coord_str in filename_regions:
            sanitized = coord_str.replace(",", "").replace(":", "_").replace("-", "_")
            region_parts.append(sanitized)
    else:
        for region_data in region_data_list:
            r = region_data["region"]
            coord_str = r.coordinate_str or f"{r.chromosome}_{r.start}_{r.end}"
            # Strip display commas and replace coordinate delimiters for filename safety.
            sanitized = coord_str.replace(",", "").replace(":", "_").replace("-", "_")
            region_parts.append(sanitized)

    # Join all regions with separator
    regions_str = "_".join(region_parts)

    # Build filename
    filename = f"{prefix}_{regions_str}_bokeh.html" if prefix else f"{regions_str}_bokeh.html"

    os.makedirs(output_dir, exist_ok=True)
    return os.path.join(output_dir, filename)
