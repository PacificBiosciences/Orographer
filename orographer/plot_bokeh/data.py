"""Data and layout helpers for Bokeh plotting (read ordering, positions, filenames)."""

import os


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


def get_read_haplotype(segments):
    """Get the haplotype tag from the first segment of a read.
    Treats None, 0, and `Unassigned` as unassigned (group 0)."""
    if not segments:
        return 0
    hp = segments[0].haplotype_tag
    if hp is None or hp == 0 or str(hp).lower() == "unassigned":
        return 0
    return hp


def sort_read_names(segments_by_read):
    """
    Sort read names by haplotype group, then by start position.
    Dynamically detects all haplotypes present in the data.

    Returns:
        tuple: (sorted_read_names, haplotype_groups, haplotype_order) where:
               - haplotype_groups: haplotype -> list of read names
               - haplotype_order: display order (phased first, then unassigned)
    """
    # First pass: discover all unique haplotypes
    haplotype_groups = {}

    for read_name in segments_by_read:
        haplotype = get_read_haplotype(segments_by_read[read_name])
        if haplotype not in haplotype_groups:
            haplotype_groups[haplotype] = []
        haplotype_groups[haplotype].append(read_name)

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
    phased_haplotypes = sorted([hp for hp in haplotype_groups if hp != 0])
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
    spacing_per_alignment=0.10,
    min_read_height=0.2,
    group_spacing=2.0,
):
    """
    Calculate y positions and heights for each read based on number of alignments.
    Also tracks haplotype group boundaries for labeling.

    Returns:
        tuple: (read_to_y, read_to_y_bottom, read_heights, total_height,
               group_boundaries)
               group_boundaries: haplotype -> (y_start, y_end)
    """
    alignments_per_read = {read_name: len(segments_by_read[read_name]) for read_name in read_names}
    read_heights = {}
    cumulative_y = 0
    read_to_y = {}
    read_to_y_bottom = {}

    # Track group boundaries
    group_boundaries = {}
    current_group = None

    for read_name in read_names:
        haplotype = get_read_haplotype(segments_by_read[read_name])

        if current_group != haplotype:
            if current_group is not None and current_group in group_boundaries:
                group_boundaries[current_group] = (
                    group_boundaries[current_group][0],
                    cumulative_y,
                )

            if current_group is not None and haplotype_groups.get(haplotype):
                cumulative_y += group_spacing

            if haplotype_groups.get(haplotype):
                group_boundaries[haplotype] = (cumulative_y, cumulative_y)
                current_group = haplotype

        num_alignments = alignments_per_read[read_name]
        read_height = max(min_read_height, num_alignments * spacing_per_alignment)
        read_heights[read_name] = read_height
        read_to_y[read_name] = cumulative_y + read_height / 2
        read_to_y_bottom[read_name] = cumulative_y
        cumulative_y += read_height

    if current_group is not None and current_group in group_boundaries:
        group_boundaries[current_group] = (
            group_boundaries[current_group][0],
            cumulative_y,
        )

    total_height = cumulative_y
    return read_to_y, read_to_y_bottom, read_heights, total_height, group_boundaries


def generate_multi_region_filename(region_data_list, output_dir, prefix):
    """
    Generate output filename for multiple regions, including all region coordinates.

    Args:
        region_data_list: List of region data dictionaries
        output_dir: Output directory path
        prefix: Optional prefix string

    Returns:
        str: Full path to output file
    """
    # Extract region coordinates and sanitize for filename
    region_parts = []
    for region_data in region_data_list:
        coord_str = region_data["region"].coordinate_str
        # Replace colons and hyphens with underscores for filename safety
        sanitized = coord_str.replace(":", "_").replace("-", "_")
        region_parts.append(sanitized)

    # Join all regions with separator
    regions_str = "_".join(region_parts)

    # Build filename
    filename = f"{prefix}_{regions_str}_bokeh.html" if prefix else f"{regions_str}_bokeh.html"

    os.makedirs(output_dir, exist_ok=True)
    return os.path.join(output_dir, filename)
