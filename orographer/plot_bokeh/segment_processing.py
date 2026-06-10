"""Read segment processing and connector data assembly for Bokeh plots."""

import math
from collections.abc import Mapping, Sequence

from ..utils import COMPLEX_SV_REGION_TYPE, ISOSEQ_REGION_TYPE, PARAPHASE_REGION_TYPE
from .utils import PLOT_CONFIG, get_haplotype_color


def get_segment_color(segment, region_type):
    """Color for a segment (YC tag or haplotype/strand; paraphase no YC -> grey)."""
    if hasattr(segment, "color_tag") and segment.color_tag:
        rgb = segment.color_tag.split(",")
        if len(rgb) == 3 and all(channel.isdigit() for channel in rgb):
            r, g, b = (int(channel) for channel in rgb)
            if all(0 <= channel <= 255 for channel in (r, g, b)):
                return f"#{r:02x}{g:02x}{b:02x}"
    # Paraphase: reads without YC tag are grey regardless of strand
    if region_type == PARAPHASE_REGION_TYPE:
        return PLOT_CONFIG["segment_paraphase_color"]
    haplotype = segment.haplotype_tag
    if region_type == ISOSEQ_REGION_TYPE:
        return get_haplotype_color(haplotype)
    if haplotype is None or haplotype == 0:
        return PLOT_CONFIG["segment_unassigned_color"]
    return (
        PLOT_CONFIG["segment_fwd_color"]
        if segment.is_fwd_strand
        else PLOT_CONFIG["segment_rev_color"]
    )


def format_alignment_coordinates(chrom: str, start: int, end: int) -> str:
    """Return alignment coordinates with the exact reference span."""
    span = max(0, end - start)
    return f"{chrom}:{start}-{end} ({span:,} bp)"


def display_block_coordinates(block_start: int, block_end: int) -> tuple[int, int]:
    """Convert a 0-based half-open reference block to displayed coordinates."""
    return block_start + 1, block_end + 1


def get_base_color(base):
    """Get IGV-style color for a nucleotide base."""
    return PLOT_CONFIG["base_colors"].get(base.upper(), PLOT_CONFIG["variant_color_unknown"])


def get_vcf_variant_color(variant):
    """Get color for a VCF variant based on its type."""
    if variant.variant_type == "SNP" and variant.alt_base:
        return get_base_color(variant.alt_base)
    elif variant.variant_type == "INSERTION":
        return PLOT_CONFIG["variant_color_insertion"]
    elif variant.variant_type == "DELETION":
        return PLOT_CONFIG["variant_color_deletion"]
    else:
        return PLOT_CONFIG["variant_color_unknown"]


def connection_key(region_idx, row_index, read_name, segment):
    """Stable key for a visible segment within a multi-region complex-SV plot."""
    return (
        region_idx,
        row_index,
        read_name,
        getattr(segment, "alignment_order", 0),
        getattr(segment, "fwd_read_start", 0),
        getattr(segment, "fwd_read_end", 0),
        segment.chrom,
        segment.pos,
        segment.end,
    )


def segment_overlaps_region(segment, region) -> bool:
    """Return True when *segment* has any visible span in *region*."""
    return segment.pos < region.end and segment.end > region.start


def alignment_identity(segment: object) -> tuple[str, int, int, int, int, bool]:
    """Return a region-independent identity for one alignment segment."""
    return (
        segment.chrom,
        segment.pos,
        segment.end,
        getattr(segment, "fwd_read_start", 0),
        getattr(segment, "fwd_read_end", 0),
        getattr(segment, "is_fwd_strand", True),
    )


def alignment_order_sort_key(
    identity: tuple[str, int, int, int, int, bool],
) -> tuple[int, int, str, int, int, bool]:
    """Sort alignments by original read coordinates, independent of display region."""
    chrom, pos, end, fwd_read_start, fwd_read_end, is_fwd_strand = identity
    return (fwd_read_start, fwd_read_end, chrom, pos, end, is_fwd_strand)


def assign_region_agnostic_alignment_orders(
    region_data_list: Sequence[Mapping[str, object]],
) -> None:
    """Assign read-level alignment numbers across all displayed regions."""
    identities_by_read = {}
    segments_by_identity = {}

    for region_data in region_data_list:
        region = region_data["region"]
        for row_index, bam_row in enumerate(region_data.get("bam_rows", [])):
            for display_read_name, segments in bam_row.get("segments_by_read", {}).items():
                for segment in segments:
                    if not segment_overlaps_region(segment, region):
                        continue
                    read_name = getattr(segment, "readname", display_read_name)
                    read_key = (row_index, read_name)
                    identity = alignment_identity(segment)
                    identities_by_read.setdefault(read_key, set()).add(identity)
                    segments_by_identity.setdefault((read_key, identity), []).append(segment)

    for read_key, identities in identities_by_read.items():
        for alignment_order, identity in enumerate(
            sorted(identities, key=alignment_order_sort_key),
            start=1,
        ):
            for segment in segments_by_identity[(read_key, identity)]:
                segment.alignment_order = alignment_order


def haplotype_label_for_segment(segment) -> str:
    """Return a compact display label for segment haplotype metadata."""
    haplotype = getattr(segment, "haplotype_tag", None)
    return f"HP{haplotype}" if haplotype else "Unassigned"


def connection_id_for_entries(row_index, read_name, source_entry, target_entry) -> str:
    """Return a stable ID shared by both endpoints of one adjacent read link."""
    left_entry, right_entry = sorted(
        [source_entry, target_entry],
        key=lambda entry: (
            entry["alignment_order"],
            entry["fwd_read_start"],
            entry["fwd_read_end"],
            entry["region_idx"],
            entry["pos"],
            entry["end"],
        ),
    )
    return (
        f"row{row_index}:{read_name}:"
        f"a{left_entry['alignment_order']}-a{right_entry['alignment_order']}:"
        f"r{left_entry['region_idx']}-r{right_entry['region_idx']}:"
        f"{left_entry['pos']}-{left_entry['end']}:"
        f"{right_entry['pos']}-{right_entry['end']}"
    )


def build_read_connection_index(region_data_list):
    """Build previous/next read-continuation metadata by segment key."""
    if not region_data_list:
        return {}
    if any(
        bam_row.get("region_type") != COMPLEX_SV_REGION_TYPE
        for region_data in region_data_list
        for bam_row in region_data.get("bam_rows", [])
    ):
        return {}

    visible_by_read = {}
    for region_idx, region_data in enumerate(region_data_list):
        region = region_data["region"]
        for row_index, bam_row in enumerate(region_data.get("bam_rows", [])):
            for display_read_name, segments in bam_row.get("segments_by_read", {}).items():
                for segment in segments:
                    if not segment_overlaps_region(segment, region):
                        continue
                    read_name = getattr(segment, "readname", display_read_name)
                    key = connection_key(region_idx, row_index, display_read_name, segment)
                    visible_by_read.setdefault((row_index, read_name), []).append(
                        {
                            "key": key,
                            "region_idx": region_idx,
                            "row_index": row_index,
                            "read_name": read_name,
                            "alignment_order": getattr(segment, "alignment_order", 0),
                            "fwd_read_start": getattr(segment, "fwd_read_start", 0),
                            "fwd_read_end": getattr(segment, "fwd_read_end", 0),
                            "pos": segment.pos,
                            "end": segment.end,
                            "haplotype_label": haplotype_label_for_segment(segment),
                        }
                    )

    connection_index = {}
    for (row_index, read_name), entries in visible_by_read.items():
        ordered_entries = sorted(
            entries,
            key=lambda entry: (
                entry["alignment_order"],
                entry["fwd_read_start"],
                entry["fwd_read_end"],
                entry["region_idx"],
                entry["pos"],
                entry["end"],
            ),
        )
        for entry_idx, entry in enumerate(ordered_entries):
            if entry_idx > 0:
                previous_entry = ordered_entries[entry_idx - 1]
                connection_id = connection_id_for_entries(
                    row_index, read_name, previous_entry, entry
                )
                connection_index.setdefault(entry["key"], {})["previous"] = (
                    connection_target_metadata(entry, previous_entry, connection_id)
                )
            if entry_idx + 1 < len(ordered_entries):
                next_entry = ordered_entries[entry_idx + 1]
                connection_id = connection_id_for_entries(row_index, read_name, entry, next_entry)
                connection_index.setdefault(entry["key"], {})["next"] = connection_target_metadata(
                    entry, next_entry, connection_id
                )
    return connection_index


def build_read_alignment_summaries(region_data_list):
    """Build rendered alignment coordinate summaries by BAM row and source read."""
    summaries = {}
    if not region_data_list:
        return summaries

    for region_idx, region_data in enumerate(region_data_list):
        region = region_data["region"]
        for row_index, bam_row in enumerate(region_data.get("bam_rows", [])):
            for display_read_name, segments in bam_row.get("segments_by_read", {}).items():
                for segment in segments:
                    if not segment_overlaps_region(segment, region):
                        continue
                    read_name = getattr(segment, "readname", display_read_name)
                    summary_key = (row_index, read_name)
                    summaries.setdefault(summary_key, []).append(
                        {
                            "alignment_number": getattr(segment, "alignment_order", 0),
                            "coordinates": format_alignment_coordinates(
                                segment.chrom,
                                segment.pos,
                                segment.end,
                            ),
                            "region_idx": region_idx,
                            "fwd_read_start": getattr(segment, "fwd_read_start", 0),
                            "fwd_read_end": getattr(segment, "fwd_read_end", 0),
                            "pos": segment.pos,
                            "end": segment.end,
                        }
                    )

    for entries in summaries.values():
        entries.sort(
            key=lambda entry: (
                entry["alignment_number"],
                entry["fwd_read_start"],
                entry["fwd_read_end"],
                entry["region_idx"],
                entry["pos"],
                entry["end"],
            )
        )
    return summaries


def connection_target_metadata(source_entry, target_entry, connection_id):
    """Return metadata describing how one segment end continues to another region."""
    target_region_idx = target_entry["region_idx"]
    source_region_idx = source_entry["region_idx"]
    same_region = target_region_idx == source_region_idx
    target_direction = "same"
    if not same_region:
        target_direction = "left" if target_region_idx < source_region_idx else "right"
    source_haplotype = source_entry["haplotype_label"]
    target_haplotype = target_entry["haplotype_label"]
    return {
        "connection_id": connection_id,
        "connection_scope": "same_region" if same_region else "cross_region",
        "target_key": target_entry["key"],
        "target_region_idx": target_region_idx,
        "target_region_number": target_region_idx + 1,
        "target_alignment_order": target_entry["alignment_order"],
        "target_direction": target_direction,
        "source_haplotype": source_haplotype,
        "target_haplotype": target_haplotype,
        "haplotype_transition": source_haplotype != target_haplotype,
    }


def connection_metadata_for_segment(connection_index, key):
    """Return previous/next connection metadata for *key*, or an empty dict."""
    return connection_index.get(key, {})


def read_filter_flags(alignment_summaries, row_index, read_name) -> tuple[bool, bool]:
    """Return split/chimeric and multi-region flags for one source read."""
    summary_entries = (alignment_summaries or {}).get((row_index, read_name), [])
    region_indices = {entry["region_idx"] for entry in summary_entries}
    return len(summary_entries) > 1, len(region_indices) > 1


def connector_endpoint_y(y_plot, role, has_both_connections):
    """Offset endpoint rings when one segment has connections from both read ends."""
    if not has_both_connections:
        return y_plot
    offset = PLOT_CONFIG["connector_endpoint_y_offset"]
    return y_plot - offset if role == "previous" else y_plot + offset


def double_connected_read_min_height() -> float:
    """Minimum row height needed to contain vertically split connector endpoint rings."""
    return (2 * PLOT_CONFIG["connector_endpoint_y_offset"]) + 0.12


def connector_offset_min_heights_by_read(
    segments_by_read,
    region,
    region_idx,
    row_index,
    connection_index,
) -> dict[str, float]:
    """Return display reads that need extra vertical room for split endpoint rings."""
    if not connection_index:
        return {}
    min_height = double_connected_read_min_height()
    read_min_heights = {}
    for display_read_name, segments in segments_by_read.items():
        for segment in segments:
            if not segment_overlaps_region(segment, region):
                continue
            key = connection_key(region_idx, row_index, display_read_name, segment)
            connections = connection_metadata_for_segment(connection_index, key)
            if "previous" in connections and "next" in connections:
                read_min_heights[display_read_name] = min_height
                break
    return read_min_heights


def add_connector_rows(
    connector_data,
    connections,
    read_name,
    layout_read_name,
    arrow_x0,
    arrow_x1,
    previous_y,
    next_y,
    coordinate_start,
    coordinate_end,
    plot_dom_class,
    plot_model_id,
    has_split_alignment,
    has_multiregion_connection,
) -> None:
    """Append connector marker/stub rows for previous and next read links."""
    if not connections:
        return
    region_span = max(1, coordinate_end - coordinate_start)
    stub_length = max(1, region_span * 0.035)
    for role, endpoint_x in (("previous", arrow_x0), ("next", arrow_x1)):
        metadata = connections.get(role)
        if not metadata:
            continue
        if metadata.get("connection_scope") == "same_region":
            continue
        endpoint_y = previous_y if role == "previous" else next_y
        target_direction = metadata["target_direction"]
        if target_direction == "left":
            stub_x1 = max(coordinate_start, endpoint_x - stub_length)
            marker_angle = math.pi / 2
        elif target_direction == "right":
            stub_x1 = min(coordinate_end, endpoint_x + stub_length)
            marker_angle = -math.pi / 2
        else:
            stub_x1 = endpoint_x
            marker_angle = 0
        tooltip = (
            f"{role}: alignment {metadata['target_alignment_order']} "
            f"in region {metadata['target_region_number']}"
        )
        if metadata["haplotype_transition"]:
            tooltip += f" {metadata['source_haplotype']} to {metadata['target_haplotype']}"
        connector_data["stub_x0"].append(endpoint_x)
        connector_data["stub_x1"].append(stub_x1)
        connector_data["y"].append(endpoint_y)
        connector_data["marker_x"].append(endpoint_x)
        connector_data["marker_angle"].append(marker_angle)
        connector_data["read_name"].append(read_name)
        connector_data["layout_read_name"].append(layout_read_name)
        connector_data["source_kind"].append("connector")
        connector_data["tooltip"].append(tooltip)
        connector_data["target_region"].append(metadata["target_region_number"])
        connector_data["target_alignment"].append(metadata["target_alignment_order"])
        connector_data["connection_role"].append(role)
        connector_data["connection_id"].append(metadata["connection_id"])
        connector_data["connection_scope"].append(metadata.get("connection_scope", "cross_region"))
        connector_data["source_haplotype"].append(metadata["source_haplotype"])
        connector_data["target_haplotype"].append(metadata["target_haplotype"])
        connector_data["haplotype_transition"].append(metadata["haplotype_transition"])
        connector_data["plot_dom_class"].append(plot_dom_class)
        connector_data["plot_model_id"].append(plot_model_id)
        connector_data["has_split_alignment"].append(has_split_alignment)
        connector_data["has_multiregion_connection"].append(has_multiregion_connection)
        connector_data["read_filter_alpha"].append(1.0)
        connector_data["connector_line_alpha"].append(0.9)
        connector_data["connector_selected_alpha"].append(1.0)
        connector_data["connector_nonselected_alpha"].append(
            PLOT_CONFIG["connector_nonselection_alpha"]
        )


def empty_same_region_connector_data():
    """Return the data container used for same-panel read-continuation lines."""
    return {
        "xs": [],
        "ys": [],
        "read_name": [],
        "layout_read_name": [],
        "source_kind": [],
        "connection_id": [],
        "has_split_alignment": [],
        "has_multiregion_connection": [],
        "read_filter_alpha": [],
        "connector_line_alpha": [],
        "connector_selected_alpha": [],
        "connector_nonselected_alpha": [],
    }


def connector_lane_spacing() -> float:
    """Return the minimum visual spacing between routed connector lanes."""
    endpoint_offset = PLOT_CONFIG["connector_endpoint_y_offset"]
    lane_offset = PLOT_CONFIG["connector_lane_offset"]
    return max(lane_offset * 2, endpoint_offset + 0.02)


def intervals_overlap(
    first_min: float,
    first_max: float,
    second_min: float,
    second_max: float,
) -> bool:
    """Return True when two horizontal connector intervals overlap."""
    return first_min <= second_max and second_min <= first_max


def point_between(value: float, first: float, second: float) -> bool:
    """Return True when value lies between first and second, inclusive."""
    lower = min(first, second)
    upper = max(first, second)
    return lower <= value <= upper


def connector_candidate_lanes(request) -> list[float]:
    """Return candidate y-lanes for one same-region connector request."""
    clearance = PLOT_CONFIG["connector_endpoint_y_offset"] + 0.08
    lane_spacing = connector_lane_spacing()
    source_y = request["source_y"]
    target_y = request["target_y"]
    lane_min = request["lane_min"]
    lane_max = request["lane_max"]
    outer_lane_min = request["outer_lane_min"]
    outer_lane_max = request["outer_lane_max"]
    above = max(source_y, target_y) + clearance
    below = min(source_y, target_y) - clearance
    candidates = []

    if lane_min <= lane_max:
        candidates.extend(
            [
                min(above, lane_max),
                max(below, lane_min),
                lane_max,
                lane_min,
            ]
        )
    else:
        clamped_above = min(max(above, outer_lane_min), outer_lane_max)
        clamped_below = min(max(below, outer_lane_min), outer_lane_max)
        candidates.extend(
            [
                clamped_above,
                clamped_below,
                source_y,
                target_y,
                (source_y + target_y) / 2,
            ]
        )
        for step in range(1, 5):
            high_lane = min(clamped_above + lane_spacing * step, outer_lane_max)
            low_lane = max(clamped_below - lane_spacing * step, outer_lane_min)
            candidates.append(high_lane)
            candidates.append(low_lane)

    unique_candidates = []
    for candidate in candidates:
        if candidate not in unique_candidates:
            unique_candidates.append(candidate)
    return unique_candidates


def connector_lane_cost(request, candidate_y: float, placed_routes) -> float:
    """Score a candidate connector lane; lower cost means less visual interference."""
    lane_spacing = connector_lane_spacing()
    x_min = request["x_min"]
    x_max = request["x_max"]
    cost = abs(candidate_y - request["source_y"]) + abs(candidate_y - request["target_y"])

    for placed in placed_routes:
        placed_y = placed["lane_y"]
        same_lane = abs(candidate_y - placed_y) <= lane_spacing * 0.5
        route_overlap = intervals_overlap(x_min, x_max, placed["x_min"], placed["x_max"])
        if same_lane and route_overlap:
            cost += 500
        if route_overlap:
            proximity = max(0.0, lane_spacing - abs(candidate_y - placed_y))
            cost += proximity * 80
        if point_between(placed_y, request["source_y"], candidate_y) and (
            placed["x_min"] <= request["source_x"] <= placed["x_max"]
        ):
            cost += 100
        if point_between(placed_y, request["target_y"], candidate_y) and (
            placed["x_min"] <= request["target_x"] <= placed["x_max"]
        ):
            cost += 100

    return cost


def route_same_region_connectors(requests) -> list[float]:
    """Choose connector lanes with interval packing and simple crossing avoidance."""
    assignments = [0.0] * len(requests)
    placed_routes = []
    routing_order = sorted(
        requests,
        key=lambda request: (-request["span"], request["index"]),
    )

    for request in routing_order:
        candidates = connector_candidate_lanes(request)
        lane_y = min(
            candidates,
            key=lambda candidate_y: connector_lane_cost(
                request,
                candidate_y,
                placed_routes,
            ),
        )
        assignments[request["index"]] = lane_y
        placed_routes.append(
            {
                "x_min": request["x_min"],
                "x_max": request["x_max"],
                "lane_y": lane_y,
            }
        )

    return assignments


def add_same_region_connector_lines(line_data, visible_segment_endpoints, connection_index) -> None:
    """Append dashed connector paths between adjacent same-region read segments."""
    if not connection_index:
        return
    requests = []
    for source_key, source_endpoints in visible_segment_endpoints.items():
        metadata = connection_index.get(source_key, {}).get("next")
        if not metadata or metadata.get("connection_scope") != "same_region":
            continue
        target_endpoints = visible_segment_endpoints.get(metadata.get("target_key"))
        if not target_endpoints:
            continue
        source_x = source_endpoints["next_x"]
        source_y = source_endpoints["next_y"]
        target_x = target_endpoints["previous_x"]
        target_y = target_endpoints["previous_y"]
        source_lane_min = source_endpoints["read_y_min"]
        source_lane_max = source_endpoints["read_y_max"]
        target_lane_min = target_endpoints["read_y_min"]
        target_lane_max = target_endpoints["read_y_max"]
        lane_min = max(source_lane_min, target_lane_min)
        lane_max = min(source_lane_max, target_lane_max)
        x_min = min(source_x, target_x)
        x_max = max(source_x, target_x)
        requests.append(
            {
                "index": len(requests),
                "source_x": source_x,
                "source_y": source_y,
                "target_x": target_x,
                "target_y": target_y,
                "x_min": x_min,
                "x_max": x_max,
                "span": x_max - x_min,
                "lane_min": lane_min,
                "lane_max": lane_max,
                "outer_lane_min": min(source_lane_min, target_lane_min),
                "outer_lane_max": max(source_lane_max, target_lane_max),
                "read_name": source_endpoints["read_name"],
                "layout_read_name": source_endpoints.get(
                    "layout_read_name",
                    source_endpoints["read_name"],
                ),
                "connection_id": metadata["connection_id"],
                "has_split_alignment": source_endpoints.get("has_split_alignment", False),
                "has_multiregion_connection": source_endpoints.get(
                    "has_multiregion_connection",
                    False,
                ),
            }
        )

    lane_assignments = route_same_region_connectors(requests)
    for request, lane_y in zip(requests, lane_assignments, strict=True):
        source_x = request["source_x"]
        source_y = request["source_y"]
        target_x = request["target_x"]
        target_y = request["target_y"]
        line_data["xs"].append([source_x, source_x, target_x, target_x])
        line_data["ys"].append([source_y, lane_y, lane_y, target_y])
        line_data["read_name"].append(request["read_name"])
        line_data["layout_read_name"].append(request["layout_read_name"])
        line_data["source_kind"].append("connector")
        line_data["connection_id"].append(request["connection_id"])
        line_data["has_split_alignment"].append(request["has_split_alignment"])
        line_data["has_multiregion_connection"].append(request["has_multiregion_connection"])
        line_data["read_filter_alpha"].append(1.0)
        line_data["connector_line_alpha"].append(PLOT_CONFIG["connector_line_alpha"])
        line_data["connector_selected_alpha"].append(PLOT_CONFIG["connector_selection_line_alpha"])
        line_data["connector_nonselected_alpha"].append(PLOT_CONFIG["connector_nonselection_alpha"])


def process_segments(
    segments_by_read,
    read_names,
    read_to_y,
    read_heights,
    coordinate_start,
    coordinate_end,
    region_type,
    sample_label=None,
    region_idx=0,
    row_index=0,
    connection_index=None,
    alignment_summaries=None,
    plot_dom_class="",
    plot_model_id="",
):
    """Process segments and extract arrow, clickable label, and variant data."""
    arrow_x0_list, arrow_x1_list, arrow_y_list, arrow_color_list, arrow_read_names = (
        [],
        [],
        [],
        [],
        [],
    )
    arrow_segment_ids, arrow_region_indices, arrow_row_indices, arrow_alignment_orders = (
        [],
        [],
        [],
        [],
    )
    arrow_fwd_read_starts, arrow_fwd_read_ends, arrow_source_kinds = [], [], []
    arrow_layout_read_names = []
    arrow_has_split, arrow_has_multiregion, arrow_filter_alpha = [], [], []
    arrow_line_alphas, arrowhead_alphas = [], []
    arrow_selected_alphas, arrow_nonselected_alphas = [], []
    arrow_y0_list, arrow_y1_list = [], []
    same_region_connector_data = empty_same_region_connector_data()
    visible_segment_endpoints = {}
    connector_data = {
        "stub_x0": [],
        "stub_x1": [],
        "y": [],
        "marker_x": [],
        "marker_angle": [],
        "read_name": [],
        "layout_read_name": [],
        "source_kind": [],
        "tooltip": [],
        "target_region": [],
        "target_alignment": [],
        "connection_role": [],
        "connection_id": [],
        "connection_scope": [],
        "source_haplotype": [],
        "target_haplotype": [],
        "haplotype_transition": [],
        "plot_dom_class": [],
        "plot_model_id": [],
        "has_split_alignment": [],
        "has_multiregion_connection": [],
        "read_filter_alpha": [],
        "connector_line_alpha": [],
        "connector_selected_alpha": [],
        "connector_nonselected_alpha": [],
    }
    clickable_x, clickable_y, clickable_customdata = [], [], []
    mismatch_x, mismatch_y, mismatch_alt, mismatch_color = [], [], [], []
    mismatch_read_names = []
    mismatch_layout_read_names = []
    mismatch_split, mismatch_multiregion, mismatch_filter_alpha = [], [], []
    mismatch_fill_alpha, mismatch_line_alpha, mismatch_text_alpha = [], [], []
    insertion_x, insertion_y, insertion_size, insertion_count, insertion_is_1bp = (
        [],
        [],
        [],
        [],
        [],
    )
    insertion_split, insertion_multiregion, insertion_filter_alpha = [], [], []
    insertion_read_names = []
    insertion_layout_read_names = []
    insertion_fill_alpha, insertion_line_alpha, insertion_text_alpha = [], [], []
    deletion_x0, deletion_x1, deletion_y, deletion_is_1bp = [], [], [], []
    deletion_read_names = []
    deletion_layout_read_names = []
    deletion_split, deletion_multiregion, deletion_filter_alpha = [], [], []
    deletion_line_alpha = []

    for read_name in read_names:
        segments = segments_by_read[read_name]
        y_pos = read_to_y[read_name]
        read_height = read_heights[read_name]
        num_segments = len(segments)
        for idx, segment in enumerate(segments):
            has_split_alignment, has_multiregion_connection = read_filter_flags(
                alignment_summaries,
                row_index,
                segment.readname,
            )
            segment_blocks = [(segment.pos, segment.end)]
            is_isoseq = region_type == ISOSEQ_REGION_TYPE
            if is_isoseq:
                segment_blocks = getattr(segment, "aligned_blocks", None) or segment_blocks
            visible_blocks = []
            for block_index, (block_start, block_end) in enumerate(segment_blocks):
                if is_isoseq:
                    display_start, display_end = display_block_coordinates(block_start, block_end)
                else:
                    display_start, display_end = block_start, block_end
                if display_start > coordinate_end or display_end < coordinate_start:
                    continue
                visible_blocks.append(
                    (
                        block_index,
                        max(display_start, coordinate_start),
                        min(display_end, coordinate_end),
                    )
                )
            if visible_blocks:
                block_count = len(segment_blocks)
                block_number_by_index = [
                    block_index + 1 if segment.is_fwd_strand else block_count - block_index
                    for block_index in range(block_count)
                ]
                block_summary = sorted(
                    [
                        (
                            block_number_by_index[block_index],
                            format_alignment_coordinates(
                                segment.chrom,
                                *display_block_coordinates(block_start, block_end),
                            ),
                        )
                        for block_index, (block_start, block_end) in enumerate(segment_blocks)
                    ],
                    key=lambda item: item[0],
                )
                block_summary_numbers = [item[0] for item in block_summary]
                block_summary_coordinates = [item[1] for item in block_summary]
                segment_start = min(start for _block_index, start, _end in visible_blocks)
                segment_end = max(end for _block_index, _start, end in visible_blocks)
                color = get_segment_color(segment, region_type=region_type)
                if num_segments > 1:
                    spacing = read_height / (num_segments + 1)
                    y_offset = -read_height / 2 + spacing * (idx + 1)
                else:
                    y_offset = 0
                y_plot = y_pos + y_offset
                mid_x = (segment_start + segment_end) / 2
                alignment_number = getattr(segment, "alignment_order", 0) or (idx + 1)
                key = connection_key(region_idx, row_index, read_name, segment)
                connections = connection_metadata_for_segment(connection_index or {}, key)
                has_both_connections = "previous" in connections and "next" in connections
                arrow_y0 = connector_endpoint_y(y_plot, "previous", has_both_connections)
                arrow_y1 = connector_endpoint_y(y_plot, "next", has_both_connections)
                segment_id = (
                    f"r{region_idx}:row{row_index}:{segment.readname}:"
                    f"{alignment_number}:{segment.chrom}:{segment.pos}-{segment.end}"
                )
                for block_index, plot_start, plot_end in visible_blocks:
                    display_alignment_number = (
                        block_number_by_index[block_index] if is_isoseq else alignment_number
                    )
                    if segment.is_fwd_strand:
                        arrow_x0, arrow_x1 = plot_start, plot_end
                    else:
                        arrow_x0, arrow_x1 = plot_end, plot_start
                    if abs(arrow_x1 - arrow_x0) <= 0:
                        continue
                    arrow_x0_list.append(arrow_x0)
                    arrow_x1_list.append(arrow_x1)
                    arrow_y_list.append(y_plot)
                    arrow_y0_list.append(arrow_y0)
                    arrow_y1_list.append(arrow_y1)
                    arrow_color_list.append(color)
                    arrow_read_names.append(segment.readname)
                    arrow_layout_read_names.append(read_name)
                    arrow_segment_ids.append(f"{segment_id}:block{display_alignment_number}")
                    arrow_region_indices.append(region_idx)
                    arrow_row_indices.append(row_index)
                    arrow_alignment_orders.append(display_alignment_number)
                    arrow_fwd_read_starts.append(getattr(segment, "fwd_read_start", 0))
                    arrow_fwd_read_ends.append(getattr(segment, "fwd_read_end", 0))
                    arrow_source_kinds.append("arrow")
                    arrow_has_split.append(has_split_alignment)
                    arrow_has_multiregion.append(has_multiregion_connection)
                    arrow_filter_alpha.append(1.0)
                    arrow_line_alphas.append(PLOT_CONFIG["arrow_line_alpha"])
                    arrowhead_alphas.append(PLOT_CONFIG["arrowhead_fill_alpha"])
                    arrow_selected_alphas.append(1.0)
                    arrow_nonselected_alphas.append(PLOT_CONFIG["arrow_nonselection_line_alpha"])
                    if is_isoseq:
                        strand_str = "Forward (+)" if segment.is_fwd_strand else "Reverse (-)"
                        haplotype_str = (
                            f"HP:{segment.haplotype_tag}" if segment.haplotype_tag else "Unassigned"
                        )
                        clickable_x.append((plot_start + plot_end) / 2)
                        clickable_y.append(y_plot)
                        clickable_customdata.append(
                            {
                                "read_name": segment.readname,
                                "layout_read_name": read_name,
                                "alignment_number": display_alignment_number,
                                "strand": strand_str,
                                "coordinates": format_alignment_coordinates(
                                    segment.chrom,
                                    plot_start,
                                    plot_end,
                                ),
                                "haplotype": haplotype_str,
                                "sample_label": sample_label or "",
                                "inclusion_reason": "",
                                "all_alignment_numbers": block_summary_numbers,
                                "all_alignment_coordinates": block_summary_coordinates,
                                "has_split_alignment": has_split_alignment,
                                "has_multiregion_connection": has_multiregion_connection,
                            }
                        )
                if visible_blocks:
                    previous_x = segment_start if segment.is_fwd_strand else segment_end
                    next_x = segment_end if segment.is_fwd_strand else segment_start
                    visible_segment_endpoints[key] = {
                        "previous_x": previous_x,
                        "previous_y": y_plot,
                        "next_x": next_x,
                        "next_y": y_plot,
                        "read_y_min": y_pos - (read_height / 2),
                        "read_y_max": y_pos + (read_height / 2),
                        "read_name": segment.readname,
                        "layout_read_name": read_name,
                        "has_split_alignment": has_split_alignment,
                        "has_multiregion_connection": has_multiregion_connection,
                    }
                    add_connector_rows(
                        connector_data,
                        connections,
                        segment.readname,
                        read_name,
                        previous_x,
                        next_x,
                        arrow_y0,
                        arrow_y1,
                        coordinate_start,
                        coordinate_end,
                        plot_dom_class,
                        plot_model_id,
                        has_split_alignment,
                        has_multiregion_connection,
                    )

                if hasattr(segment, "mismatches") and segment.mismatches:
                    for ref_pos, alt_base in segment.mismatches:
                        ref_pos_1based = ref_pos + 1
                        if coordinate_start <= ref_pos_1based <= coordinate_end:
                            mismatch_x.append(ref_pos_1based)
                            mismatch_y.append(y_plot)
                            mismatch_alt.append(alt_base)
                            mismatch_color.append(get_base_color(alt_base))
                            mismatch_read_names.append(segment.readname)
                            mismatch_layout_read_names.append(read_name)
                            mismatch_split.append(has_split_alignment)
                            mismatch_multiregion.append(has_multiregion_connection)
                            mismatch_filter_alpha.append(1.0)
                            mismatch_fill_alpha.append(PLOT_CONFIG["mismatch_fill_alpha"])
                            mismatch_line_alpha.append(PLOT_CONFIG["mismatch_line_alpha"])
                            mismatch_text_alpha.append(1.0)
                if hasattr(segment, "insertions") and segment.insertions:
                    for ref_pos, inserted_bases in segment.insertions:
                        ref_pos_1based = ref_pos + 1
                        if coordinate_start <= ref_pos_1based <= coordinate_end:
                            insertion_x.append(ref_pos_1based)
                            insertion_y.append(y_plot)
                            insertion_size.append(
                                min(
                                    PLOT_CONFIG["insertion_size_max"],
                                    max(
                                        PLOT_CONFIG["insertion_size_min"],
                                        len(inserted_bases),
                                    ),
                                )
                            )
                            insertion_count.append(len(inserted_bases))
                            insertion_is_1bp.append(len(inserted_bases) == 1)
                            insertion_read_names.append(segment.readname)
                            insertion_layout_read_names.append(read_name)
                            insertion_split.append(has_split_alignment)
                            insertion_multiregion.append(has_multiregion_connection)
                            insertion_filter_alpha.append(1.0)
                            insertion_fill_alpha.append(PLOT_CONFIG["insertion_fill_alpha"])
                            insertion_line_alpha.append(PLOT_CONFIG["insertion_line_alpha"])
                            insertion_text_alpha.append(1.0)
                if hasattr(segment, "deletions") and segment.deletions:
                    for del_start, del_end in segment.deletions:
                        del_start_1based = del_start + 1
                        del_end_1based = del_end - 1
                        if (
                            del_start_1based <= coordinate_end
                            and del_end_1based >= coordinate_start
                        ):
                            visible_start = max(del_start_1based, coordinate_start)
                            visible_end = min(del_end_1based, coordinate_end)
                            deletion_x0.append(visible_start)
                            deletion_x1.append(visible_end)
                            deletion_y.append(y_plot)
                            deletion_is_1bp.append((del_end - del_start) == 1)
                            deletion_read_names.append(segment.readname)
                            deletion_layout_read_names.append(read_name)
                            deletion_split.append(has_split_alignment)
                            deletion_multiregion.append(has_multiregion_connection)
                            deletion_filter_alpha.append(1.0)
                            deletion_line_alpha.append(PLOT_CONFIG["deletion_line_alpha"])

                strand_str = "Forward (+)" if segment.is_fwd_strand else "Reverse (-)"
                haplotype_str = (
                    f"HP:{segment.haplotype_tag}" if segment.haplotype_tag else "Unassigned"
                )
                inclusion_reason = ""
                if region_type == COMPLEX_SV_REGION_TYPE:
                    inclusion_reason = getattr(segment, "inclusion_reason", "") or ""
                if not is_isoseq:
                    clickable_x.append(mid_x)
                    clickable_y.append(y_plot)
                    alignment_number = getattr(segment, "alignment_order", 0) or (idx + 1)
                    summary_entries = (alignment_summaries or {}).get(
                        (row_index, segment.readname),
                        [],
                    )
                    clickable_customdata.append(
                        {
                            "read_name": segment.readname,
                            "layout_read_name": read_name,
                            "alignment_number": alignment_number,
                            "strand": strand_str,
                            "coordinates": format_alignment_coordinates(
                                segment.chrom,
                                segment.pos,
                                segment.end,
                            ),
                            "haplotype": haplotype_str,
                            "sample_label": sample_label or "",
                            "inclusion_reason": inclusion_reason,
                            "all_alignment_numbers": [
                                entry["alignment_number"] for entry in summary_entries
                            ],
                            "all_alignment_coordinates": [
                                entry["coordinates"] for entry in summary_entries
                            ],
                            "has_split_alignment": has_split_alignment,
                            "has_multiregion_connection": has_multiregion_connection,
                        }
                    )

    add_same_region_connector_lines(
        same_region_connector_data,
        visible_segment_endpoints,
        connection_index or {},
    )

    return (
        {
            "x0": arrow_x0_list,
            "x1": arrow_x1_list,
            "y": arrow_y_list,
            "y0": arrow_y0_list,
            "y1": arrow_y1_list,
            "color": arrow_color_list,
            "read_name": arrow_read_names,
            "layout_read_name": arrow_layout_read_names,
            "segment_id": arrow_segment_ids,
            "region_idx": arrow_region_indices,
            "row_index": arrow_row_indices,
            "alignment_order": arrow_alignment_orders,
            "fwd_read_start": arrow_fwd_read_starts,
            "fwd_read_end": arrow_fwd_read_ends,
            "source_kind": arrow_source_kinds,
            "has_split_alignment": arrow_has_split,
            "has_multiregion_connection": arrow_has_multiregion,
            "read_filter_alpha": arrow_filter_alpha,
            "arrow_line_alpha": arrow_line_alphas,
            "arrowhead_alpha": arrowhead_alphas,
            "arrow_selected_alpha": arrow_selected_alphas,
            "arrow_nonselected_alpha": arrow_nonselected_alphas,
        },
        {"x": clickable_x, "y": clickable_y, "customdata": clickable_customdata},
        {
            "mismatch": {
                "x": mismatch_x,
                "y": mismatch_y,
                "alt": mismatch_alt,
                "color": mismatch_color,
                "read_name": mismatch_read_names,
                "layout_read_name": mismatch_layout_read_names,
                "has_split_alignment": mismatch_split,
                "has_multiregion_connection": mismatch_multiregion,
                "read_filter_alpha": mismatch_filter_alpha,
                "marker_fill_alpha": mismatch_fill_alpha,
                "marker_line_alpha": mismatch_line_alpha,
                "text_alpha": mismatch_text_alpha,
            },
            "insertion": {
                "x": insertion_x,
                "y": insertion_y,
                "size": insertion_size,
                "count": insertion_count,
                "is_1bp": insertion_is_1bp,
                "read_name": insertion_read_names,
                "layout_read_name": insertion_layout_read_names,
                "has_split_alignment": insertion_split,
                "has_multiregion_connection": insertion_multiregion,
                "read_filter_alpha": insertion_filter_alpha,
                "marker_fill_alpha": insertion_fill_alpha,
                "marker_line_alpha": insertion_line_alpha,
                "text_alpha": insertion_text_alpha,
            },
            "deletion": {
                "x0": deletion_x0,
                "x1": deletion_x1,
                "y": deletion_y,
                "is_1bp": deletion_is_1bp,
                "read_name": deletion_read_names,
                "layout_read_name": deletion_layout_read_names,
                "has_split_alignment": deletion_split,
                "has_multiregion_connection": deletion_multiregion,
                "read_filter_alpha": deletion_filter_alpha,
                "line_alpha": deletion_line_alpha,
            },
        },
        connector_data,
        same_region_connector_data,
    )
