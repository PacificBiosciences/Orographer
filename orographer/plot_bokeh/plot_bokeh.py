"""
Orchestration and track content for Bokeh plots
(colors, segments, arrows, variants, entry point).
"""

import json
import logging
import math
import os
import uuid
from collections.abc import Mapping, Sequence

from bokeh.layouts import column, row
from bokeh.models import Button, ColumnDataSource, CustomJS, Div, Select, Spacer, TapTool

from ..utils import COMPLEX_SV_REGION_TYPE, ISOSEQ_REGION_TYPE, PARAPHASE_REGION_TYPE
from . import coverage as coverage_data
from . import isoseq as isoseq_data
from . import isoseq_store, isoseq_tracks, segment_processing
from . import renderers as plot_renderers
from . import sources as source_data
from .callbacks import (
    get_arrow_tap_callback,
    get_arrow_tap_callback_multi_region,
    get_number_click_callback,
    get_vcf_variant_click_callback,
    save_plot_with_modal,
)
from .data import (
    READ_FILTER_LAYOUT_MODES,
    ReadFilterLayoutRequest,
    add_read_filter_visibility_columns,
    build_read_filter_layouts,
    calculate_read_positions,
    generate_multi_region_filename,
    sort_read_names,
)
from .figures import (
    _build_gene_track_raw_data,
    _build_hp_label_raw_data,
    _build_repeat_density_source_data,
    _build_separator_raw_data,
    _insertion_raw_source_data,
    add_coverage_to_figure,
    add_cursor_guide_callbacks,
    add_cursor_guide_to_figures,
    add_gene_track,
    add_haplotype_labels,
    add_insertion_markers_to_figure,
    add_separator_lines,
    build_coverage_source_data,
    create_bokeh_figure,
    create_bokeh_figure_shared_x,
    create_coordinate_display,
    create_coverage_figure,
    create_dotplot_thumbnail,
    create_gene_track_figure,
    create_genomic_x_axis_strip,
    create_global_checkbox_controls,
    create_repeat_density_figure,
    create_vcf_track_figure,
    format_region_size,
)
from .lazy_store import write_json_gzip
from .utils import (
    PLOT_CONFIG,
    RegionBuildState,
    build_haplotype_color_map,
    get_haplotype_color,
    load_javascript,
)

logger = logging.getLogger(__name__)
logging.getLogger("bokeh").setLevel(logging.WARNING)

ISOSEQ_READ_PAGE_SIZE = isoseq_data.ISOSEQ_READ_PAGE_SIZE
ISOSEQ_READ_SHARD_SIZE = isoseq_data.ISOSEQ_READ_SHARD_SIZE
add_isoseq_transcripts_to_plot = isoseq_tracks.add_isoseq_transcripts_to_plot
create_isoseq_controls = isoseq_tracks.create_isoseq_controls
_build_isoseq_bam_row_track = isoseq_tracks.build_isoseq_bam_row_track

get_segment_color = segment_processing.get_segment_color
format_alignment_coordinates = segment_processing.format_alignment_coordinates
get_base_color = segment_processing.get_base_color
get_vcf_variant_color = segment_processing.get_vcf_variant_color
connection_key = segment_processing.connection_key
segment_overlaps_region = segment_processing.segment_overlaps_region
alignment_identity = segment_processing.alignment_identity
alignment_order_sort_key = segment_processing.alignment_order_sort_key
assign_region_agnostic_alignment_orders = segment_processing.assign_region_agnostic_alignment_orders
haplotype_label_for_segment = segment_processing.haplotype_label_for_segment
connection_id_for_entries = segment_processing.connection_id_for_entries
build_read_connection_index = segment_processing.build_read_connection_index
build_read_alignment_summaries = segment_processing.build_read_alignment_summaries
connection_target_metadata = segment_processing.connection_target_metadata
connection_metadata_for_segment = segment_processing.connection_metadata_for_segment
read_filter_flags = segment_processing.read_filter_flags
connector_endpoint_y = segment_processing.connector_endpoint_y
double_connected_read_min_height = segment_processing.double_connected_read_min_height
connector_offset_min_heights_by_read = segment_processing.connector_offset_min_heights_by_read
add_connector_rows = segment_processing.add_connector_rows
empty_same_region_connector_data = segment_processing.empty_same_region_connector_data
connector_lane_spacing = segment_processing.connector_lane_spacing
intervals_overlap = segment_processing.intervals_overlap
point_between = segment_processing.point_between
connector_candidate_lanes = segment_processing.connector_candidate_lanes
connector_lane_cost = segment_processing.connector_lane_cost
route_same_region_connectors = segment_processing.route_same_region_connectors
add_same_region_connector_lines = segment_processing.add_same_region_connector_lines
process_segments = segment_processing.process_segments


def _layout_mode_column(column_name: str, mode: str) -> str:
    """Return the alternate source-column name for a read filter layout mode."""
    return f"{column_name}_{mode}"


def _read_filter_flags_by_layout_read(
    segments_by_read: Mapping[str, Sequence[object]],
    row_index: int,
    alignment_summaries: Mapping[tuple[int, str], Sequence[Mapping[str, object]]] | None,
) -> dict[str, tuple[bool, bool]]:
    """Return split and multiregion evidence flags keyed by displayed read name."""
    flags_by_read: dict[str, tuple[bool, bool]] = {}
    for display_read_name, segments in segments_by_read.items():
        has_split = False
        has_multiregion = False
        for segment in segments:
            segment_read_name = getattr(segment, "readname", display_read_name)
            segment_split, segment_multiregion = read_filter_flags(
                alignment_summaries,
                row_index,
                str(segment_read_name),
            )
            has_split = has_split or segment_split
            has_multiregion = has_multiregion or segment_multiregion
        flags_by_read[display_read_name] = (has_split, has_multiregion)
    return flags_by_read


def _build_read_filter_layout_modes(
    request: ReadFilterLayoutRequest,
) -> dict[str, dict[str, object]]:
    """Build compact read layouts for all checkbox combinations."""
    return build_read_filter_layouts(request)


def _layout_mode_delta(
    read_name: str,
    mode: str,
    layout_modes: Mapping[str, Mapping[str, object]],
) -> float | None:
    """Return the y shift for a read in a layout mode, or None if hidden."""
    all_layout = layout_modes.get("all", {})
    mode_layout = layout_modes.get(mode, {})
    all_positions = all_layout.get("read_to_y", {})
    mode_positions = mode_layout.get("read_to_y", {})
    if not isinstance(all_positions, Mapping) or not isinstance(mode_positions, Mapping):
        return None
    if read_name not in all_positions or read_name not in mode_positions:
        return None
    return float(mode_positions[read_name]) - float(all_positions[read_name])


def _add_scalar_layout_mode_columns(
    source_data: dict,
    layout_modes: Mapping[str, Mapping[str, object]],
    column_names: Sequence[str],
) -> None:
    """Add alternate scalar y columns to a read-indexed source dictionary."""
    layout_read_names = source_data.get("layout_read_name") or source_data.get("read_name")
    if not layout_read_names:
        return
    for mode in READ_FILTER_LAYOUT_MODES:
        for column_name in column_names:
            if column_name not in source_data:
                continue
            mode_values = []
            for row_index, base_value in enumerate(source_data[column_name]):
                delta = _layout_mode_delta(
                    str(layout_read_names[row_index]),
                    mode,
                    layout_modes,
                )
                mode_values.append(base_value if delta is None else base_value + delta)
            source_data[_layout_mode_column(column_name, mode)] = mode_values


def _add_line_layout_mode_columns(
    source_data: dict,
    layout_modes: Mapping[str, Mapping[str, object]],
    column_name: str,
) -> None:
    """Add alternate y polyline columns to a read-indexed source dictionary."""
    layout_read_names = source_data.get("layout_read_name") or source_data.get("read_name")
    if not layout_read_names or column_name not in source_data:
        return
    for mode in READ_FILTER_LAYOUT_MODES:
        mode_values = []
        for row_index, base_values in enumerate(source_data[column_name]):
            delta = _layout_mode_delta(
                str(layout_read_names[row_index]),
                mode,
                layout_modes,
            )
            if delta is None:
                mode_values.append(list(base_values))
            else:
                mode_values.append([value + delta for value in base_values])
        source_data[_layout_mode_column(column_name, mode)] = mode_values


def _add_clickable_layout_mode_columns(
    clickable_data: dict,
    layout_modes: Mapping[str, Mapping[str, object]],
) -> None:
    """Add alternate y columns to clickable-label data before source creation."""
    clickable_data["read_name"] = [
        str(info.get("read_name", "")) for info in clickable_data["customdata"]
    ]
    clickable_data["layout_read_name"] = [
        str(info.get("layout_read_name", info.get("read_name", "")))
        for info in clickable_data["customdata"]
    ]
    _add_scalar_layout_mode_columns(clickable_data, layout_modes, ("y",))


def _add_read_layout_mode_columns(
    arrow_data: dict,
    clickable_data: dict,
    variant_data: dict,
    connector_data: dict,
    same_region_connector_data: dict,
    layout_modes: Mapping[str, Mapping[str, object]],
) -> None:
    """Annotate read-derived Bokeh source dictionaries with compact y layouts."""
    _add_scalar_layout_mode_columns(arrow_data, layout_modes, ("y", "y0", "y1"))
    _add_clickable_layout_mode_columns(clickable_data, layout_modes)
    for data in variant_data.values():
        _add_scalar_layout_mode_columns(data, layout_modes, ("y",))
    _add_scalar_layout_mode_columns(connector_data, layout_modes, ("y",))
    _add_line_layout_mode_columns(same_region_connector_data, layout_modes, "ys")


def _layout_mode_y_bounds(
    layout_modes: Mapping[str, Mapping[str, object]],
) -> dict[str, list[float]]:
    """Return y-range bounds for each compact read layout mode."""
    bounds = {}
    for mode in READ_FILTER_LAYOUT_MODES:
        layout = layout_modes.get(mode, {})
        alignments_height = float(layout.get("alignments_height", 0.0))
        top_padding = max(alignments_height * 0.02, 0.5)
        bounds[mode] = [alignments_height, -top_padding]
    return bounds


def _group_end_y_by_layout_mode(
    haplotype: int,
    layout_modes: Mapping[str, Mapping[str, object]],
    fallback_y: float,
) -> dict[str, float]:
    """Return insertion-summary marker y positions for all layout modes."""
    y_by_mode = {}
    for mode in READ_FILTER_LAYOUT_MODES:
        group_boundaries = layout_modes.get(mode, {}).get("group_boundaries", {})
        if isinstance(group_boundaries, Mapping) and haplotype in group_boundaries:
            y_by_mode[mode] = float(group_boundaries[haplotype][1]) + 0.3
        else:
            y_by_mode[mode] = fallback_y
    return y_by_mode


def phase_block_boundaries_by_haplotype(
    segments_by_read: Mapping[str, Sequence[object]],
    coordinate_start: int,
    coordinate_end: int,
) -> dict[int, list[int]]:
    """Return visible phase-set start/end markers grouped by haplotype."""
    block_ends: dict[tuple[int, int], int] = {}
    for segments in segments_by_read.values():
        for segment in segments:
            haplotype_tag = getattr(segment, "haplotype_tag", None)
            phaseset_tag = getattr(segment, "phaseset_tag", None)
            if haplotype_tag is None or phaseset_tag is None:
                continue
            try:
                haplotype = int(haplotype_tag)
                phaseset_start = int(phaseset_tag)
            except (TypeError, ValueError):
                continue
            if haplotype <= 0:
                continue

            segment_start = int(getattr(segment, "pos", coordinate_end + 1))
            segment_end = int(getattr(segment, "end", coordinate_start - 1))
            if segment_start > coordinate_end or segment_end < coordinate_start:
                continue

            block_key = (haplotype, phaseset_start)
            block_ends[block_key] = max(block_ends.get(block_key, segment_end), segment_end)

    boundaries: dict[int, set[int]] = {}
    for (haplotype, phaseset_start), block_end in block_ends.items():
        for position in (phaseset_start, block_end):
            if coordinate_start <= position <= coordinate_end:
                boundaries.setdefault(haplotype, set()).add(position)

    return {haplotype: sorted(positions) for haplotype, positions in boundaries.items()}


def _phase_block_spans(
    group_boundaries: Mapping[int, tuple[float, float]],
) -> dict[int, tuple[float, float]]:
    """Return y spans for phase-block markers by haplotype."""
    ordered_groups = sorted(group_boundaries.items(), key=lambda item: item[1][0])
    block_spans = {}
    for group_index, (haplotype, (y_start, y_end)) in enumerate(ordered_groups):
        block_y_start = y_start
        block_y_end = y_end
        if group_index > 0:
            previous_y_end = ordered_groups[group_index - 1][1][1]
            block_y_start = (previous_y_end + y_start) / 2
        if group_index + 1 < len(ordered_groups):
            next_y_start = ordered_groups[group_index + 1][1][0]
            block_y_end = (y_end + next_y_start) / 2
        block_spans[haplotype] = (block_y_start, block_y_end)
    return block_spans


def _build_phase_boundaries_raw_data(
    phase_boundaries: Mapping[int, Sequence[int]],
    group_boundaries: Mapping[int, tuple[float, float]],
    coordinate_start: int,
    coordinate_end: int,
    layout_modes: Mapping[str, Mapping[str, object]] | None = None,
) -> dict | None:
    """Compute phase-block boundary source data without creating Bokeh models.

    Used by add_phase_block_boundaries_to_plot (rendering) and _serialize_region_for_swap.
    Returns None when there are no boundaries to draw.
    """
    xs: list = []
    ys: list = []
    colors: list = []
    row_haplotypes: list = []
    block_spans = _phase_block_spans(group_boundaries)
    for haplotype, positions in phase_boundaries.items():
        if haplotype not in block_spans:
            continue
        y_start, y_end = block_spans[haplotype]
        if y_start == y_end:
            continue
        for position in positions:
            if coordinate_start <= position <= coordinate_end:
                xs.append([position, position])
                ys.append([y_start, y_end])
                colors.append(PLOT_CONFIG["sample_label_color"])
                row_haplotypes.append(haplotype)
    if not xs:
        return None
    source_data: dict = {
        "xs": xs,
        "ys": ys,
        "color": colors,
        "read_layout_mode": ["all"] * len(xs),
    }
    for mode in READ_FILTER_LAYOUT_MODES:
        mode_boundaries: dict = {}
        if layout_modes:
            mode_boundaries = layout_modes.get(mode, {}).get("group_boundaries", {})
        mode_spans = _phase_block_spans(mode_boundaries) if mode_boundaries else block_spans
        source_data[_layout_mode_column("ys", mode)] = [
            [
                mode_spans.get(haplotype, block_spans[haplotype])[0],
                mode_spans.get(haplotype, block_spans[haplotype])[1],
            ]
            for haplotype in row_haplotypes
        ]
    add_read_filter_visibility_columns(source_data)
    return source_data


def add_phase_block_boundaries_to_plot(
    plot_figure,
    phase_boundaries: Mapping[int, Sequence[int]],
    group_boundaries: Mapping[int, tuple[float, float]],
    coordinate_start: int,
    coordinate_end: int,
    layout_modes: Mapping[str, Mapping[str, object]] | None = None,
):
    """Render phase-block boundary rules inside each haplotype lane."""
    source_data = _build_phase_boundaries_raw_data(
        phase_boundaries, group_boundaries, coordinate_start, coordinate_end, layout_modes
    )
    if source_data is None:
        return None
    source = ColumnDataSource(source_data)
    return plot_figure.multi_line(
        xs="xs",
        ys="ys",
        source=source,
        line_color="color",
        line_width=PLOT_CONFIG["phase_block_line_width"],
        line_alpha=PLOT_CONFIG["phase_block_line_alpha"],
        line_dash=PLOT_CONFIG["phase_block_line_dash"],
        level="underlay",
    )


def add_vcf_variants(vcf_figure, vcf_variants, coordinate_start, coordinate_end, sample_label=None):
    """Add VCF variants to a separate figure; make them clickable for modal."""
    if not vcf_variants:
        return
    variant_xs, variant_ys, variant_colors, variant_angles = [], [], [], []
    (
        variant_coordinates,
        variant_types,
        variant_alt_alleles,
        variant_alt_bases,
        variant_haplotypes,
    ) = ([], [], [], [], [])
    vcf_y = 0.5
    sample_label_str = sample_label or ""

    for variant in vcf_variants:
        if variant.variant_type == "DELETION":
            variant_center = variant.pos + (len(variant.ref) - 1) / 2.0
            variant_end = variant.pos + len(variant.ref) - 1
            coord_str = f"{variant.chrom}:{variant.pos:,}-{variant_end:,}"
        else:
            variant_center = variant.pos
            coord_str = f"{variant.chrom}:{variant.pos:,}"
        if variant_center < coordinate_start or variant_center > coordinate_end:
            continue
        variant_xs.append(variant_center)
        variant_ys.append(vcf_y)
        variant_colors.append(get_vcf_variant_color(variant))
        variant_angles.append(math.pi)
        variant_coordinates.append(coord_str)
        variant_types.append(variant.variant_type)
        variant_alt_alleles.append(variant.alt)
        variant_alt_bases.append(variant.alt_base if variant.alt_base else "")
        variant_haplotypes.append(", ".join(variant.haplotypes) if variant.haplotypes else "None")

    variant_source = None
    if variant_xs:
        variant_source = ColumnDataSource(
            data={
                "x": variant_xs,
                "y": variant_ys,
                "color": variant_colors,
                "angle": variant_angles,
                "coordinates": variant_coordinates,
                "variant_type": variant_types,
                "alt_allele": variant_alt_alleles,
                "alt_base": variant_alt_bases,
                "haplotypes": variant_haplotypes,
                "sample_label": [sample_label_str] * len(variant_xs),
            }
        )
        variant_renderer = vcf_figure.scatter(
            x="x",
            y="y",
            source=variant_source,
            marker="triangle",
            size=PLOT_CONFIG["vcf_triangle_size"],
            angle="angle",
            fill_color="color",
            fill_alpha=PLOT_CONFIG["vcf_triangle_fill_alpha"],
            line_color="color",
            line_alpha=PLOT_CONFIG["vcf_triangle_line_alpha"],
            line_width=PLOT_CONFIG["vcf_triangle_line_width"],
        )
        tap_tool = TapTool()
        tap_tool.renderers = [variant_renderer]
        vcf_figure.add_tools(tap_tool)
        vcf_figure.toolbar.active_tap = tap_tool
        variant_click_callback = get_vcf_variant_click_callback(variant_source)
        variant_source.selected.js_on_change("indices", variant_click_callback)
    return variant_source


def add_arrows_to_plot(plot_figure, arrow_data, allow_empty=False):
    """Add arrow segments and arrowheads to the plot using batched glyphs."""
    return plot_renderers.add_arrows_to_plot(plot_figure, arrow_data, allow_empty)


def add_read_connection_markers_to_plot(plot_figure, connector_data):
    """Render in-panel read-continuation endpoint markers."""
    if not connector_data["stub_x0"]:
        return None
    cfg = PLOT_CONFIG
    row_count = len(connector_data["stub_x0"])
    connector_data = dict(connector_data)
    connector_data.setdefault("read_filter_alpha", [1.0] * row_count)
    connector_data.setdefault("has_split_alignment", [False] * row_count)
    connector_data.setdefault("has_multiregion_connection", [False] * row_count)
    connector_data.setdefault("connector_line_alpha", [0.9] * row_count)
    connector_data.setdefault("connector_selected_alpha", [1.0] * row_count)
    connector_data.setdefault(
        "connector_nonselected_alpha",
        [cfg["connector_nonselection_alpha"]] * row_count,
    )
    add_read_filter_visibility_columns(connector_data)
    connector_source = ColumnDataSource(data=connector_data)
    marker_renderer = plot_figure.scatter(
        x="marker_x",
        y="y",
        source=connector_source,
        marker="circle",
        size=cfg["connector_marker_size"],
        fill_color=cfg["connector_marker_fill_color"],
        fill_alpha=0.0,
        line_color=cfg["connector_marker_line_color"],
        line_alpha="connector_line_alpha",
        line_width=2.0,
        selection_fill_color=cfg["connector_selection_color"],
        selection_fill_alpha=0.0,
        selection_line_color=cfg["connector_selection_color"],
        selection_line_alpha="connector_selected_alpha",
        nonselection_fill_alpha=cfg["connector_nonselection_alpha"],
        nonselection_line_alpha="connector_nonselected_alpha",
        name="orographer_read_connector_marker",
    )
    return connector_source, None, marker_renderer


def add_same_region_read_connection_lines_to_plot(plot_figure, connector_data):
    """Render same-region read-continuation lines as discrete selectable paths."""
    if not connector_data["xs"]:
        return None
    cfg = PLOT_CONFIG
    row_count = len(connector_data["xs"])
    connector_data = dict(connector_data)
    connector_data.setdefault("read_filter_alpha", [1.0] * row_count)
    connector_data.setdefault("has_split_alignment", [False] * row_count)
    connector_data.setdefault("has_multiregion_connection", [False] * row_count)
    connector_data.setdefault("connector_line_alpha", [cfg["connector_line_alpha"]] * row_count)
    connector_data.setdefault(
        "connector_selected_alpha",
        [cfg["connector_selection_line_alpha"]] * row_count,
    )
    connector_data.setdefault(
        "connector_nonselected_alpha",
        [cfg["connector_nonselection_alpha"]] * row_count,
    )
    add_read_filter_visibility_columns(connector_data)
    connector_source = ColumnDataSource(data=connector_data)
    renderer = plot_figure.multi_line(
        xs="xs",
        ys="ys",
        source=connector_source,
        line_color=cfg["connector_line_color"],
        line_width=cfg["connector_line_width"],
        line_alpha="connector_line_alpha",
        line_dash="dashed",
        selection_line_color=cfg["connector_selection_line_color"],
        selection_line_alpha="connector_selected_alpha",
        selection_line_width=cfg["connector_selection_line_width"],
        nonselection_line_alpha="connector_nonselected_alpha",
    )
    plot_figure.multi_line(
        xs="xs",
        ys="ys",
        source=connector_source,
        line_color=cfg["connector_selection_line_color"],
        line_width=cfg["connector_hit_line_width"],
        line_alpha=0.001,
        selection_line_alpha=0.001,
        selection_line_width=cfg["connector_hit_line_width"],
        nonselection_line_alpha=0.001,
        nonselection_line_width=cfg["connector_hit_line_width"],
        name="orographer_same_region_connector_hit_area",
    )
    return connector_source, renderer


def add_variants_to_plot(plot_figure, variant_data):
    """Add variant markers (mismatches, insertions, deletions); return LOD/1bp."""
    return plot_renderers.add_variants_to_plot(plot_figure, variant_data)


def setup_variant_lod_callback(
    plot_figure, variant_renderers, one_bp_renderers=None, hide_1bp_checkbox=None
):
    """Set up LOD callback to toggle marker vs text views.
    If one_bp_renderers and hide_1bp_checkbox are provided, keep 1bp hidden when active.
    """
    if not variant_renderers["marker"] and not variant_renderers["text"]:
        return
    args = {
        "marker_renderers": variant_renderers["marker"],
        "text_renderers": variant_renderers["text"],
        "x_range": plot_figure.x_range,
    }
    if one_bp_renderers is not None and hide_1bp_checkbox is not None:
        args["one_bp_renderers"] = one_bp_renderers
        args["hide_1bp_checkbox"] = hide_1bp_checkbox
    callback = CustomJS(
        args=args,
        code=load_javascript("variant_lod_callback.js"),
        module=True,
    )
    plot_figure.x_range.js_on_change("start", callback)
    plot_figure.x_range.js_on_change("end", callback)


def clickable_label_source_data(clickable_data):
    """Return ColumnDataSource data for clickable alignment labels."""
    return source_data.clickable_label_source_data(clickable_data)


def empty_clickable_label_source_data():
    """Return empty label source columns matching clickable_label_source_data()."""
    return source_data.empty_clickable_label_source_data()


def add_clickable_labels(
    plot_figure,
    tap_tool,
    clickable_data,
    arrow_source=None,
    arrow_renderer=None,
    arrow_tap_callback=None,
    reset_callback=None,
    allow_empty=False,
):
    """Add clickable number labels with modal; optional arrow tap/reset callbacks.
    Returns list of renderers for visibility toggling, or None if no labels.
    """
    return plot_renderers.add_clickable_labels(
        plot_figure,
        tap_tool,
        clickable_data,
        arrow_source=arrow_source,
        arrow_renderer=arrow_renderer,
        arrow_tap_callback=arrow_tap_callback,
        reset_callback=reset_callback,
        allow_empty=allow_empty,
    )


def add_alignment_number_modal_callbacks(alignment_label_sources) -> None:
    """Attach alignment-number modal callbacks with access to every label source."""
    for source in alignment_label_sources:
        source.selected.js_on_change(
            "indices",
            get_number_click_callback(source, alignment_label_sources),
        )


def add_vcf_track_to_region(
    plot_figure, vcf_variants, coordinate_start, coordinate_end, sample_label=None
) -> tuple[object | None, object | None]:
    """Add VCF variant track to a region plot if variants are provided.

    Returns (vcf_figure, vcf_source); both are None when vcf_variants is absent.
    """
    if not vcf_variants:
        return None, None
    vcf_figure = create_vcf_track_figure(plot_figure)
    vcf_source = add_vcf_variants(
        vcf_figure,
        vcf_variants,
        coordinate_start,
        coordinate_end,
        sample_label=sample_label,
    )
    return vcf_figure, vcf_source


def add_gene_track_to_region(
    plot_figure, gene_annotations, coordinate_start, coordinate_end
) -> tuple[object | None, object, dict]:
    """Add gene annotation track to a region plot if annotations are provided.

    Returns (gene_figure, plot_figure, gene_sources_dict) where gene_sources_dict
    maps "body"/"exon"/"intron"/"arrow"/"label" to ColumnDataSource or None.
    """
    if not gene_annotations:
        return None, plot_figure, {}
    gene_track_height = 4.0
    gene_figure = create_gene_track_figure(plot_figure, gene_track_height)
    actual_gene_track_height, gene_sources = add_gene_track(
        gene_figure, gene_annotations, 0, coordinate_start, coordinate_end
    )
    label_padding = 0.5
    if actual_gene_track_height > gene_track_height:
        gene_figure.y_range.start = actual_gene_track_height
    gene_figure.y_range.end = -label_padding
    return gene_figure, plot_figure, gene_sources


def _hide_xaxis_on_track_figures(region_state, gene_figure):
    """All track figures hide x-axis; sole axis is ``create_genomic_x_axis_strip``."""
    for fig in region_state.plot_figures:
        fig.xaxis.visible = False
        fig.xaxis.axis_label = None
    for fig in region_state.vcf_figures:
        if fig is not None:
            fig.xaxis.visible = False
            fig.xaxis.axis_label = None
    if gene_figure:
        gene_figure.xaxis.visible = False
        gene_figure.xaxis.axis_label = None


def add_multi_region_callbacks(
    all_plot_figures,
    all_arrow_sources,
    all_arrow_renderers,
    all_selectable_sources=None,
    all_region_reset_bounds=None,
):
    """Add cross-region highlighting callbacks and per-region range reset."""
    if len(all_arrow_sources) == 0:
        return
    selection_sources = all_selectable_sources or all_arrow_sources
    for plot_fig, arrow_source, _arrow_renderer in zip(
        all_plot_figures, all_arrow_sources, all_arrow_renderers, strict=False
    ):
        if len(selection_sources) == 1:
            single_region_callback = get_arrow_tap_callback(arrow_source)
            plot_fig.js_on_event("tap", single_region_callback)
        else:
            multi_region_callback = get_arrow_tap_callback_multi_region(
                arrow_source, selection_sources
            )
            plot_fig.js_on_event("tap", multi_region_callback)
    if all_region_reset_bounds:
        for (
            region_plot_figures,
            _region_sources,
            _region_renderers,
            reset_bounds,
        ) in all_region_reset_bounds:
            range_reset_callback = CustomJS(
                args={
                    "x_range": reset_bounds["x_range"],
                    "x_start": reset_bounds["x_start"],
                    "x_end": reset_bounds["x_end"],
                    "y_ranges": reset_bounds["y_ranges"],
                    "y_bounds": reset_bounds["y_bounds"],
                    "all_sources": reset_bounds["all_sources"],
                },
                code=load_javascript("range_reset_callback.js"),
            )
            for plot_fig in region_plot_figures:
                plot_fig.js_on_event("reset", range_reset_callback)
    else:
        shared_reset_callback = CustomJS(
            args={"all_sources": selection_sources},
            code=load_javascript("shared_reset_callback.js"),
        )
        for plot_fig in all_plot_figures:
            plot_fig.js_on_event("reset", shared_reset_callback)


def add_read_connection_overlay_callbacks(
    all_plot_figures,
    all_connector_sources,
    all_selectable_sources=None,
) -> None:
    """Redraw explicit connection overlays when selected endpoints or ranges move."""
    if not all_connector_sources:
        return

    overlay_callback = CustomJS(
        args={"all_sources": all_selectable_sources or []},
        code=load_javascript("read_connection_overlay_callback.js"),
    )
    seen_ranges = set()
    for plot_fig in all_plot_figures:
        for axis_range in (plot_fig.x_range, plot_fig.y_range):
            if axis_range is None:
                continue
            range_id = id(axis_range)
            if range_id in seen_ranges:
                continue
            seen_ranges.add(range_id)
            axis_range.js_on_change("start", overlay_callback)
            axis_range.js_on_change("end", overlay_callback)


def add_alignment_label_selection_callbacks(selection_sources, alignment_label_sources) -> None:
    """Fade alignment number labels for unselected reads while any read is selected."""
    if not selection_sources or not alignment_label_sources:
        return

    label_callback = CustomJS(
        args={
            "selection_sources": selection_sources,
            "alignment_label_sources": alignment_label_sources,
        },
        code=load_javascript("alignment_label_selection_callback.js"),
    )
    for source in selection_sources:
        source.selected.js_on_change("indices", label_callback)


def create_read_search_button(all_selectable_sources):
    """Create a compact button that opens read-name selection controls in the modal."""
    button_width = 90
    button_height = 22
    button_spacing = 4
    button_render_padding = 14
    compact_button_styles = {
        "font-size": "11px",
        "line-height": "1",
        "min-width": "0",
        "padding": "0",
        "text-align": "center",
        "width": f"{button_width}px",
    }
    read_button = Button(
        label="Select reads",
        button_type="default",
        width=button_width,
        height=button_height,
        margin=(0, 0, 0, 0),
        styles=compact_button_styles.copy(),
    )
    read_button.js_on_click(
        CustomJS(
            args={"all_sources": all_selectable_sources},
            code=load_javascript("read_search_modal_callback.js"),
        )
    )
    clear_button = Button(
        label="Clear selected",
        button_type="default",
        width=button_width,
        height=button_height,
        disabled=True,
        margin=(0, 0, 0, 0),
        styles=compact_button_styles.copy(),
    )
    clear_button.js_on_click(
        CustomJS(
            args={"all_sources": all_selectable_sources, "clear_button": clear_button},
            code=load_javascript("clear_read_selection_callback.js"),
        )
    )
    clear_state_callback = CustomJS(
        args={"all_sources": all_selectable_sources, "clear_button": clear_button},
        code=load_javascript("clear_read_selection_state_callback.js"),
    )
    for source in all_selectable_sources:
        source.selected.js_on_change("indices", clear_state_callback)
    return row(
        read_button,
        clear_button,
        width=(button_width + button_render_padding) * 2 + button_spacing,
        height=button_height,
        spacing=button_spacing,
        sizing_mode="fixed",
    )


def _stable_gene_color(gene_id: str) -> str:
    """Return a stable color for transcript glyphs from the gene of origin."""
    return isoseq_data.stable_gene_color(gene_id)


def _format_span(start: int, end: int) -> str:
    return isoseq_data.format_span(start, end)


def _intron_direction_positions(
    intron_start: int,
    intron_end: int,
    region_span: int,
) -> list[float]:
    """Return capped, evenly spaced x positions for transcript strand markers."""
    return isoseq_data.intron_direction_positions(intron_start, intron_end, region_span)


def _safe_chunk_token(value) -> str:
    """Return a filesystem-safe token for lazy-loaded chunk files."""
    return isoseq_data.safe_chunk_token(value)


def _natural_sort_key(value):
    """Return a case-insensitive sort key that treats digit runs numerically."""
    return isoseq_data.natural_sort_key(value)


def _empty_arrow_source_data():
    """Return empty read-glyph columns compatible with add_arrows_to_plot()."""
    return source_data.empty_arrow_source_data()


def _empty_isoseq_coverage_data():
    return source_data.empty_isoseq_coverage_data()


def _coverage_blocks_for_read(read_name, segments_by_read, coordinate_start, coordinate_end):
    yield from coverage_data.coverage_blocks_for_read(
        read_name,
        segments_by_read,
        coordinate_start,
        coordinate_end,
    )


def _isoseq_coverage_block_cache(segments_by_read, coordinate_start, coordinate_end):
    """Cache per-read coverage blocks once for transcript coverage assembly."""
    return coverage_data.coverage_block_cache(segments_by_read, coordinate_start, coordinate_end)


def _isoseq_coverage_for_cached_reads(
    read_names,
    coverage_blocks,
    coordinate_start,
    coordinate_end,
):
    return coverage_data.coverage_for_cached_reads(
        read_names,
        coverage_blocks,
        coordinate_start,
        coordinate_end,
    )


def _isoseq_coverage_for_blocks(blocks, coordinate_start, coordinate_end):
    return coverage_data.coverage_for_blocks(blocks, coordinate_start, coordinate_end)


def _isoseq_read_page_url_template(chunk_url_prefix, chunk_key):
    """Return a browser-loadable static manifest URL for IsoSeq read pages."""
    return isoseq_store.read_manifest_url(chunk_url_prefix, chunk_key)


def _isoseq_coverage_url(chunk_url_prefix, chunk_key):
    """Return a browser-loadable static JSON URL for one transcript coverage payload."""
    return isoseq_store.coverage_url(chunk_url_prefix, chunk_key)


def _write_isoseq_read_page(
    chunk_dir,
    chunk_key,
    page,
    payload,
):
    """Write one static IsoSeq read page payload."""
    if not chunk_dir:
        return
    os.makedirs(chunk_dir, exist_ok=True)
    page_path = os.path.join(chunk_dir, f"{chunk_key}_page{page}.json.gz")
    compact_payload = _compact_isoseq_read_page_payload(payload)
    write_json_gzip(page_path, compact_payload)


def _write_isoseq_coverage_page(
    chunk_dir,
    chunk_key,
    coverage_data,
):
    """Write one static IsoSeq transcript coverage payload."""
    isoseq_store.write_coverage_page(chunk_dir, chunk_key, coverage_data)


def _compact_isoseq_read_page_payload(payload):
    """Return a lean payload that browser code expands into Bokeh source columns."""
    arrow_data = payload["arrow_data"]
    label_data = payload["label_data"]
    first_coordinate = next(iter(label_data.get("coordinates", [])), "")
    chrom = first_coordinate.split(":", 1)[0] if ":" in first_coordinate else ""
    read_index_by_name = {}
    reads = {
        "name": [],
        "gene_id": [],
        "gene_name": [],
        "transcript_id": [],
        "annotation_label": [],
        "group_id": [],
        "strand": [],
        "haplotype": [],
        "sample_label": [],
        "all_alignment_numbers": [],
        "all_alignment_coordinates": [],
    }
    block_read_indices = []
    for row_index, read_name in enumerate(arrow_data["read_name"]):
        if read_name not in read_index_by_name:
            read_index_by_name[read_name] = len(reads["name"])
            reads["name"].append(read_name)
            reads["gene_id"].append(arrow_data["gene_id"][row_index])
            reads["gene_name"].append(arrow_data["gene_name"][row_index])
            reads["transcript_id"].append(arrow_data["transcript_id"][row_index])
            reads["annotation_label"].append(
                label_data.get("annotation_label", [""])[row_index]
                if label_data.get("annotation_label") else ""
            )
            reads["group_id"].append(arrow_data["group_id"][row_index])
            reads["strand"].append(label_data["strand"][row_index])
            reads["haplotype"].append(label_data["haplotype"][row_index])
            reads["sample_label"].append(label_data["sample_label"][row_index])
            reads["all_alignment_numbers"].append(label_data["all_alignment_numbers"][row_index])
            reads["all_alignment_coordinates"].append(
                label_data["all_alignment_coordinates"][row_index]
            )
        block_read_indices.append(read_index_by_name[read_name])

    blocks = {
        "read_index": block_read_indices,
        "x0": arrow_data["x0"],
        "x1": arrow_data["x1"],
        "y": arrow_data["y"],
        "y0": arrow_data["y0"],
        "y1": arrow_data["y1"],
        "color": arrow_data["color"],
        "segment_id": arrow_data["segment_id"],
        "alignment_order": arrow_data["alignment_order"],
        "fwd_read_start": arrow_data["fwd_read_start"],
        "fwd_read_end": arrow_data["fwd_read_end"],
    }
    return {
        "schema": "isoseq_compact_v1",
        "chrom": chrom,
        "reads": reads,
        "blocks": blocks,
        "page": payload["page"],
        "page_size": payload["page_size"],
        "loaded_read_count": payload["loaded_read_count"],
        "assigned_read_count": payload["assigned_read_count"],
    }


def _isoseq_block_number_by_index(segment_blocks, is_fwd_strand):
    return isoseq_data.block_number_by_index(segment_blocks, is_fwd_strand)


def _isoseq_read_shard_record(read_name, segments, read_metadata):
    return isoseq_data.read_shard_record(read_name, segments, read_metadata)


def _write_isoseq_read_store(
    chunk_dir,
    chunk_url_prefix,
    manifest_key,
    groups,
    read_records,
    read_id_by_name,
    manifest_context,
):
    isoseq_store.write_read_store(
        chunk_dir,
        manifest_key,
        groups,
        read_records,
        read_id_by_name,
        manifest_context,
    )


def _alignment_summaries_for_segments(row_index, segments_by_read):
    summaries = {}
    for display_read_name, segments in segments_by_read.items():
        for segment in segments:
            read_name = getattr(segment, "readname", display_read_name)
            summaries.setdefault((row_index, read_name), []).append(
                {
                    "alignment_number": getattr(segment, "alignment_order", 0),
                    "coordinates": format_alignment_coordinates(
                        segment.chrom,
                        segment.pos,
                        segment.end,
                    ),
                    "region_idx": 0,
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


def _normalise_isoseq_arrow_source_data(arrow_data):
    """Return read glyph data with every browser source column populated."""
    return isoseq_data.normalise_arrow_source_data(arrow_data)


def _render_isoseq_read_page_payload(
    group,
    page_read_names,
    segments_by_read,
    read_metadata,
    coordinate_start,
    coordinate_end,
    region_idx,
    row_index,
    sample_label,
    plot_dom_class,
    plot_model_id,
    selected_read_y_start,
    page,
):
    """Render one portable IsoSeq read page without any runtime BAM dependency."""
    read_height = 0.24
    read_to_y = {
        read_name: selected_read_y_start + index * read_height + read_height / 2
        for index, read_name in enumerate(page_read_names)
    }
    read_heights = dict.fromkeys(page_read_names, read_height)
    page_segments_by_read = {
        read_name: segments_by_read[read_name]
        for read_name in page_read_names
        if read_name in segments_by_read
    }
    alignment_summaries = _alignment_summaries_for_segments(row_index, page_segments_by_read)
    arrow_data, clickable_data, _variants, _connectors, _same_region = process_segments(
        page_segments_by_read,
        page_read_names,
        read_to_y,
        read_heights,
        coordinate_start,
        coordinate_end,
        ISOSEQ_REGION_TYPE,
        sample_label=sample_label or "",
        region_idx=region_idx,
        row_index=row_index,
        connection_index={},
        alignment_summaries=alignment_summaries,
        plot_dom_class=plot_dom_class,
        plot_model_id=plot_model_id,
    )
    page_read_metadata = {
        read_name: read_metadata.get(read_name, {}) for read_name in page_read_names
    }
    _add_isoseq_metadata_to_sources(arrow_data, clickable_data, page_read_metadata)
    return {
        "arrow_data": _normalise_isoseq_arrow_source_data(arrow_data),
        "label_data": clickable_label_source_data(clickable_data),
        "page": page,
        "page_size": ISOSEQ_READ_PAGE_SIZE,
        "loaded_read_count": len(page_read_names),
        "assigned_read_count": int(group.get("assigned_read_count", len(page_read_names))),
    }


def _prepare_isoseq_lazy_chunks(
    groups,
    segments_by_read,
    layout,
    coordinate_start,
    coordinate_end,
    region_idx,
    row_index,
    sample_label,
    plot_dom_class,
    plot_model_id,
    annotation_id="primary",
    chunk_dir=None,
    chunk_url_prefix=None,
    read_page_size=ISOSEQ_READ_PAGE_SIZE,
):
    """Write static sharded read data for browser-side IsoSeq lazy loading."""
    return isoseq_store.prepare_lazy_chunks(
        groups,
        segments_by_read,
        layout,
        coordinate_start,
        coordinate_end,
        region_idx,
        row_index,
        sample_label,
        annotation_id=annotation_id,
        chunk_dir=chunk_dir,
        chunk_url_prefix=chunk_url_prefix,
        read_page_size=read_page_size,
    )


def _build_isoseq_layout(groups, segments_by_read):
    """Return row positions for IsoSeq transcript/read groups."""
    return isoseq_data.build_layout(groups, segments_by_read)


def _add_isoseq_metadata_to_sources(arrow_data, clickable_data, read_metadata):
    """Attach transcript/gene metadata to read glyph and label sources."""
    isoseq_data.add_metadata_to_sources(arrow_data, clickable_data, read_metadata)


def _jsonable(data: dict) -> dict:
    """Convert ColumnDataSource data dict values to JSON-serialisable Python types."""
    result: dict = {}
    for key, values in data.items():
        if hasattr(values, "tolist"):
            result[key] = values.tolist()
        elif isinstance(values, list):
            cleaned = []
            for v in values:
                if hasattr(v, "tolist"):
                    cleaned.append(v.tolist())
                elif isinstance(v, list):
                    cleaned.append([x.item() if hasattr(x, "item") else x for x in v])
                else:
                    cleaned.append(v.item() if hasattr(v, "item") else v)
            result[key] = cleaned
        else:
            result[key] = values.item() if hasattr(values, "item") else values
    return result


def _vcf_data_for_region(
    vcf_variants: list | None,
    coordinate_start: int,
    coordinate_end: int,
    sample_label: str | None,
) -> dict | None:
    """Build the raw VCF data dict (same columns as add_vcf_variants) without Bokeh."""
    if not vcf_variants:
        return None
    xs, ys, colors, angles = [], [], [], []
    coords_col, types_col, alts_col, alt_bases_col, hps_col, labels_col = [], [], [], [], [], []
    sample_label_str = sample_label or ""
    for variant in vcf_variants:
        if variant.variant_type == "DELETION":
            center = variant.pos + (len(variant.ref) - 1) / 2.0
            end = variant.pos + len(variant.ref) - 1
            coord_str = f"{variant.chrom}:{variant.pos:,}-{end:,}"
        else:
            center = variant.pos
            coord_str = f"{variant.chrom}:{variant.pos:,}"
        if center < coordinate_start or center > coordinate_end:
            continue
        xs.append(center)
        ys.append(0.5)
        colors.append(get_vcf_variant_color(variant))
        angles.append(math.pi)
        coords_col.append(coord_str)
        types_col.append(variant.variant_type)
        alts_col.append(variant.alt)
        alt_bases_col.append(variant.alt_base if variant.alt_base else "")
        hps_col.append(", ".join(variant.haplotypes) if variant.haplotypes else "None")
        labels_col.append(sample_label_str)
    if not xs:
        return None
    return {
        "x": xs,
        "y": ys,
        "color": colors,
        "angle": angles,
        "coordinates": coords_col,
        "variant_type": types_col,
        "alt_allele": alts_col,
        "alt_base": alt_bases_col,
        "haplotypes": hps_col,
        "sample_label": labels_col,
    }


def _annotate_connector_source_data(connector_data: dict) -> dict:
    """Add default alpha/visibility columns matching add_read_connection_markers_to_plot."""
    if not connector_data.get("stub_x0"):
        return connector_data
    cfg = PLOT_CONFIG
    row_count = len(connector_data["stub_x0"])
    result = dict(connector_data)
    result.setdefault("read_filter_alpha", [1.0] * row_count)
    result.setdefault("has_split_alignment", [False] * row_count)
    result.setdefault("has_multiregion_connection", [False] * row_count)
    result.setdefault("connector_line_alpha", [0.9] * row_count)
    result.setdefault("connector_selected_alpha", [1.0] * row_count)
    result.setdefault(
        "connector_nonselected_alpha", [cfg["connector_nonselection_alpha"]] * row_count
    )
    add_read_filter_visibility_columns(result)
    return result


def _annotate_arc_connector_source_data(arc_data: dict) -> dict:
    """Add default alpha/visibility columns matching add_same_region_read_connection_lines."""
    if not arc_data.get("xs"):
        return arc_data
    cfg = PLOT_CONFIG
    row_count = len(arc_data["xs"])
    result = dict(arc_data)
    result.setdefault("read_filter_alpha", [1.0] * row_count)
    result.setdefault("has_split_alignment", [False] * row_count)
    result.setdefault("has_multiregion_connection", [False] * row_count)
    result.setdefault("connector_line_alpha", [cfg["connector_line_alpha"]] * row_count)
    result.setdefault(
        "connector_selected_alpha", [cfg["connector_selection_line_alpha"]] * row_count
    )
    result.setdefault(
        "connector_nonselected_alpha", [cfg["connector_nonselection_alpha"]] * row_count
    )
    add_read_filter_visibility_columns(result)
    return result


def _serialize_region_for_swap(
    region_data: dict,
    region_idx: int,
    connection_index: dict,
    alignment_summaries: dict,
) -> dict:
    """Return a compact JSON-safe dict for one region's swap data.

    Only complex-SV (non-IsoSeq) rows are serialised; IsoSeq rows produce None
    placeholders so indices stay aligned with slot source lists.
    """
    region = region_data["region"]
    coordinate_start = region.start
    coordinate_end = region.end
    bam_rows_data = []

    for row_index, bam_row in enumerate(region_data.get("bam_rows", [])):
        if bam_row.get("region_type") == ISOSEQ_REGION_TYPE:
            bam_rows_data.append(None)
            continue

        segments_by_read = bam_row.get("segments_by_read") or {}
        if not segments_by_read:
            bam_rows_data.append(None)
            continue

        region_type = bam_row.get("region_type", COMPLEX_SV_REGION_TYPE)
        sample_label = bam_row.get("sample_label")
        vcf_variants = bam_row.get("vcf_variants")
        coverage_tracks = bam_row.get("coverage_tracks")

        expected_haplotypes = _get_expected_haplotypes(bam_row)
        read_names, haplotype_groups, haplotype_order = sort_read_names(
            segments_by_read, expected_haplotypes
        )
        min_read_heights = connector_offset_min_heights_by_read(
            segments_by_read, region, region_idx, row_index, connection_index
        )
        read_to_y, _read_to_y_bottom, read_heights, alignments_height, _group_boundaries = (
            calculate_read_positions(
                read_names,
                segments_by_read,
                haplotype_groups,
                haplotype_order,
                min_read_heights_by_read=min_read_heights,
            )
        )
        read_filter_flags_by_read = _read_filter_flags_by_layout_read(
            segments_by_read, row_index, alignment_summaries
        )
        layout_modes = _build_read_filter_layout_modes(
            ReadFilterLayoutRequest(
                read_names=read_names,
                segments_by_read=segments_by_read,
                haplotype_groups=haplotype_groups,
                haplotype_order=haplotype_order,
                read_filter_flags_by_read=read_filter_flags_by_read,
                min_read_heights_by_read=min_read_heights,
            )
        )

        arrow_data, clickable_data, variant_data, connector_data, same_region_connector_data = (
            process_segments(
                segments_by_read,
                read_names,
                read_to_y,
                read_heights,
                coordinate_start,
                coordinate_end,
                region_type,
                sample_label=sample_label or "",
                region_idx=region_idx,
                row_index=row_index,
                connection_index=connection_index,
                alignment_summaries=alignment_summaries,
                # DOM class / model ID filled by JS after swap
                plot_dom_class="",
                plot_model_id="",
            )
        )
        _add_read_layout_mode_columns(
            arrow_data, clickable_data, variant_data, connector_data, same_region_connector_data,
            layout_modes,
        )
        add_read_filter_visibility_columns(arrow_data)
        # Populate chevron columns so the MultiLine renderer has valid 2-d data on swap.
        # The chevron_callback refines them to the correct pixel length when x_range updates.
        source_data.add_read_chevron_data(arrow_data)
        connector_data = _annotate_connector_source_data(connector_data)
        same_region_connector_data = _annotate_arc_connector_source_data(
            same_region_connector_data
        )

        color_haplotypes = list(haplotype_order)
        if coverage_tracks:
            for hp in coverage_tracks:
                if hp != -1 and hp not in color_haplotypes:
                    color_haplotypes.append(hp)
        haplotype_color_map = build_haplotype_color_map(color_haplotypes)
        cov_raw = build_coverage_source_data(coverage_tracks or {}, haplotype_color_map)
        vcf_raw = _vcf_data_for_region(vcf_variants, coordinate_start, coordinate_end, sample_label)

        def _hp_color_fn(hp: object, _cmap: dict = haplotype_color_map) -> str:
            return _cmap.get(hp, get_haplotype_color(hp))

        sep_data, lbl_data = _build_hp_label_raw_data(
            _group_boundaries,
            haplotype_order,
            coordinate_start,
            coordinate_end,
            haplotype_color_fn=_hp_color_fn,
            layout_modes=layout_modes,
        )
        separator_raw = _build_separator_raw_data(
            read_names,
            _read_to_y_bottom,
            read_heights,
            coordinate_start,
            coordinate_end,
            layout_modes,
            read_filter_flags_by_read,
        )
        phase_bnd_data = _build_phase_boundaries_raw_data(
            phase_block_boundaries_by_haplotype(segments_by_read, coordinate_start, coordinate_end),
            _group_boundaries,
            coordinate_start,
            coordinate_end,
            layout_modes,
        )
        variant_srcs_data = [
            _jsonable(d) for d in plot_renderers._build_variant_sources_data_list(variant_data)
        ]
        insertion_raw_data: list = []
        insertion_summary = bam_row.get("insertion_summary", {})
        for _ins_hp in sorted(insertion_summary.keys()):
            _ins_sites = insertion_summary[_ins_hp]
            if not _ins_sites or _ins_hp not in _group_boundaries:
                continue
            _ins_y_pos = _group_boundaries[_ins_hp][1] + 0.3
            _ins_y_by_mode = _group_end_y_by_layout_mode(_ins_hp, layout_modes, _ins_y_pos)
            insertion_raw_data.append(
                _jsonable(
                    _insertion_raw_source_data(_ins_sites, _ins_hp, _ins_y_pos, _ins_y_by_mode)
                )
            )
        top_padding = max(alignments_height * 0.02, 0.5)
        bam_rows_data.append({
            "arrow_data": _jsonable(arrow_data),
            "connector_data": _jsonable(connector_data),
            "same_region_connector_data": _jsonable(same_region_connector_data),
            "cov_total": _jsonable(cov_raw["total"]) if cov_raw.get("total") else None,
            "cov_hp": _jsonable(cov_raw["hp"]) if cov_raw.get("hp") else None,
            "cov_y_max": (
                float(max(cov_raw["total"]["y"]) * 1.05) if cov_raw.get("total") else None
            ),
            "vcf": _jsonable(vcf_raw) if vcf_raw else None,
            "y_bounds": [alignments_height, -top_padding],
            "hp_sep_data": _jsonable(sep_data) if sep_data else None,
            "hp_lbl_data": _jsonable(lbl_data) if lbl_data else None,
            "alignment_label_data": _jsonable(clickable_label_source_data(clickable_data)),
            "separator_data": _jsonable(separator_raw),
            "phase_boundary_data": _jsonable(phase_bnd_data) if phase_bnd_data else None,
            "variant_sources_data": variant_srcs_data if variant_srcs_data else None,
            "insertion_raw_data": insertion_raw_data if insertion_raw_data else None,
        })

    chromosome = region.chromosome
    gene_annotations = region_data.get("gene_annotations") or []
    gene_track_data = None
    if gene_annotations:
        _gene_h, body_d, exon_d, intron_d, arrow_d, label_d = _build_gene_track_raw_data(
            gene_annotations, 0, coordinate_start, coordinate_end
        )
        gene_track_data = {
            "body": _jsonable(body_d) if body_d.get("x0") else None,
            "exon": _jsonable(exon_d) if exon_d.get("left") else None,
            "intron": _jsonable(intron_d) if intron_d.get("xs") else None,
            "arrow": _jsonable(arrow_d) if arrow_d.get("x") else None,
            "label": _jsonable(label_d) if label_d.get("x") else None,
            "y_range_start": float(max(_gene_h, 4.0)),
        }
    repeat_density = region_data.get("repeat_density")
    repeat_density_data = None
    if repeat_density is not None:
        _dens_raw = _build_repeat_density_source_data(
            repeat_density, coordinate_start, coordinate_end
        )
        repeat_density_data = {
            "x": _dens_raw["x"],
            "top": _dens_raw["top"],
            "width": _dens_raw["width"],
            "y_max": float(max(repeat_density.max(), 1.0)),
        }
    return {
        "bam_rows": bam_rows_data,
        "x_start": coordinate_start,
        "x_end": coordinate_end,
        "chrom": chromosome,
        "orig_size_label": _format_original_region_size_label(coordinate_start, coordinate_end),
        "orig_coord_label": f"{chromosome}:{coordinate_start:,}-{coordinate_end:,}",
        "gene_track": gene_track_data,
        "repeat_density": repeat_density_data,
    }


def _get_expected_haplotypes(bam_row):
    """Return haplotypes that should get reserved slots for this BAM row."""
    coverage_tracks = bam_row.get("coverage_tracks")
    region_type = bam_row.get("region_type")
    if coverage_tracks and region_type == COMPLEX_SV_REGION_TYPE:
        return {hp for hp in coverage_tracks if hp != -1}
    return None


def _build_bam_row_track(
    bam_row,
    region,
    coordinate_start,
    coordinate_end,
    region_state,
    region_idx,
    row_index,
    connection_index=None,
    alignment_summaries=None,
    isoseq_chunk_dir=None,
    isoseq_chunk_url_prefix=None,
):
    """Add one BAM row to the region. Mutates region_state. Skips if no segments."""
    if bam_row.get("region_type") == ISOSEQ_REGION_TYPE:
        _build_isoseq_bam_row_track(
            bam_row,
            region,
            coordinate_start,
            coordinate_end,
            region_state,
            region_idx,
            row_index,
            alignment_summaries=alignment_summaries,
            add_vcf_track_to_region_fn=add_vcf_track_to_region,
            isoseq_chunk_dir=isoseq_chunk_dir,
            isoseq_chunk_url_prefix=isoseq_chunk_url_prefix,
        )
        return

    segments_by_read = bam_row["segments_by_read"]
    if not segments_by_read:
        logger.debug(f"No segments to plot for region {region_idx + 1}, row {row_index + 1}.")
        return

    vcf_variants = bam_row.get("vcf_variants")
    coverage_tracks = bam_row.get("coverage_tracks")
    region_type = bam_row["region_type"]
    sample_label = bam_row.get("sample_label")

    expected_haplotypes = _get_expected_haplotypes(bam_row)
    read_names, haplotype_groups, haplotype_order = sort_read_names(
        segments_by_read, expected_haplotypes
    )
    draw_read_names = read_names
    min_read_heights = connector_offset_min_heights_by_read(
        segments_by_read,
        region,
        region_idx,
        row_index,
        connection_index,
    )
    read_to_y, read_to_y_bottom, read_heights, alignments_height, group_boundaries = (
        calculate_read_positions(
            read_names,
            segments_by_read,
            haplotype_groups,
            haplotype_order,
            min_read_heights_by_read=min_read_heights,
        )
    )
    read_filter_flags_by_read = _read_filter_flags_by_layout_read(
        segments_by_read,
        row_index,
        alignment_summaries,
    )
    layout_modes = _build_read_filter_layout_modes(
        ReadFilterLayoutRequest(
            read_names=read_names,
            segments_by_read=segments_by_read,
            haplotype_groups=haplotype_groups,
            haplotype_order=haplotype_order,
            read_filter_flags_by_read=read_filter_flags_by_read,
            min_read_heights_by_read=min_read_heights,
        )
    )

    if region_state.shared_x_range is None:
        plot_figure, tap_tool = create_bokeh_figure(
            coordinate_start, coordinate_end, alignments_height
        )
        region_state.shared_x_range = plot_figure.x_range
        region_state.first_plot_figure = plot_figure
        region_state.region_type = region_type
    else:
        plot_figure, tap_tool = create_bokeh_figure_shared_x(
            region_state.shared_x_range, alignments_height
        )
    plot_dom_class = f"orographer-alignment-plot-r{region_idx}-row{row_index}"
    plot_figure.css_classes = [*plot_figure.css_classes, plot_dom_class]
    color_haplotypes = list(haplotype_order)
    if coverage_tracks:
        for hp in coverage_tracks:
            if hp != -1 and hp not in color_haplotypes:
                color_haplotypes.append(hp)
    haplotype_color_map = build_haplotype_color_map(color_haplotypes)

    def haplotype_color_fn(hp):
        return haplotype_color_map.get(hp, get_haplotype_color(hp))

    if sample_label:
        region_state.row_components.append(
            Div(
                text=sample_label,
                sizing_mode="stretch_width",
                styles={
                    "font-size": PLOT_CONFIG["sample_label_font_size"],
                    "font-family": "Arial, sans-serif",
                    "font-weight": "bold",
                    "padding-bottom": PLOT_CONFIG["sample_label_padding_bottom"],
                    "text-align": "center",
                    "width": "100%",
                    "box-sizing": "border-box",
                    "color": PLOT_CONFIG["sample_label_color"],
                },
            )
        )

    if coverage_tracks:
        cov_figure = create_coverage_figure(plot_figure)
        cov_sources = add_coverage_to_figure(cov_figure, coverage_tracks, haplotype_color_map)
        region_state.row_components.append(cov_figure)
        region_state.cursor_guide_figures.append(cov_figure)
        region_state.coverage_total_sources.append(cov_sources.get("total"))
        region_state.coverage_hp_sources.append(cov_sources.get("hp"))
        region_state.coverage_y_ranges_per_row.append(cov_figure.y_range)
    else:
        region_state.coverage_total_sources.append(None)
        region_state.coverage_hp_sources.append(None)
        region_state.coverage_y_ranges_per_row.append(None)

    vcf_figure, vcf_source = add_vcf_track_to_region(
        plot_figure,
        vcf_variants,
        coordinate_start,
        coordinate_end,
        sample_label=sample_label,
    )
    region_state.vcf_sources.append(vcf_source)
    if vcf_figure:
        region_state.row_components.append(vcf_figure)
        region_state.cursor_guide_figures.append(vcf_figure)
        region_state.vcf_figures.append(vcf_figure)
    else:
        region_state.vcf_figures.append(None)

    separator_source = add_separator_lines(
        plot_figure,
        read_names,
        read_to_y_bottom,
        read_heights,
        coordinate_start,
        coordinate_end,
        layout_modes=layout_modes,
        read_filter_flags_by_read=read_filter_flags_by_read,
    )
    region_state.read_filter_sources.append(separator_source)
    if separator_source is not None:
        region_state.clearable_sources.append(separator_source)
    region_state.separator_source_per_row.append(separator_source)
    phase_boundaries = phase_block_boundaries_by_haplotype(
        segments_by_read,
        coordinate_start,
        coordinate_end,
    )
    phase_set_marker_renderer = add_phase_block_boundaries_to_plot(
        plot_figure,
        phase_boundaries,
        group_boundaries,
        coordinate_start,
        coordinate_end,
        layout_modes=layout_modes,
    )
    if phase_set_marker_renderer is not None:
        region_state.phase_set_marker_renderers.append(phase_set_marker_renderer)
        region_state.read_filter_sources.append(phase_set_marker_renderer.data_source)
        region_state.clearable_sources.append(phase_set_marker_renderer.data_source)
    region_state.phase_set_source_per_row.append(
        phase_set_marker_renderer.data_source if phase_set_marker_renderer is not None else None
    )
    (
        arrow_data,
        clickable_data,
        variant_data,
        connector_data,
        same_region_connector_data,
    ) = process_segments(
        segments_by_read,
        draw_read_names,
        read_to_y,
        read_heights,
        coordinate_start,
        coordinate_end,
        region_type,
        sample_label=sample_label,
        region_idx=region_idx,
        row_index=row_index,
        connection_index=connection_index,
        alignment_summaries=alignment_summaries,
        plot_dom_class=plot_dom_class,
        plot_model_id=plot_figure.id,
    )
    _add_read_layout_mode_columns(
        arrow_data,
        clickable_data,
        variant_data,
        connector_data,
        same_region_connector_data,
        layout_modes,
    )
    arrow_result = add_arrows_to_plot(plot_figure, arrow_data)
    arrow_source = arrow_result[0] if arrow_result else None
    arrow_renderer = arrow_result[1] if arrow_result else None
    connector_result = add_read_connection_markers_to_plot(plot_figure, connector_data)
    connector_source = connector_result[0] if connector_result else None
    same_region_connector_result = add_same_region_read_connection_lines_to_plot(
        plot_figure,
        same_region_connector_data,
    )
    same_region_connector_source = (
        same_region_connector_result[0] if same_region_connector_result else None
    )
    if arrow_source:
        region_state.arrow_sources.append(arrow_source)
        region_state.read_filter_sources.append(arrow_source)
    if arrow_renderer:
        region_state.arrow_renderers.append(arrow_renderer)
    region_state.marker_connector_sources.append(connector_source)
    if connector_source:
        region_state.connector_sources.append(connector_source)
        region_state.read_filter_sources.append(connector_source)
    region_state.arc_connector_sources.append(same_region_connector_source)
    if same_region_connector_source:
        region_state.connector_sources.append(same_region_connector_source)
        region_state.read_filter_sources.append(same_region_connector_source)
    top_padding = max(alignments_height * 0.02, 0.5)
    region_state.y_bounds.append((plot_figure.y_range, alignments_height, -top_padding))
    region_state.read_filter_y_bounds.append(
        (plot_figure.y_range, _layout_mode_y_bounds(layout_modes))
    )
    region_state.slot_dom_classes.append(plot_dom_class)
    region_state.plot_model_ids.append(plot_figure.id)
    # Vertical scroll when zoomed is handled by WheelPanTool(dimensions="height") on the figure
    variant_renderers = add_variants_to_plot(plot_figure, variant_data)
    variant_sources = variant_renderers.get("sources", [])
    region_state.read_filter_sources.extend(variant_sources)
    region_state.clearable_sources.extend(variant_sources)
    region_state.variant_sources_per_row.append(list(variant_sources))
    row_one_bp = variant_renderers.get("one_bp", [])
    region_state.one_bp_renderers.extend(row_one_bp)
    region_state.one_bp_markers.extend(variant_renderers.get("one_bp_markers", []))
    region_state.one_bp_texts.extend(variant_renderers.get("one_bp_texts", []))
    region_state.one_bp_segments.extend(variant_renderers.get("one_bp_segments", []))
    region_state.variant_renderers.append(variant_renderers)
    region_state.one_bp_by_row.append(row_one_bp)
    label_result = add_clickable_labels(
        plot_figure,
        tap_tool,
        clickable_data,
        arrow_source,
        arrow_renderer,
        arrow_tap_callback=CustomJS(code=""),
    )
    if label_result:
        label_source, label_renderers = label_result
        region_state.alignment_label_sources.append(label_source)
        region_state.alignment_label_renderers.extend(label_renderers)
        region_state.read_filter_sources.append(label_source)
        region_state.clearable_sources.append(label_source)
        region_state.alignment_label_source_per_row.append(label_source)
    else:
        region_state.alignment_label_source_per_row.append(None)
    hp_label_sources_dict = add_haplotype_labels(
        plot_figure,
        group_boundaries,
        haplotype_order,
        coordinate_start,
        coordinate_end,
        haplotype_color_fn=haplotype_color_fn,
        layout_modes=layout_modes,
    )
    hp_label_sources = [v for v in hp_label_sources_dict.values() if v is not None]
    region_state.read_filter_sources.extend(hp_label_sources)
    region_state.clearable_sources.extend(hp_label_sources)
    region_state.hp_label_sources_per_row.append(hp_label_sources_dict)
    region_state.row_components.append(plot_figure)
    region_state.plot_figures.append(plot_figure)
    region_state.cursor_guide_figures.append(plot_figure)

    insertion_summary = bam_row.get("insertion_summary", {})
    row_insertion_raw_sources: list = []
    for hp in sorted(insertion_summary.keys()):
        sites = insertion_summary[hp]
        if not sites or hp not in group_boundaries:
            continue
        y_pos = group_boundaries[hp][1] + 0.3
        ins_renderer, insertion_sources = add_insertion_markers_to_figure(
            plot_figure,
            sites,
            hp,
            y_pos,
            coordinate_start,
            coordinate_end,
            y_by_layout_mode=_group_end_y_by_layout_mode(hp, layout_modes, y_pos),
            color=haplotype_color_fn(hp),
        )
        region_state.read_filter_sources.extend(insertion_sources)
        region_state.clearable_sources.extend(insertion_sources)
        row_insertion_raw_sources.append(insertion_sources[0])
        current_renderers = tap_tool.renderers
        if isinstance(current_renderers, list):
            tap_tool.renderers = [*current_renderers, ins_renderer]
    region_state.insertion_raw_sources_per_row.append(row_insertion_raw_sources)


def _format_region_coord_label(chromosome: str, start: int, end: int) -> str:
    """Return a plain-text region label for use in swappable slot coord divs."""
    size_str = format_region_size(end - start + 1)
    return f"{chromosome}:{start:,}-{end:,} ({size_str})"


def _format_original_region_size_label(start: int, end: int) -> str:
    """Return the 'Original region (size):' string for the navigation bar."""
    return f"Original region ({format_region_size(end - start + 1)}):"


def _build_region_layout_column(
    region_data,
    region_state,
    region_idx,
    hide_1bp_checkbox_override=None,
    is_swappable: bool = False,
):
    """Build column layout for one region and reset bounds. Returns (column, bounds)."""
    region = region_data["region"]
    coordinate_start = region.start
    coordinate_end = region.end
    chromosome = region.chromosome
    gene_annotations = region_data["gene_annotations"]

    hide_1bp_checkbox = None
    _orig_label_styles = {
        "font-size": "14px",
        "font-family": "Arial, sans-serif",
        "color": PLOT_CONFIG["sample_label_color"],
        "padding-right": "20px",
        "text-align": "right",
        "margin": "0",
        "padding-top": "0",
        "padding-bottom": "2px",
    }
    if is_swappable:
        orig_size_div = Div(
            text=_format_original_region_size_label(coordinate_start, coordinate_end),
            styles=_orig_label_styles,
        )
        orig_coord_div = Div(
            text=f"{chromosome}:{coordinate_start:,}-{coordinate_end:,}",
            styles=_orig_label_styles,
        )
        region_state.orig_size_div = orig_size_div
        region_state.orig_coord_div = orig_coord_div
        orig_region_col = column(orig_size_div, orig_coord_div, spacing=0)
        nav_chrom_div = Div(text=chromosome, visible=False)
        nav_orig_start_div = Div(text=str(coordinate_start), visible=False)
        nav_orig_end_div = Div(text=str(coordinate_end), visible=False)
        region_state.nav_chrom_div = nav_chrom_div
        region_state.nav_orig_start_div = nav_orig_start_div
        region_state.nav_orig_end_div = nav_orig_end_div
        coord_display, hide_1bp_checkbox, coord_input = create_coordinate_display(
            region_state.first_plot_figure,
            chromosome,
            coordinate_start,
            coordinate_end,
            read_search_button=None,
            show_checkbox_controls=False,
            model_name_suffix="" if region_idx == 0 else f"region_{region_idx}",
            original_region_widget=orig_region_col,
            nav_chrom_div=nav_chrom_div,
            nav_orig_start_div=nav_orig_start_div,
            nav_orig_end_div=nav_orig_end_div,
        )
        region_state.coord_input = coord_input
    else:
        coord_display, hide_1bp_checkbox, _ = create_coordinate_display(
            region_state.first_plot_figure,
            chromosome,
            coordinate_start,
            coordinate_end,
            read_search_button=None,
            show_checkbox_controls=False,
            model_name_suffix="" if region_idx == 0 else f"region_{region_idx}",
        )
    active_hide_1bp_checkbox = hide_1bp_checkbox or hide_1bp_checkbox_override
    for plot_fig, var_rend, one_bp in zip(
        region_state.plot_figures,
        region_state.variant_renderers,
        region_state.one_bp_by_row,
        strict=False,
    ):
        setup_variant_lod_callback(
            plot_fig,
            var_rend,
            one_bp_renderers=one_bp,
            hide_1bp_checkbox=active_hide_1bp_checkbox,
        )
    layout_components = [coord_display, *region_state.row_components]

    gene_figure = None
    if region_state.region_type != ISOSEQ_REGION_TYPE:
        gene_figure, _, gene_sources = add_gene_track_to_region(
            region_state.first_plot_figure,
            gene_annotations,
            coordinate_start,
            coordinate_end,
        )
        for gene_src in gene_sources.values():
            if gene_src is not None:
                region_state.clearable_sources.append(gene_src)
        if is_swappable:
            region_state.gene_track_sources = gene_sources
            if gene_figure is not None:
                region_state.gene_track_y_range = gene_figure.y_range
    if gene_figure:
        layout_components.append(gene_figure)
        region_state.cursor_guide_figures.append(gene_figure)

    _hide_xaxis_on_track_figures(region_state, gene_figure)

    repeat_density = region_data.get("repeat_density")
    if repeat_density is not None:
        density_figure, density_source = create_repeat_density_figure(
            region_state.shared_x_range,
            repeat_density,
            coordinate_start,
            coordinate_end,
        )
        region_state.repeat_density_source = density_source
        region_state.repeat_density_figure = density_figure
        if is_swappable:
            region_state.clearable_sources.append(density_source)
        layout_components.append(density_figure)
        region_state.cursor_guide_figures.append(density_figure)

    axis_strip = create_genomic_x_axis_strip(region_state.shared_x_range)
    region_state.cursor_guide_figures.append(axis_strip)
    layout_components.append(axis_strip)

    region_layout = column(*layout_components, sizing_mode="stretch_both", spacing=0)
    region_reset_bounds = {
        "x_range": region_state.shared_x_range,
        "x_start": coordinate_start,
        "x_end": coordinate_end,
        "y_ranges": [b[0] for b in region_state.y_bounds],
        "y_bounds": [[b[1], b[2]] for b in region_state.y_bounds],
        "all_sources": [*region_state.arrow_sources, *region_state.connector_sources],
    }
    return region_layout, region_reset_bounds, hide_1bp_checkbox


def _global_control_renderers(region_builds):
    """Collect renderers controlled by the single multi-region checkbox block."""
    one_bp_renderers = []
    one_bp_markers = []
    one_bp_texts = []
    one_bp_segments = []
    alignment_label_renderers = []
    phase_set_marker_renderers = []
    read_filter_sources = []
    read_filter_y_ranges = []
    read_filter_y_bounds = []
    region_types = set()

    for _region_idx, _region_data, region_state in region_builds:
        one_bp_renderers.extend(region_state.one_bp_renderers)
        one_bp_markers.extend(region_state.one_bp_markers)
        one_bp_texts.extend(region_state.one_bp_texts)
        one_bp_segments.extend(region_state.one_bp_segments)
        alignment_label_renderers.extend(region_state.alignment_label_renderers)
        phase_set_marker_renderers.extend(region_state.phase_set_marker_renderers)
        read_filter_sources.extend(region_state.read_filter_sources)
        read_filter_y_ranges.extend(bounds[0] for bounds in region_state.read_filter_y_bounds)
        read_filter_y_bounds.extend(bounds[1] for bounds in region_state.read_filter_y_bounds)
        region_types.add(region_state.region_type)

    return {
        "one_bp_renderers": one_bp_renderers,
        "one_bp_markers": one_bp_markers,
        "one_bp_texts": one_bp_texts,
        "one_bp_segments": one_bp_segments,
        "alignment_label_renderers": alignment_label_renderers,
        "phase_set_marker_renderers": phase_set_marker_renderers,
        "read_filter_sources": read_filter_sources,
        "read_filter_y_ranges": read_filter_y_ranges,
        "read_filter_y_bounds": read_filter_y_bounds,
        "default_hide_alignment_numbers": PARAPHASE_REGION_TYPE in region_types,
        "default_hide_1bp_indels": bool(
            {COMPLEX_SV_REGION_TYPE, ISOSEQ_REGION_TYPE} & region_types
        ),
        "enable_read_filter_checkboxes": COMPLEX_SV_REGION_TYPE in region_types,
    }


def _region_select_label(region_data: dict, index: int) -> str:
    """Return a compact display label for one region Select option."""
    region = region_data["region"]
    span = region.end - region.start + 1
    if span >= 1_000_000:
        size_str = f"{span / 1_000_000:.1f} Mb"
    elif span >= 1_000:
        size_str = f"{span / 1_000:.1f} kb"
    else:
        size_str = f"{span:,} bp"
    return f"Region {index + 1}: {region.chromosome}:{region.start:,}-{region.end:,} ({size_str})"


def _build_region_select_row(
    region_data_list: list[dict],
    region_builds: list[tuple],
    connection_index: dict,
    alignment_summaries: dict,
) -> tuple[object, object, object]:
    """Wire region-swap Select widgets and return (left_select, right_select, hidden_col).

    Called only when N > 2 regions are present. region_builds has exactly 2 entries
    (slot 0 and slot 1). All N regions are serialised to a hidden JSON div.
    Callers embed each select into its panel column for correct centering.
    """
    n_regions = len(region_data_list)
    all_options = [(str(i), _region_select_label(rd, i)) for i, rd in enumerate(region_data_list)]
    # Each panel's initial options exclude the other panel's initially-shown region.
    left_init_options = [(v, lbl) for v, lbl in all_options if v != "1"]
    right_init_options = [(v, lbl) for v, lbl in all_options if v != "0"]

    left_select = Select(
        value="0",
        options=left_init_options,
        width=340,
        styles={"font-family": "Arial, sans-serif"},
    )
    right_select = Select(
        value="1",
        options=right_init_options,
        width=340,
        styles={"font-family": "Arial, sans-serif"},
    )
    left_idx_div = Div(text="0", visible=False, name="orographer_left_idx")
    right_idx_div = Div(text="1", visible=False, name="orographer_right_idx")

    logger.debug(
        "Serialising %d regions to JSON for slot-swap (N=%d).", n_regions, n_regions
    )
    serialized_regions = {
        str(i): _serialize_region_for_swap(rd, i, connection_index, alignment_summaries)
        for i, rd in enumerate(region_data_list)
    }
    json_div = Div(
        text=json.dumps(serialized_regions, separators=(",", ":")),
        visible=False,
        name="orographer_region_json",
    )

    _region_idx0, _region_data0, state0 = region_builds[0]
    _region_idx1, _region_data1, state1 = region_builds[1]

    swap_args = {
        "left_select": left_select,
        "right_select": right_select,
        "left_idx_div": left_idx_div,
        "right_idx_div": right_idx_div,
        "json_div": json_div,
        "all_options": all_options,
        "s0_arrow_sources": state0.arrow_sources,
        "s0_marker_connector_sources": state0.marker_connector_sources,
        "s0_arc_connector_sources": state0.arc_connector_sources,
        "s0_coverage_total_sources": state0.coverage_total_sources,
        "s0_coverage_hp_sources": state0.coverage_hp_sources,
        "s0_vcf_sources": state0.vcf_sources,
        "s0_clearable_sources": state0.clearable_sources,
        "s0_x_range": state0.shared_x_range,
        "s0_y_ranges": [b[0] for b in state0.y_bounds],
        "s0_orig_size_div": state0.orig_size_div,
        "s0_orig_coord_div": state0.orig_coord_div,
        "s0_dom_classes": state0.slot_dom_classes,
        "s0_plot_model_ids": state0.plot_model_ids,
        "s1_arrow_sources": state1.arrow_sources,
        "s1_marker_connector_sources": state1.marker_connector_sources,
        "s1_arc_connector_sources": state1.arc_connector_sources,
        "s1_coverage_total_sources": state1.coverage_total_sources,
        "s1_coverage_hp_sources": state1.coverage_hp_sources,
        "s1_vcf_sources": state1.vcf_sources,
        "s1_clearable_sources": state1.clearable_sources,
        "s1_x_range": state1.shared_x_range,
        "s1_y_ranges": [b[0] for b in state1.y_bounds],
        "s1_orig_size_div": state1.orig_size_div,
        "s1_orig_coord_div": state1.orig_coord_div,
        "s1_dom_classes": state1.slot_dom_classes,
        "s1_plot_model_ids": state1.plot_model_ids,
        # Alignment (read number) label sources per bam-row (None for rows without labels)
        "s0_alignment_label_sources": state0.alignment_label_source_per_row,
        "s1_alignment_label_sources": state1.alignment_label_source_per_row,
        # HP label sources per bam-row (None for rows without labels/separators)
        "s0_hp_sep_sources": [r.get("separator") for r in state0.hp_label_sources_per_row],
        "s0_hp_lbl_sources": [r.get("label") for r in state0.hp_label_sources_per_row],
        "s1_hp_sep_sources": [r.get("separator") for r in state1.hp_label_sources_per_row],
        "s1_hp_lbl_sources": [r.get("label") for r in state1.hp_label_sources_per_row],
        # Gene track ColumnDataSources (one set per slot, keyed by type)
        "s0_gene_body": state0.gene_track_sources.get("body"),
        "s0_gene_exon": state0.gene_track_sources.get("exon"),
        "s0_gene_intron": state0.gene_track_sources.get("intron"),
        "s0_gene_arrow": state0.gene_track_sources.get("arrow"),
        "s0_gene_label": state0.gene_track_sources.get("label"),
        "s1_gene_body": state1.gene_track_sources.get("body"),
        "s1_gene_exon": state1.gene_track_sources.get("exon"),
        "s1_gene_intron": state1.gene_track_sources.get("intron"),
        "s1_gene_arrow": state1.gene_track_sources.get("arrow"),
        "s1_gene_label": state1.gene_track_sources.get("label"),
        # Coordinate-navigation models (updated after swap so nav bounds stay correct)
        "s0_coord_input": state0.coord_input,
        "s0_nav_chrom_div": state0.nav_chrom_div,
        "s0_nav_orig_start_div": state0.nav_orig_start_div,
        "s0_nav_orig_end_div": state0.nav_orig_end_div,
        "s1_coord_input": state1.coord_input,
        "s1_nav_chrom_div": state1.nav_chrom_div,
        "s1_nav_orig_start_div": state1.nav_orig_start_div,
        "s1_nav_orig_end_div": state1.nav_orig_end_div,
        # Separator-line source per bam-row (horizontal dotted lines between reads)
        "s0_separator_sources": state0.separator_source_per_row,
        "s1_separator_sources": state1.separator_source_per_row,
        # Phase-set boundary source per bam-row (vertical boundary lines; None if no phase data)
        "s0_phase_set_sources": state0.phase_set_source_per_row,
        "s1_phase_set_sources": state1.phase_set_source_per_row,
        # Variant sources per bam-row (nested list; variable length per row)
        "s0_variant_sources_per_row": state0.variant_sources_per_row,
        "s1_variant_sources_per_row": state1.variant_sources_per_row,
        # Insertion-summary raw sources per bam-row (nested list; one per HP)
        "s0_insertion_raw_sources_per_row": state0.insertion_raw_sources_per_row,
        "s1_insertion_raw_sources_per_row": state1.insertion_raw_sources_per_row,
        # Repeat-density (Ident) track source and figure y_range per slot
        "s0_repeat_density_source": state0.repeat_density_source,
        "s0_repeat_density_y_range": (
            state0.repeat_density_figure.y_range if state0.repeat_density_figure else None
        ),
        "s1_repeat_density_source": state1.repeat_density_source,
        "s1_repeat_density_y_range": (
            state1.repeat_density_figure.y_range if state1.repeat_density_figure else None
        ),
        # Gene track figure y_range per slot (None when region has no gene annotations)
        "s0_gene_track_y_range": state0.gene_track_y_range,
        "s1_gene_track_y_range": state1.gene_track_y_range,
        # Coverage figure y_range per bam-row per slot (None for rows without coverage)
        "s0_coverage_y_ranges": state0.coverage_y_ranges_per_row,
        "s1_coverage_y_ranges": state1.coverage_y_ranges_per_row,
    }
    swap_callback = CustomJS(
        args=swap_args,
        code=load_javascript("region_swap_callback.js"),
    )
    left_select.js_on_change("value", swap_callback)
    right_select.js_on_change("value", swap_callback)

    hidden_col = column(
        json_div,
        left_idx_div,
        right_idx_div,
        visible=False,
        sizing_mode="fixed",
        width=0,
        height=0,
    )
    return left_select, right_select, hidden_col


def plot_reads_bokeh(region_data_list, output_config):
    """
    Create Bokeh HTML plot for one or more regions with cross-region read highlighting.
    Each region has bam_rows (1-3 BAMs); stacked top to bottom (other, other, primary).
    output_config has .output_dir and .prefix attributes.
    """
    if not region_data_list:
        logger.debug("No regions to plot.")
        return None
    output_file = generate_multi_region_filename(
        region_data_list,
        output_config.output_dir,
        output_config.prefix,
        filename_regions=output_config.filename_regions,
    )
    has_isoseq = any(
        bam_row.get("region_type") == ISOSEQ_REGION_TYPE
        for region_data in region_data_list
        for bam_row in region_data.get("bam_rows", [])
    )
    all_isoseq = has_isoseq and all(
        bam_row.get("region_type") == ISOSEQ_REGION_TYPE
        for region_data in region_data_list
        for bam_row in region_data.get("bam_rows", [])
    )
    isoseq_chunk_dir = None
    isoseq_chunk_url_prefix = None
    if has_isoseq:
        output_stem = os.path.splitext(os.path.basename(output_file))[0]
        chunk_root = f"{output_stem}_chunks"
        chunk_token = uuid.uuid4().hex[:12]
        isoseq_chunk_url_prefix = f"{chunk_root}/{chunk_token}"
        isoseq_chunk_dir = os.path.join(os.path.dirname(output_file), chunk_root, chunk_token)
    all_arrow_sources = []
    all_arrow_renderers = []
    all_selectable_sources = []
    all_connector_sources = []
    all_transcript_sources = []
    all_isoseq_filter_sources = []
    all_isoseq_transcript_label_sources = []
    all_isoseq_intron_sources = []
    all_isoseq_intron_arrow_sources = []
    all_isoseq_feature_label_sources = []
    all_isoseq_comparison_components = []
    all_alignment_label_sources = []
    all_region_layouts = []
    all_plot_figures = []
    all_cursor_guide_figures = []
    all_cursor_guide_spans = []
    all_region_reset_bounds = []
    region_builds = []
    connection_index = {}
    alignment_summaries = {}
    if not all_isoseq:
        assign_region_agnostic_alignment_orders(region_data_list)
        connection_index = build_read_connection_index(region_data_list)
        alignment_summaries = build_read_alignment_summaries(region_data_list)

    n_regions = len(region_data_list)
    use_slot_swap = n_regions > 2 and not has_isoseq
    build_limit = 2 if use_slot_swap else n_regions

    for region_idx, region_data in enumerate(region_data_list):
        if region_idx >= build_limit:
            break
        region = region_data["region"]
        coordinate_start = region.start
        coordinate_end = region.end
        bam_rows = region_data["bam_rows"]

        region_state = RegionBuildState()

        for row_index, bam_row in enumerate(bam_rows):
            _build_bam_row_track(
                bam_row,
                region,
                coordinate_start,
                coordinate_end,
                region_state,
                region_idx,
                row_index,
                connection_index=connection_index,
                alignment_summaries=alignment_summaries,
                isoseq_chunk_dir=isoseq_chunk_dir,
                isoseq_chunk_url_prefix=isoseq_chunk_url_prefix,
            )

        if not region_state.first_plot_figure:
            logger.debug(f"No segments to plot for region {region_idx + 1}.")
            continue

        all_arrow_sources.extend(region_state.arrow_sources)
        all_arrow_renderers.extend(region_state.arrow_renderers)
        all_selectable_sources.extend(region_state.arrow_sources)
        all_selectable_sources.extend(region_state.connector_sources)
        all_connector_sources.extend(region_state.connector_sources)
        all_transcript_sources.extend(region_state.transcript_sources)
        all_isoseq_filter_sources.extend(region_state.isoseq_filter_sources)
        all_isoseq_transcript_label_sources.extend(region_state.isoseq_transcript_label_sources)
        all_isoseq_intron_sources.extend(region_state.isoseq_intron_sources)
        all_isoseq_intron_arrow_sources.extend(region_state.isoseq_intron_arrow_sources)
        all_isoseq_feature_label_sources.extend(region_state.isoseq_feature_label_sources)
        all_isoseq_comparison_components.extend(region_state.isoseq_comparison_components)
        all_alignment_label_sources.extend(region_state.alignment_label_sources)
        all_plot_figures.extend(region_state.plot_figures)
        region_builds.append((region_idx, region_data, region_state))

    global_control_renderers = _global_control_renderers(region_builds)
    global_checkbox_controls = None
    global_hide_1bp_checkbox = None
    global_cursor_guide_checkbox = None
    if region_builds:
        _first_region_idx, first_region_data, first_region_state = region_builds[0]
        shared_dotplot_payload = region_data_list[0].get("dotplot_payload")
        dotplot_thumbnail = None
        if shared_dotplot_payload is not None:
            dotplot_thumbnail = create_dotplot_thumbnail(
                first_region_state.first_plot_figure,
                shared_dotplot_payload["matrix"],
                first_region_data["region"],
                size=18,
                blocks=shared_dotplot_payload["blocks"],
                region_label=shared_dotplot_payload["label"],
                modal_title=shared_dotplot_payload["title"],
                individual_payloads=shared_dotplot_payload.get("individual_payloads"),
            )
        read_search_controls = (
            create_read_search_button(all_selectable_sources) if all_selectable_sources else None
        )
        isoseq_controls = create_isoseq_controls(
            all_transcript_sources,
            [
                *all_arrow_sources,
                *all_alignment_label_sources,
                *all_isoseq_filter_sources,
            ],
            transcript_label_sources=all_isoseq_transcript_label_sources,
            intron_sources=all_isoseq_intron_sources,
            intron_arrow_sources=all_isoseq_intron_arrow_sources,
            feature_label_sources=all_isoseq_feature_label_sources,
            read_arrow_sources=all_arrow_sources,
            read_label_sources=all_alignment_label_sources,
            plot_figures=all_plot_figures,
            comparison_components=all_isoseq_comparison_components,
            dotplot_thumbnail=dotplot_thumbnail,
        )
        if read_search_controls is not None and isoseq_controls is not None:
            read_search_width = getattr(read_search_controls, "width", 0) or 0
            isoseq_width = getattr(isoseq_controls, "width", 0) or 0
            read_search_controls = row(
                read_search_controls,
                isoseq_controls,
                spacing=14,
                sizing_mode="fixed",
                width=read_search_width + isoseq_width + 14,
                height=30,
            )
        elif isoseq_controls is not None:
            read_search_controls = isoseq_controls
        (
            global_checkbox_controls,
            global_hide_1bp_checkbox,
            global_cursor_guide_checkbox,
        ) = create_global_checkbox_controls(
            first_region_state.first_plot_figure,
            dotplot_thumbnail=None if isoseq_controls is not None else dotplot_thumbnail,
            read_search_controls=read_search_controls,
            one_bp_renderers=global_control_renderers["one_bp_renderers"],
            one_bp_markers=global_control_renderers["one_bp_markers"],
            one_bp_texts=global_control_renderers["one_bp_texts"],
            one_bp_segments=global_control_renderers["one_bp_segments"],
            alignment_label_renderers=global_control_renderers["alignment_label_renderers"],
            phase_set_marker_renderers=global_control_renderers["phase_set_marker_renderers"],
            default_hide_alignment_numbers=global_control_renderers[
                "default_hide_alignment_numbers"
            ],
            default_hide_1bp_indels=global_control_renderers["default_hide_1bp_indels"],
            enable_read_filter_checkboxes=global_control_renderers["enable_read_filter_checkboxes"],
            read_filter_sources=global_control_renderers["read_filter_sources"],
            read_filter_y_ranges=global_control_renderers["read_filter_y_ranges"],
            read_filter_y_bounds=global_control_renderers["read_filter_y_bounds"],
        )
    for _layout_idx, (region_idx, region_data, region_state) in enumerate(region_builds):
        region_layout, region_reset_bounds, _hide_1bp_checkbox = _build_region_layout_column(
            region_data,
            region_state,
            region_idx,
            hide_1bp_checkbox_override=global_hide_1bp_checkbox,
            is_swappable=use_slot_swap,
        )
        all_region_layouts.append(region_layout)
        all_cursor_guide_figures.extend(region_state.cursor_guide_figures)
        all_cursor_guide_spans.extend(
            add_cursor_guide_to_figures(region_state.cursor_guide_figures)
        )
        all_region_reset_bounds.append(
            (
                region_state.plot_figures,
                region_state.arrow_sources,
                region_state.arrow_renderers,
                region_reset_bounds,
            )
        )

    if not all_region_layouts:
        no_data_div = Div(
            text="No alignments in the requested region(s).",
            styles={"color": "#666", "font-size": "14px", "padding": "20px"},
        )
        all_region_layouts = [column(no_data_div, sizing_mode="stretch_both")]

    hidden_divs_col = None
    if use_slot_swap and len(region_builds) == 2:
        left_select, right_select, hidden_divs_col = _build_region_select_row(
            region_data_list, region_builds, connection_index, alignment_summaries
        )
        # Embed each select centered inside its panel column.
        for sel, layout in zip(
            (left_select, right_select), all_region_layouts[:2], strict=False
        ):
            centered = row(
                Spacer(sizing_mode="stretch_width"),
                sel,
                Spacer(sizing_mode="stretch_width"),
                sizing_mode="stretch_width",
            )
            all_region_layouts[all_region_layouts.index(layout)] = column(
                centered, layout, sizing_mode="stretch_both", spacing=4
            )

    region_row = row(*all_region_layouts, sizing_mode="stretch_both", spacing=20)

    layout_parts = []
    if global_checkbox_controls is not None:
        layout_parts.append(global_checkbox_controls)
    if hidden_divs_col is not None:
        layout_parts.append(hidden_divs_col)
    layout_parts.append(region_row)

    if len(layout_parts) == 1:
        final_layout = layout_parts[0]
    else:
        final_layout = column(*layout_parts, sizing_mode="stretch_both", spacing=0)
    add_multi_region_callbacks(
        all_plot_figures,
        all_arrow_sources,
        all_arrow_renderers,
        all_selectable_sources,
        all_region_reset_bounds,
    )
    add_read_connection_overlay_callbacks(
        all_plot_figures,
        all_connector_sources,
        all_selectable_sources,
    )
    add_alignment_number_modal_callbacks(all_alignment_label_sources)
    add_alignment_label_selection_callbacks(all_selectable_sources, all_alignment_label_sources)
    add_cursor_guide_callbacks(
        all_cursor_guide_figures,
        all_cursor_guide_spans,
        global_cursor_guide_checkbox,
    )
    save_plot_with_modal(final_layout, output_file, output_config.prefix)
    return output_file
