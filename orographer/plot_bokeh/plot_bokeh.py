"""
Orchestration and track content for Bokeh plots
(colors, segments, arrows, variants, entry point).
"""

import logging
import math

from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Div, HoverTool, TapTool

from ..utils import PARAPHASE_REGION_TYPE
from .callbacks import (
    get_arrow_tap_callback,
    get_arrow_tap_callback_multi_region,
    get_number_click_callback,
    get_vcf_variant_click_callback,
    save_plot_with_modal,
)
from .data import (
    calculate_read_positions,
    generate_multi_region_filename,
    sort_read_names,
)
from .figures import (
    add_gene_track,
    add_haplotype_labels,
    add_separator_lines,
    create_bokeh_figure,
    create_bokeh_figure_shared_x,
    create_coordinate_display,
    create_gene_track_figure,
    create_genomic_x_axis_strip,
    create_vcf_track_figure,
)
from .utils import PLOT_CONFIG, RegionBuildState, load_javascript

logger = logging.getLogger(__name__)
logging.getLogger("bokeh").setLevel(logging.WARNING)


def get_segment_color(segment, region_type):
    """Color for a segment (YC tag or haplotype/strand; paraphase no YC -> grey)."""
    if hasattr(segment, "color_tag") and segment.color_tag:
        try:
            rgb = segment.color_tag.split(",")
            if len(rgb) == 3:
                r, g, b = int(rgb[0]), int(rgb[1]), int(rgb[2])
                return f"#{r:02x}{g:02x}{b:02x}"
        except (ValueError, IndexError):
            pass
    # Paraphase: reads without YC tag are grey regardless of strand
    if region_type == PARAPHASE_REGION_TYPE:
        return PLOT_CONFIG["segment_paraphase_color"]
    haplotype = segment.haplotype_tag
    if haplotype is None or haplotype == 0:
        return PLOT_CONFIG["segment_unassigned_color"]
    return (
        PLOT_CONFIG["segment_fwd_color"]
        if segment.is_fwd_strand
        else PLOT_CONFIG["segment_rev_color"]
    )


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


def process_segments(
    segments_by_read,
    read_names,
    read_to_y,
    read_heights,
    coordinate_start,
    coordinate_end,
    region_type,
    sample_label=None,
):
    """Process segments and extract arrow, clickable label, and variant data."""
    arrow_x0_list, arrow_x1_list, arrow_y_list, arrow_color_list, arrow_read_names = (
        [],
        [],
        [],
        [],
        [],
    )
    clickable_x, clickable_y, clickable_customdata = [], [], []
    mismatch_x, mismatch_y, mismatch_alt, mismatch_color = [], [], [], []
    insertion_x, insertion_y, insertion_size, insertion_count, insertion_is_1bp = (
        [],
        [],
        [],
        [],
        [],
    )
    deletion_x0, deletion_x1, deletion_y, deletion_is_1bp = [], [], [], []

    for read_name in read_names:
        segments = segments_by_read[read_name]
        y_pos = read_to_y[read_name]
        read_height = read_heights[read_name]
        num_segments = len(segments)
        for idx, segment in enumerate(segments):
            start, end = segment.pos, segment.end
            if start < coordinate_end and end > coordinate_start:
                plot_start = max(start, coordinate_start)
                plot_end = min(end, coordinate_end)
                color = get_segment_color(segment, region_type=region_type)
                if num_segments > 1:
                    spacing = read_height / (num_segments + 1)
                    y_offset = -read_height / 2 + spacing * (idx + 1)
                else:
                    y_offset = 0
                y_plot = y_pos + y_offset
                if segment.is_fwd_strand:
                    arrow_x0, arrow_x1 = plot_start, plot_end
                else:
                    arrow_x0, arrow_x1 = plot_end, plot_start
                mid_x = (plot_start + plot_end) / 2
                if abs(arrow_x1 - arrow_x0) > 0:
                    arrow_x0_list.append(arrow_x0)
                    arrow_x1_list.append(arrow_x1)
                    arrow_y_list.append(y_plot)
                    arrow_color_list.append(color)
                    arrow_read_names.append(segment.readname)

                if hasattr(segment, "mismatches") and segment.mismatches:
                    for ref_pos, alt_base in segment.mismatches:
                        ref_pos_1based = ref_pos + 1
                        if coordinate_start <= ref_pos_1based <= coordinate_end:
                            mismatch_x.append(ref_pos_1based)
                            mismatch_y.append(y_plot)
                            mismatch_alt.append(alt_base)
                            mismatch_color.append(get_base_color(alt_base))
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

                strand_str = "Forward (+)" if segment.is_fwd_strand else "Reverse (-)"
                haplotype_str = (
                    f"HP:{segment.haplotype_tag}" if segment.haplotype_tag else "Unassigned"
                )
                clickable_x.append(mid_x)
                clickable_y.append(y_plot)
                alignment_number = getattr(segment, "alignment_order", 0) or (idx + 1)
                clickable_customdata.append(
                    {
                        "read_name": segment.readname,
                        "alignment_number": alignment_number,
                        "strand": strand_str,
                        "coordinates": f"{segment.chrom}:{segment.pos}-{segment.end}",
                        "haplotype": haplotype_str,
                        "sample_label": sample_label or "",
                    }
                )

    return (
        {
            "x0": arrow_x0_list,
            "x1": arrow_x1_list,
            "y": arrow_y_list,
            "color": arrow_color_list,
            "read_name": arrow_read_names,
        },
        {"x": clickable_x, "y": clickable_y, "customdata": clickable_customdata},
        {
            "mismatch": {
                "x": mismatch_x,
                "y": mismatch_y,
                "alt": mismatch_alt,
                "color": mismatch_color,
            },
            "insertion": {
                "x": insertion_x,
                "y": insertion_y,
                "size": insertion_size,
                "count": insertion_count,
                "is_1bp": insertion_is_1bp,
            },
            "deletion": {
                "x0": deletion_x0,
                "x1": deletion_x1,
                "y": deletion_y,
                "is_1bp": deletion_is_1bp,
            },
        },
    )


def add_arrows_to_plot(plot_figure, arrow_data):
    """Add arrow segments and arrowheads to the plot using batched glyphs."""
    if not arrow_data["x0"]:
        return None
    arrow_source = ColumnDataSource(data=arrow_data)
    arrow_renderer = plot_figure.segment(
        x0="x0",
        y0="y",
        x1="x1",
        y1="y",
        source=arrow_source,
        line_color="color",
        line_width=PLOT_CONFIG["arrow_line_width"],
        line_alpha=PLOT_CONFIG["arrow_line_alpha"],
        selection_line_color=PLOT_CONFIG["arrow_selection_line_color"],
        selection_line_alpha=1,
        selection_line_width=PLOT_CONFIG["arrow_selection_line_width"],
        nonselection_line_alpha=PLOT_CONFIG["arrow_nonselection_line_alpha"],
        nonselection_line_width=PLOT_CONFIG["arrow_nonselection_line_width"],
    )
    arrowhead_angles = []
    for i in range(len(arrow_data["x0"])):
        x0, x1 = arrow_data["x0"][i], arrow_data["x1"][i]
        arrowhead_angles.append(-math.pi / 2 if x1 > x0 else math.pi / 2)
    arrowhead_source = ColumnDataSource(
        data={
            "x": arrow_data["x1"],
            "y": arrow_data["y"],
            "color": arrow_data["color"],
            "angle": arrowhead_angles,
        }
    )
    plot_figure.scatter(
        x="x",
        y="y",
        source=arrowhead_source,
        marker="triangle",
        size=PLOT_CONFIG["arrowhead_size"],
        angle="angle",
        fill_color="color",
        fill_alpha=PLOT_CONFIG["arrowhead_fill_alpha"],
        line_color="color",
        line_alpha=PLOT_CONFIG["arrowhead_line_alpha"],
        line_width=PLOT_CONFIG["arrowhead_line_width"],
    )
    return arrow_source, arrow_renderer


def add_variants_to_plot(plot_figure, variant_data):
    """Add variant markers (mismatches, insertions, deletions); return LOD/1bp."""
    renderers = {
        "marker": [],
        "text": [],
        "one_bp": [],
        "one_bp_markers": [],
        "one_bp_texts": [],
        "one_bp_segments": [],
    }
    mismatch_data = variant_data["mismatch"]
    if mismatch_data["x"]:
        mismatch_source = ColumnDataSource(data=mismatch_data)
        mismatch_marker = plot_figure.scatter(
            x="x",
            y="y",
            source=mismatch_source,
            marker="square",
            size=PLOT_CONFIG["mismatch_size"],
            fill_color="color",
            fill_alpha=PLOT_CONFIG["mismatch_fill_alpha"],
            line_color="color",
            line_alpha=PLOT_CONFIG["mismatch_line_alpha"],
            line_width=PLOT_CONFIG["mismatch_line_width"],
        )
        renderers["marker"].append(mismatch_marker)
        mismatch_text = plot_figure.text(
            x="x",
            y="y",
            text="alt",
            source=mismatch_source,
            text_font_size=PLOT_CONFIG["mismatch_text_font_size"],
            text_color="color",
            text_align="center",
            text_baseline="middle",
            text_font_style="bold",
        )
        mismatch_text.visible = False
        renderers["text"].append(mismatch_text)
    insertion_data = variant_data["insertion"]
    if insertion_data["x"]:
        insertion_is_1bp = insertion_data["is_1bp"]
        idx_1bp = [i for i, b in enumerate(insertion_is_1bp) if b]
        idx_other = [i for i, b in enumerate(insertion_is_1bp) if not b]
        insertion_color = PLOT_CONFIG["variant_color_insertion"]
        for idx_list, is_one_bp in [(idx_1bp, True), (idx_other, False)]:
            if not idx_list:
                continue
            sub = {k: [v[i] for i in idx_list] for k, v in insertion_data.items() if k != "is_1bp"}
            sub_source = ColumnDataSource(data=sub)
            marker = plot_figure.scatter(
                x="x",
                y="y",
                source=sub_source,
                marker="square",
                size="size",
                fill_color=insertion_color,
                fill_alpha=PLOT_CONFIG["insertion_fill_alpha"],
                line_color=insertion_color,
                line_alpha=PLOT_CONFIG["insertion_line_alpha"],
                line_width=PLOT_CONFIG["insertion_line_width"],
            )
            renderers["marker"].append(marker)
            if is_one_bp:
                renderers["one_bp"].append(marker)
                renderers["one_bp_markers"].append(marker)
            text_data = {
                "x": sub["x"],
                "y": sub["y"],
                "text": [f"{sub['count'][j]}I" for j in range(len(sub["x"]))],
            }
            text_source = ColumnDataSource(data=text_data)
            text = plot_figure.text(
                x="x",
                y="y",
                text="text",
                source=text_source,
                text_font_size=PLOT_CONFIG["insertion_text_font_size"],
                text_color=insertion_color,
                text_align="center",
                text_baseline="middle",
                text_font_style="bold",
            )
            text.visible = False
            renderers["text"].append(text)
            if is_one_bp:
                renderers["one_bp"].append(text)
                renderers["one_bp_texts"].append(text)
    deletion_data = variant_data["deletion"]
    if deletion_data["x0"]:
        deletion_is_1bp = deletion_data["is_1bp"]
        idx_1bp = [i for i, b in enumerate(deletion_is_1bp) if b]
        idx_other = [i for i, b in enumerate(deletion_is_1bp) if not b]
        for idx_list, is_one_bp in [(idx_1bp, True), (idx_other, False)]:
            if not idx_list:
                continue
            sub = {
                "x0": [deletion_data["x0"][i] for i in idx_list],
                "x1": [deletion_data["x1"][i] for i in idx_list],
                "y": [deletion_data["y"][i] for i in idx_list],
            }
            sub_source = ColumnDataSource(data=sub)
            seg = plot_figure.segment(
                x0="x0",
                y0="y",
                x1="x1",
                y1="y",
                source=sub_source,
                line_color=PLOT_CONFIG["variant_color_deletion"],
                line_width=PLOT_CONFIG["deletion_line_width"],
                line_alpha=PLOT_CONFIG["deletion_line_alpha"],
            )
            if is_one_bp:
                renderers["one_bp"].append(seg)
                renderers["one_bp_segments"].append(seg)
    return renderers


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


def add_clickable_labels(
    plot_figure,
    tap_tool,
    clickable_data,
    arrow_source=None,
    arrow_renderer=None,
    arrow_tap_callback=None,
    reset_callback=None,
):
    """Add clickable number labels with modal; optional arrow tap/reset callbacks.
    Returns list of renderers for visibility toggling, or None if no labels.
    """
    if not clickable_data["x"]:
        return None
    source = ColumnDataSource(
        data={
            "x": clickable_data["x"],
            "y": clickable_data["y"],
            "text": [str(info["alignment_number"]) for info in clickable_data["customdata"]],
            "read_name": [info["read_name"] for info in clickable_data["customdata"]],
            "alignment_number": [info["alignment_number"] for info in clickable_data["customdata"]],
            "strand": [info["strand"] for info in clickable_data["customdata"]],
            "coordinates": [info["coordinates"] for info in clickable_data["customdata"]],
            "haplotype": [info["haplotype"] for info in clickable_data["customdata"]],
            "sample_label": [info.get("sample_label", "") for info in clickable_data["customdata"]],
        }
    )
    cfg = PLOT_CONFIG
    visible_circles = plot_figure.scatter(
        "x",
        "y",
        source=source,
        size=cfg["alignment_label_visible_size"],
        marker="circle",
        fill_color=cfg["alignment_label_fill_color"],
        line_color=cfg["alignment_label_fill_color"],
        alpha=0.8,
    )
    circles = plot_figure.scatter(
        "x",
        "y",
        source=source,
        size=cfg["alignment_label_hit_size"],
        marker="circle",
        fill_color=cfg["alignment_label_fill_color"],
        line_color=cfg["alignment_label_fill_color"],
        alpha=0,
    )
    tap_tool.renderers = [circles]
    if arrow_source is not None and arrow_renderer is not None:
        if arrow_tap_callback is None:
            arrow_tap_callback = get_arrow_tap_callback(arrow_source)
            plot_figure.js_on_event("tap", arrow_tap_callback)
        if reset_callback is None:
            reset_callback = CustomJS(
                args={"source": arrow_source},
                code=load_javascript("arrow_reset_callback.js"),
            )
        plot_figure.js_on_event("reset", reset_callback)
    label_text = plot_figure.text(
        "x",
        "y",
        text="text",
        source=source,
        text_font_size=cfg["alignment_label_text_font_size"],
        text_color=cfg["alignment_label_text_color"],
        text_align="center",
        text_baseline="middle",
        text_font_style="bold",
    )
    number_click_callback = get_number_click_callback(source)
    source.selected.js_on_change("indices", number_click_callback)
    plot_figure.add_tools(
        HoverTool(renderers=[circles], tooltips=[("Alignment #", "@text")], visible=False)
    )
    return [visible_circles, circles, label_text]


def add_vcf_track_to_region(
    plot_figure, vcf_variants, coordinate_start, coordinate_end, sample_label=None
):
    """Add VCF variant track to a region plot if variants are provided."""
    if not vcf_variants:
        return None
    vcf_figure = create_vcf_track_figure(plot_figure)
    add_vcf_variants(
        vcf_figure,
        vcf_variants,
        coordinate_start,
        coordinate_end,
        sample_label=sample_label,
    )
    return vcf_figure


def add_gene_track_to_region(plot_figure, gene_annotations, coordinate_start, coordinate_end):
    """Add gene annotation track to a region plot if annotations are provided."""
    if not gene_annotations:
        return None, plot_figure
    gene_track_height = 4.0
    gene_figure = create_gene_track_figure(plot_figure, gene_track_height)
    actual_gene_track_height = add_gene_track(
        gene_figure, gene_annotations, 0, coordinate_start, coordinate_end
    )
    label_padding = 0.5
    if actual_gene_track_height > gene_track_height:
        gene_figure.y_range.start = actual_gene_track_height
    gene_figure.y_range.end = -label_padding
    return gene_figure, plot_figure


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
    all_region_reset_bounds=None,
):
    """Add cross-region highlighting callbacks and per-region range reset."""
    if len(all_arrow_sources) == 0:
        return
    for plot_fig, arrow_source, _arrow_renderer in zip(
        all_plot_figures, all_arrow_sources, all_arrow_renderers, strict=False
    ):
        if len(all_arrow_sources) == 1:
            single_region_callback = get_arrow_tap_callback(arrow_source)
            plot_fig.js_on_event("tap", single_region_callback)
        else:
            multi_region_callback = get_arrow_tap_callback_multi_region(
                arrow_source, all_arrow_sources
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
            args={"all_sources": all_arrow_sources},
            code=load_javascript("shared_reset_callback.js"),
        )
        for plot_fig in all_plot_figures:
            plot_fig.js_on_event("reset", shared_reset_callback)


def _build_bam_row_track(
    bam_row, coordinate_start, coordinate_end, region_state, region_idx, row_index
):
    """Add one BAM row to the region. Mutates region_state. Skips if no segments."""
    segments_by_read = bam_row["segments_by_read"]
    if not segments_by_read:
        logger.debug(f"No segments to plot for region {region_idx + 1}, row {row_index + 1}.")
        return

    vcf_variants = bam_row.get("vcf_variants")
    region_type = bam_row["region_type"]
    sample_label = bam_row.get("sample_label")

    read_names, haplotype_groups, haplotype_order = sort_read_names(segments_by_read)
    read_to_y, read_to_y_bottom, read_heights, alignments_height, group_boundaries = (
        calculate_read_positions(read_names, segments_by_read, haplotype_groups)
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

    vcf_figure = add_vcf_track_to_region(
        plot_figure,
        vcf_variants,
        coordinate_start,
        coordinate_end,
        sample_label=sample_label,
    )
    if vcf_figure:
        region_state.row_components.append(vcf_figure)
        region_state.vcf_figures.append(vcf_figure)
    else:
        region_state.vcf_figures.append(None)

    add_separator_lines(
        plot_figure,
        read_names,
        read_to_y_bottom,
        read_heights,
        coordinate_start,
        coordinate_end,
    )
    arrow_data, clickable_data, variant_data = process_segments(
        segments_by_read,
        read_names,
        read_to_y,
        read_heights,
        coordinate_start,
        coordinate_end,
        region_type,
        sample_label=sample_label,
    )
    arrow_result = add_arrows_to_plot(plot_figure, arrow_data)
    arrow_source = arrow_result[0] if arrow_result else None
    arrow_renderer = arrow_result[1] if arrow_result else None
    if arrow_source:
        region_state.arrow_sources.append(arrow_source)
    if arrow_renderer:
        region_state.arrow_renderers.append(arrow_renderer)
    top_padding = max(alignments_height * 0.02, 0.5)
    region_state.y_bounds.append((plot_figure.y_range, alignments_height, -top_padding))
    # Vertical scroll when zoomed is handled by WheelPanTool(dimensions="height") on the figure
    variant_renderers = add_variants_to_plot(plot_figure, variant_data)
    row_one_bp = variant_renderers.get("one_bp", [])
    region_state.one_bp_renderers.extend(row_one_bp)
    region_state.one_bp_markers.extend(variant_renderers.get("one_bp_markers", []))
    region_state.one_bp_texts.extend(variant_renderers.get("one_bp_texts", []))
    region_state.one_bp_segments.extend(variant_renderers.get("one_bp_segments", []))
    region_state.variant_renderers.append(variant_renderers)
    region_state.one_bp_by_row.append(row_one_bp)
    label_renderers = add_clickable_labels(
        plot_figure,
        tap_tool,
        clickable_data,
        arrow_source,
        arrow_renderer,
        arrow_tap_callback=CustomJS(code=""),
    )
    if label_renderers:
        region_state.alignment_label_renderers.extend(label_renderers)
    add_haplotype_labels(
        plot_figure, group_boundaries, haplotype_order, coordinate_start, coordinate_end
    )
    region_state.row_components.append(plot_figure)
    region_state.plot_figures.append(plot_figure)


def _build_region_layout_column(region_data, region_state):
    """Build column layout for one region and reset bounds. Returns (column, bounds)."""
    region = region_data["region"]
    coordinate_start = region.start
    coordinate_end = region.end
    chromosome = region.chromosome
    gene_annotations = region_data["gene_annotations"]

    coord_div, hide_1bp_checkbox = create_coordinate_display(
        region_state.first_plot_figure,
        chromosome,
        coordinate_start,
        coordinate_end,
        one_bp_renderers=region_state.one_bp_renderers,
        one_bp_markers=region_state.one_bp_markers,
        one_bp_texts=region_state.one_bp_texts,
        one_bp_segments=region_state.one_bp_segments,
        alignment_label_renderers=region_state.alignment_label_renderers,
        default_hide_alignment_numbers=(region_state.region_type == PARAPHASE_REGION_TYPE),
    )
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
            hide_1bp_checkbox=hide_1bp_checkbox,
        )
    layout_components = [coord_div, *region_state.row_components]

    gene_figure, _ = add_gene_track_to_region(
        region_state.first_plot_figure,
        gene_annotations,
        coordinate_start,
        coordinate_end,
    )
    if gene_figure:
        layout_components.append(gene_figure)

    _hide_xaxis_on_track_figures(region_state, gene_figure)
    layout_components.append(create_genomic_x_axis_strip(region_state.shared_x_range))

    region_layout = column(*layout_components, sizing_mode="stretch_both", spacing=0)
    region_reset_bounds = {
        "x_range": region_state.shared_x_range,
        "x_start": coordinate_start,
        "x_end": coordinate_end,
        "y_ranges": [b[0] for b in region_state.y_bounds],
        "y_bounds": [[b[1], b[2]] for b in region_state.y_bounds],
        "all_sources": region_state.arrow_sources,
    }
    return region_layout, region_reset_bounds


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
        region_data_list, output_config.output_dir, output_config.prefix
    )
    all_arrow_sources = []
    all_arrow_renderers = []
    all_region_layouts = []
    all_plot_figures = []
    all_region_reset_bounds = []

    for region_idx, region_data in enumerate(region_data_list):
        region = region_data["region"]
        coordinate_start = region.start
        coordinate_end = region.end
        bam_rows = region_data["bam_rows"]

        region_state = RegionBuildState()

        for row_index, bam_row in enumerate(bam_rows):
            _build_bam_row_track(
                bam_row,
                coordinate_start,
                coordinate_end,
                region_state,
                region_idx,
                row_index,
            )

        if not region_state.first_plot_figure:
            logger.debug(f"No segments to plot for region {region_idx + 1}.")
            continue

        region_layout, region_reset_bounds = _build_region_layout_column(region_data, region_state)
        all_region_layouts.append(region_layout)
        all_region_reset_bounds.append(
            (
                region_state.plot_figures,
                region_state.arrow_sources,
                region_state.arrow_renderers,
                region_reset_bounds,
            )
        )
        all_arrow_sources.extend(region_state.arrow_sources)
        all_arrow_renderers.extend(region_state.arrow_renderers)
        all_plot_figures.extend(region_state.plot_figures)

    if not all_region_layouts:
        no_data_div = Div(
            text="No alignments in the requested region(s).",
            styles={"color": "#666", "font-size": "14px", "padding": "20px"},
        )
        all_region_layouts = [column(no_data_div, sizing_mode="stretch_both")]

    final_layout = row(*all_region_layouts, sizing_mode="stretch_both", spacing=20)
    add_multi_region_callbacks(
        all_plot_figures,
        all_arrow_sources,
        all_arrow_renderers,
        all_region_reset_bounds,
    )
    save_plot_with_modal(final_layout, output_file, output_config.prefix)
    return output_file
