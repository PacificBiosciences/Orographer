"""Bokeh figure creation and UI scaffolding (coords, separators, haplotype labels)."""

import math

from bokeh.layouts import column, row
from bokeh.models import (
    BoxZoomTool,
    Button,
    Checkbox,
    ColumnDataSource,
    CustomJS,
    Div,
    PanTool,
    Spacer,
    TapTool,
    TextInput,
    WheelPanTool,
)
from bokeh.plotting import figure

from .callbacks import get_exon_click_callback
from .utils import COORD_INPUT_CENTER_OFFSET_PX, PLOT_CONFIG, load_javascript


def _gene_strand_modal_label(strand_token: str) -> str:
    """Human-readable strand for the exon detail modal (not raw GTF +/-)."""
    if strand_token == "+":
        return "Forward (+)"
    if strand_token == "-":
        return "Reverse (-)"
    return strand_token


def create_bokeh_figure(coordinate_start, coordinate_end, total_height):
    """Create and configure the Bokeh figure with axes and tools."""
    top_padding = max(total_height * 0.02, 0.5)
    plot_figure = figure(
        width=(18 * 80),
        height=int(max(6, total_height * 0.6) * 80),
        x_range=(coordinate_start, coordinate_end),
        y_range=(total_height, -top_padding),
        x_axis_label="Genomic Position",
        y_axis_label="Reads",
        tools="reset",
        toolbar_location="right",
        sizing_mode="stretch_both",
    )
    plot_figure.add_tools(BoxZoomTool())
    plot_figure.add_tools(PanTool(dimensions="height"))
    plot_figure.toolbar.active_drag = plot_figure.select_one(BoxZoomTool)
    wheel_pan = WheelPanTool(dimension="height")
    wheel_pan.visible = False
    plot_figure.add_tools(wheel_pan)
    plot_figure.toolbar.active_scroll = wheel_pan
    tap_tool = TapTool()
    plot_figure.add_tools(tap_tool)
    tap_tool.visible = False
    plot_figure.toolbar.active_tap = tap_tool
    plot_figure.xaxis[0].formatter.use_scientific = False
    plot_figure.xaxis[0].formatter.power_limit_low = 0
    plot_figure.xaxis[0].formatter.power_limit_high = 0
    plot_figure.xgrid.grid_line_color = None
    plot_figure.ygrid.grid_line_color = None
    plot_figure.yaxis.major_tick_line_color = None
    plot_figure.yaxis.minor_tick_line_color = None
    plot_figure.yaxis.major_label_text_color = None
    plot_figure.yaxis.axis_line_color = None
    return plot_figure, tap_tool


def create_bokeh_figure_shared_x(shared_x_range, total_height):
    """Create a figure that shares x_range with another (multi-BAM stacked layout)."""
    top_padding = max(total_height * 0.02, 0.5)
    plot_figure = figure(
        width=(18 * 80),
        height=int(max(6, total_height * 0.6) * 80),
        x_range=shared_x_range,
        y_range=(total_height, -top_padding),
        x_axis_label="Genomic Position",
        y_axis_label="Reads",
        tools="box_zoom,reset",
        toolbar_location=None,
        sizing_mode="stretch_both",
    )
    plot_figure.add_tools(PanTool(dimensions="height"))
    # Shared-x rows have no visible toolbar; force box zoom as active drag so drag
    # interactions keep working on non-top samples (e.g. trio layouts).
    plot_figure.toolbar.active_drag = plot_figure.select_one(BoxZoomTool)
    wheel_pan = WheelPanTool(dimension="height")
    wheel_pan.visible = False
    plot_figure.add_tools(wheel_pan)
    plot_figure.toolbar.active_scroll = wheel_pan
    tap_tool = TapTool()
    plot_figure.add_tools(tap_tool)
    tap_tool.visible = False
    plot_figure.toolbar.active_tap = tap_tool
    plot_figure.xaxis[0].formatter.use_scientific = False
    plot_figure.xaxis[0].formatter.power_limit_low = 0
    plot_figure.xaxis[0].formatter.power_limit_high = 0
    plot_figure.xgrid.grid_line_color = None
    plot_figure.ygrid.grid_line_color = None
    plot_figure.yaxis.major_tick_line_color = None
    plot_figure.yaxis.minor_tick_line_color = None
    plot_figure.yaxis.major_label_text_color = None
    plot_figure.yaxis.axis_line_color = None
    return plot_figure, tap_tool


def create_vcf_track_figure(main_figure):
    """Create a figure for the VCF variant track; links x_range to main figure."""
    vcf_figure = figure(
        width=(18 * 80),
        height=20,
        x_range=main_figure.x_range,
        y_range=(1.0, 0),
        tools="box_zoom,reset",
        toolbar_location=None,
        sizing_mode="stretch_width",
    )
    vcf_figure.xaxis[0].formatter.use_scientific = False
    vcf_figure.xaxis[0].formatter.power_limit_low = 0
    vcf_figure.xaxis[0].formatter.power_limit_high = 0
    vcf_figure.xgrid.grid_line_color = None
    vcf_figure.ygrid.grid_line_color = None
    vcf_figure.yaxis.major_tick_line_color = None
    vcf_figure.yaxis.minor_tick_line_color = None
    vcf_figure.yaxis.major_label_text_color = None
    vcf_figure.yaxis.axis_line_color = None
    return vcf_figure


def create_gene_track_figure(main_figure, gene_track_height):
    """Create a separate figure for the gene track; links x_range to main figure."""
    label_padding = 0.5
    gene_figure = figure(
        width=(18 * 80),
        height=120,
        x_range=main_figure.x_range,
        y_range=(gene_track_height, -label_padding),
        tools="box_zoom,reset",
        toolbar_location=None,
        sizing_mode="stretch_width",
    )
    gene_figure.xaxis[0].formatter.use_scientific = False
    gene_figure.xaxis[0].formatter.power_limit_low = 0
    gene_figure.xaxis[0].formatter.power_limit_high = 0
    gene_figure.xgrid.grid_line_color = None
    gene_figure.ygrid.grid_line_color = None
    gene_figure.yaxis.major_tick_line_color = None
    gene_figure.yaxis.minor_tick_line_color = None
    gene_figure.yaxis.major_label_text_color = None
    gene_figure.yaxis.axis_line_color = None
    gene_figure.yaxis.axis_label = "Genes"
    return gene_figure


def create_genomic_x_axis_strip(shared_x_range):
    """Thin bottom row: sole visible genomic x-axis (linked ``shared_x_range``)."""
    fig = figure(
        width=(18 * 80),
        height=48,
        x_range=shared_x_range,
        y_range=(1, 0),
        toolbar_location=None,
        tools=[],
        sizing_mode="stretch_width",
        outline_line_color=None,
    )
    ax = fig.xaxis[0]
    ax.formatter.use_scientific = False
    ax.formatter.power_limit_low = 0
    ax.formatter.power_limit_high = 0
    fig.xaxis.axis_label = "Genomic Position"
    fig.xgrid.grid_line_color = None
    fig.ygrid.grid_line_color = None
    fig.yaxis.visible = False
    fig.min_border_top = 0
    fig.min_border_bottom = 36
    # One invisible glyph silences W-1000 MISSING_RENDERERS (axis-only figure).
    xm = (shared_x_range.start + shared_x_range.end) / 2
    fig.scatter(
        [xm],
        [0.5],
        size=0,
        alpha=0,
        line_alpha=0,
        fill_alpha=0,
        visible=False,
    )
    return fig


def format_region_size(size_bp):
    """Format region size in a human-readable format (bp, kb, Mb, etc.)."""
    if size_bp < 1000:
        return f"{size_bp} bp"
    elif size_bp < 1000000:
        return f"{size_bp / 1000.0:.1f} kb"
    else:
        return f"{size_bp / 1000000.0:.1f} Mb"


def create_coordinate_display(
    plot_figure,
    chrom,
    coordinate_start,
    coordinate_end,
    one_bp_renderers=None,
    one_bp_markers=None,
    one_bp_texts=None,
    one_bp_segments=None,
    alignment_label_renderers=None,
    default_hide_alignment_numbers=False,
):
    """Create coordinate displays: static full region at top + editable view row.
    If one_bp_renderers is provided, add a checkbox to hide/show 1bp indels.
    one_bp_markers, one_bp_texts, one_bp_segments restore LOD when re-showing.
    If alignment_label_renderers is provided, add checkbox to hide alignment numbers.
    default_hide_alignment_numbers: if True (e.g. paraphase), numbers hidden initially.
    """
    start_str = f"{coordinate_start:,}"
    end_str = f"{coordinate_end:,}"
    region_size = coordinate_end - coordinate_start + 1
    size_str = format_region_size(region_size)

    _region_label_styles = {
        "font-size": "14px",
        "font-family": "Arial, sans-serif",
        "color": PLOT_CONFIG["sample_label_color"],
        "padding-right": "20px",
        "text-align": "right",
        "margin": "0",
        "padding-top": "0",
        "padding-bottom": "2px",
    }
    region_range_str = f"Original region ({size_str}):"
    original_region_block_wide = column(
        Div(text=region_range_str, styles=_region_label_styles),
        Div(text=f"{chrom}:{start_str}-{end_str}", styles=_region_label_styles),
        spacing=0,
        name="orographer_original_region_wide",
    )
    original_region_block_narrow = column(
        Div(text=region_range_str, styles=_region_label_styles),
        Div(text=f"{chrom}:{start_str}-{end_str}", styles=_region_label_styles),
        spacing=0,
    )

    coord_input = TextInput(
        value=f"{chrom}:{start_str}-{end_str}",
        title="",
        width=280,
        height=33,
        styles={
            "font-size": "16px",
            "font-family": "Arial, sans-serif",
            "font-weight": "bold",
        },
    )
    go_button = Button(label="Go", button_type="primary", width=50, height=33)
    error_div = Div(
        text="",
        width=260,
        styles={
            "font-size": "12px",
            "font-family": "Arial, sans-serif",
            "color": "#ff0000",
            "text-align": "left",
            "padding-left": "8px",
            "min-height": "0",
            "height": "0",
            "overflow": "hidden",
            "margin": "0",
            "padding-top": "0",
            "padding-bottom": "0",
        },
    )
    orig_start = coordinate_start
    orig_end = coordinate_end

    go_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "coord_input": coord_input,
            "error_div": error_div,
            "orig_start": orig_start,
            "orig_end": orig_end,
        },
        code=load_javascript("coord_go_callback.js"),
    )
    go_button.js_on_click(go_callback)
    coord_input.js_on_change("value", go_callback)

    view_size_div = Div(
        text=size_str,
        width=64,
        css_classes=["orographer-view-size"],
        styles={
            "font-size": "12px",
            "font-family": "Arial, sans-serif",
            "color": PLOT_CONFIG["sample_label_color"],
            "text-align": "center",
            "margin": "0",
            "padding": "2px 0 0 0",
            "line-height": "1",
        },
    )
    view_callback = CustomJS(
        args={
            "coord_input": coord_input,
            "error_div": error_div,
            "chrom": chrom,
            "view_size_div": view_size_div,
        },
        code=load_javascript("view_callback.js"),
        module=True,
    )
    plot_figure.x_range.js_on_change("start", view_callback)
    plot_figure.x_range.js_on_change("end", view_callback)
    plot_figure.x_range.js_on_change("change", view_callback)

    zoom_out_btn = Button(
        label="\u2212",
        button_type="default",
        width=28,
        height=22,
        margin=(0, 1, 0, 0),
        styles={
            "font-size": "16px",
            "font-weight": "bold",
            "line-height": "1",
            "padding": "0",
            "text-align": "center",
        },
    )
    zoom_in_btn = Button(
        label="+",
        button_type="default",
        width=28,
        height=22,
        margin=(0, 0, 0, 1),
        styles={
            "font-size": "16px",
            "font-weight": "bold",
            "line-height": "1",
            "padding": "0",
            "text-align": "center",
        },
    )
    zoom_in_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "factor": 0.5,
            "orig_start": orig_start,
            "orig_end": orig_end,
        },
        code=load_javascript("zoom_buttons_callback.js"),
    )
    zoom_out_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "factor": 2,
            "orig_start": orig_start,
            "orig_end": orig_end,
        },
        code=load_javascript("zoom_buttons_callback.js"),
    )
    zoom_in_btn.js_on_click(zoom_in_callback)
    zoom_out_btn.js_on_click(zoom_out_callback)
    view_size_with_zoom = column(
        view_size_div,
        row(zoom_out_btn, zoom_in_btn, spacing=0),
        spacing=0,
        align="start",
    )

    hide_1bp_checkbox = None
    if one_bp_renderers is not None and len(one_bp_renderers) > 0:
        hide_1bp_checkbox = Checkbox(
            label="Hide 1bp INDELs", active=False, width=140, margin=(5, 5, 0, 5)
        )
        hide_1bp_callback = CustomJS(
            args={
                "one_bp_renderers": one_bp_renderers,
                "one_bp_markers": one_bp_markers or [],
                "one_bp_texts": one_bp_texts or [],
                "one_bp_segments": one_bp_segments or [],
                "x_range": plot_figure.x_range,
            },
            code=load_javascript("hide_1bp_callback.js"),
        )
        hide_1bp_checkbox.js_on_change("active", hide_1bp_callback)

    hide_alignment_numbers_checkbox = None
    if alignment_label_renderers is not None and len(alignment_label_renderers) > 0:
        hide_alignment_numbers_checkbox = Checkbox(
            label="Hide alignment numbers",
            active=default_hide_alignment_numbers,
            width=180,
            margin=(0, 5, 5, 5),
        )
        if default_hide_alignment_numbers:
            for r in alignment_label_renderers:
                r.visible = False
        hide_alignment_numbers_callback = CustomJS(
            args={"alignment_label_renderers": alignment_label_renderers},
            code=load_javascript("hide_alignment_numbers_callback.js"),
        )
        hide_alignment_numbers_checkbox.js_on_change("active", hide_alignment_numbers_callback)

    checkbox_column_items = []
    if hide_1bp_checkbox is not None:
        checkbox_column_items.append(hide_1bp_checkbox)
    if hide_alignment_numbers_checkbox is not None:
        checkbox_column_items.append(hide_alignment_numbers_checkbox)
    # Center the coordinate input: offset = right_column_width - left_column_width.
    # Tune COORD_INPUT_CENTER_OFFSET_PX in utils.py: decrease to shift input left.
    input_center_offset = COORD_INPUT_CENTER_OFFSET_PX
    center_offset_spacer = Spacer(
        width=input_center_offset,
        height=1,
        sizing_mode="fixed",
        name="orographer_center_offset_spacer",
    )
    center_group = row(
        center_offset_spacer,
        Spacer(sizing_mode="stretch_width"),
        go_button,
        coord_input,
        view_size_with_zoom,
        Spacer(sizing_mode="stretch_width"),
        sizing_mode="stretch_width",
        align="start",
        spacing=0,
    )
    error_row = row(
        Spacer(width=input_center_offset, height=0, sizing_mode="fixed"),
        Spacer(sizing_mode="stretch_width"),
        error_div,
        Spacer(sizing_mode="stretch_width"),
        sizing_mode="stretch_width",
        align="start",
        spacing=0,
    )
    coord_row_items = []
    if checkbox_column_items:
        coord_row_items.append(column(*checkbox_column_items, spacing=0))
    coord_row_items.append(center_group)
    coord_row_items.append(original_region_block_wide)
    coord_row_wide = row(
        *coord_row_items,
        sizing_mode="stretch_width",
        align="start",
    )
    original_region_row_narrow = row(
        Spacer(sizing_mode="stretch_width"),
        original_region_block_narrow,
        sizing_mode="stretch_width",
        align="end",
        visible=False,
        name="orographer_original_region_narrow_row",
    )
    coord_controls = column(
        coord_row_wide,
        error_row,
        original_region_row_narrow,
        sizing_mode="stretch_width",
        align="start",
    )
    return (
        coord_controls,
        hide_1bp_checkbox,
    )


def add_separator_lines(
    plot_figure,
    read_names,
    read_to_y_bottom,
    read_heights,
    coordinate_start,
    coordinate_end,
):
    """Add horizontal dotted lines to separate reads (single multi_line glyph)."""
    xs = []
    ys = []
    xs.append([coordinate_start, coordinate_end])
    ys.append([0, 0])
    for read_name in read_names:
        y_bottom = read_to_y_bottom[read_name] + read_heights[read_name]
        xs.append([coordinate_start, coordinate_end])
        ys.append([y_bottom, y_bottom])
    plot_figure.multi_line(
        xs=xs,
        ys=ys,
        line_color="grey",
        line_width=0.5,
        line_dash="dotted",
        line_alpha=0.3,
        level="underlay",
    )


def get_haplotype_label(haplotype):
    """Generate a label for a haplotype value."""
    if haplotype == 0:
        return "Unassigned"
    return f"{haplotype}"


def add_haplotype_labels(
    plot_figure, group_boundaries, haplotype_order, coordinate_start, coordinate_end
):
    """Add text labels and separator lines for each haplotype group."""
    label_x = []
    label_y = []
    label_texts = []
    separator_xs = []
    separator_ys = []
    label_x_pos = coordinate_start + (coordinate_end - coordinate_start) * 0.01

    for i, haplotype in enumerate(haplotype_order):
        if haplotype not in group_boundaries:
            continue
        y_start, y_end = group_boundaries[haplotype]
        if y_start == y_end:
            continue
        y_center = (y_start + y_end) / 2
        label_x.append(label_x_pos)
        label_y.append(y_center)
        label_texts.append(get_haplotype_label(haplotype))

        next_group_idx = i + 1
        while (
            next_group_idx < len(haplotype_order)
            and haplotype_order[next_group_idx] not in group_boundaries
        ):
            next_group_idx += 1
        if next_group_idx < len(haplotype_order):
            next_haplotype = haplotype_order[next_group_idx]
            if next_haplotype in group_boundaries:
                next_y_start = group_boundaries[next_haplotype][0]
                separator_y = (y_end + next_y_start) / 2
                separator_xs.append([coordinate_start, coordinate_end])
                separator_ys.append([separator_y, separator_y])

    if separator_xs:
        plot_figure.multi_line(
            xs=separator_xs,
            ys=separator_ys,
            line_color=PLOT_CONFIG["sample_label_color"],
            line_width=1.0,
            line_alpha=1,
            level="underlay",
        )
    if label_x:
        label_source = ColumnDataSource(data={"x": label_x, "y": label_y, "text": label_texts})
        plot_figure.text(
            x="x",
            y="y",
            text="text",
            source=label_source,
            text_font_size="11pt",
            text_color="black",
            text_font_style="bold",
            text_align="left",
            text_baseline="middle",
            text_alpha=10,
            text_outline_color="black",
        )


def add_gene_track(
    plot_figure, gene_annotations, gene_track_y_start, coordinate_start, coordinate_end
):
    """Add gene annotation track (IGV-style). Returns height of the gene track."""
    if not gene_annotations:
        return 0

    gene_color = "#3366CC"
    exon_height = 0.4
    gene_spacing = 0.8
    arrow_spacing = 500

    exon_lefts, exon_rights, exon_tops, exon_bottoms = [], [], [], []
    exon_gene_names, exon_strands, exon_numbers = [], [], []
    intron_xs, intron_ys = [], []
    arrow_xs, arrow_ys, arrow_angles = [], [], []
    label_xs, label_ys, label_texts = [], [], []

    gene_rows = []
    gene_row_assignments = {}
    for gene in gene_annotations:
        row_idx = 0
        for existing_row in gene_rows:
            can_use = True
            for end_x, assigned_row in existing_row:
                if gene.start < end_x and assigned_row == row_idx:
                    can_use = False
                    break
            if can_use:
                break
            row_idx += 1
        while len(gene_rows) <= row_idx:
            gene_rows.append([])
        gene_rows[row_idx].append((gene.end + 100, row_idx))
        gene_row_assignments[gene.gene_id] = row_idx

    num_rows = max(len(gene_rows), 1)
    row_height = exon_height + gene_spacing

    for gene in gene_annotations:
        row_idx = gene_row_assignments.get(gene.gene_id, 0)
        gene_y_center = gene_track_y_start + (row_idx + 0.5) * row_height
        sorted_exons = sorted(gene.exons, key=lambda exon: exon[0])

        for exon_start, exon_end, exon_number in sorted_exons:
            vis_start = max(exon_start, coordinate_start)
            vis_end = min(exon_end, coordinate_end)
            if vis_start < vis_end:
                exon_lefts.append(vis_start)
                exon_rights.append(vis_end)
                exon_tops.append(gene_y_center - exon_height / 2)
                exon_bottoms.append(gene_y_center + exon_height / 2)
                exon_gene_names.append(gene.gene_name)
                strand_token = gene.strand or ""
                if strand_token in ("", "."):
                    exon_strands.append("unknown")
                else:
                    exon_strands.append(_gene_strand_modal_label(strand_token))
                exon_numbers.append(str(exon_number))

        for i in range(len(sorted_exons) - 1):
            intron_start = sorted_exons[i][1]
            intron_end = sorted_exons[i + 1][0]
            if intron_end > coordinate_start and intron_start < coordinate_end:
                vis_start = max(intron_start, coordinate_start)
                vis_end = min(intron_end, coordinate_end)
                if vis_start < vis_end:
                    intron_xs.append([vis_start, vis_end])
                    intron_ys.append([gene_y_center, gene_y_center])
                    intron_length = vis_end - vis_start
                    num_arrows = max(1, int(intron_length / arrow_spacing))
                    for j in range(num_arrows):
                        arrow_x = vis_start + (j + 0.5) * (intron_length / num_arrows)
                        arrow_xs.append(arrow_x)
                        arrow_ys.append(gene_y_center)
                        arrow_angles.append(-math.pi / 2 if gene.strand == "+" else math.pi / 2)

        if len(sorted_exons) == 1:
            exon_start, exon_end, _ = sorted_exons[0]
            vis_start = max(exon_start, coordinate_start)
            vis_end = min(exon_end, coordinate_end)
            exon_length = vis_end - vis_start
            if exon_length > 100:
                num_arrows = max(1, int(exon_length / arrow_spacing))
                for j in range(num_arrows):
                    arrow_x = vis_start + (j + 0.5) * (exon_length / num_arrows)
                    arrow_xs.append(arrow_x)
                    arrow_ys.append(gene_y_center)
                    arrow_angles.append(-math.pi / 2 if gene.strand == "+" else math.pi / 2)

        gene_vis_start = max(gene.start, coordinate_start)
        gene_vis_end = min(gene.end, coordinate_end)
        label_xs.append((gene_vis_start + gene_vis_end) / 2)
        label_ys.append(gene_y_center - exon_height / 2)
        label_texts.append(gene.gene_name)

    if intron_xs:
        plot_figure.multi_line(
            xs=intron_xs,
            ys=intron_ys,
            line_color=gene_color,
            line_width=1.5,
            line_alpha=0.8,
        )

    if exon_lefts:
        exon_source = ColumnDataSource(
            data={
                "left": exon_lefts,
                "right": exon_rights,
                "top": exon_tops,
                "bottom": exon_bottoms,
                "gene_name": exon_gene_names,
                "gene_strand": exon_strands,
                "exon_number": exon_numbers,
            }
        )
        exon_renderer = plot_figure.quad(
            left="left",
            right="right",
            top="top",
            bottom="bottom",
            source=exon_source,
            fill_color=gene_color,
            line_color=gene_color,
            fill_alpha=0.8,
            line_alpha=1.0,
        )
        tap_tool = TapTool()
        tap_tool.renderers = [exon_renderer]
        plot_figure.add_tools(tap_tool)
        plot_figure.toolbar.active_tap = tap_tool
        exon_click_callback = get_exon_click_callback(exon_source)
        exon_source.selected.js_on_change("indices", exon_click_callback)

    if arrow_xs:
        arrow_source = ColumnDataSource(data={"x": arrow_xs, "y": arrow_ys, "angle": arrow_angles})
        plot_figure.scatter(
            x="x",
            y="y",
            source=arrow_source,
            marker="triangle",
            size=8,
            angle="angle",
            fill_color="white",
            line_color=gene_color,
            fill_alpha=0.9,
            line_alpha=1.0,
        )

    if label_xs:
        label_source = ColumnDataSource(data={"x": label_xs, "y": label_ys, "text": label_texts})
        plot_figure.text(
            x="x",
            y="y",
            text="text",
            source=label_source,
            text_font_size="8pt",
            text_color=gene_color,
            text_font_style="italic",
            text_align="center",
            text_baseline="bottom",
            text_alpha=0.9,
        )

    return num_rows * row_height + 1.0
