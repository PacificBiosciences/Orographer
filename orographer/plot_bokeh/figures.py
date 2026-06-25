"""Bokeh figure creation and UI scaffolding (coords, separators, haplotype labels)."""

import base64
import math
from dataclasses import dataclass
from io import BytesIO

import numpy as np
from bokeh.events import MouseLeave, MouseMove
from bokeh.layouts import column, row
from bokeh.models import (
    BasicTicker,
    BoxZoomTool,
    Button,
    Checkbox,
    ColumnDataSource,
    CustomJS,
    CustomJSHover,
    CustomJSTicker,
    CustomJSTickFormatter,
    Div,
    HoverTool,
    PanTool,
    Range1d,
    Spacer,
    Span,
    TapTool,
    TextInput,
    WheelPanTool,
)
from bokeh.plotting import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter

from .callbacks import get_exon_click_callback
from .data import READ_FILTER_LAYOUT_MODES, add_read_filter_visibility_columns
from .utils import (
    COORD_INPUT_CENTER_OFFSET_PX,
    PLOT_CONFIG,
    build_haplotype_color_map,
    get_haplotype_color,
    load_javascript,
)

GENE_DETAIL_SPAN_BP = 200_000
GENE_MEDIUM_SPAN_BP = 1_000_000
GENE_LABEL_SCREEN_PX_PER_CHAR = 7
GENE_LABEL_PADDING_SCREEN_PX = 16
GENE_INTRON_ARROW_SPACING_BP = 500
GENE_MAX_ARROWS_PER_FEATURE = 12
GENE_MIN_EXON_ARROW_BP = 100
GENE_OVERVIEW_BODY_HEIGHT = 0.12
GENE_EXON_HEIGHT = 0.4
GENE_ROW_SPACING = 0.8
GENE_TRACK_PLOT_WIDTH_PX = 18 * 80
GENE_DETAIL_EXON_ALPHA = 0.8
GENE_MEDIUM_EXON_ALPHA = 0.55
GENE_HIDDEN_ALPHA = 0.0
GENE_OVERVIEW_BODY_ALPHA = 0.75
GENE_LABEL_ALPHA = 0.9


def enable_button_click_hold(
    button: Button,
    *,
    initial_delay: int = 260,
    repeat_delay: int = 95,
) -> None:
    """Tag a button for pointer-hold repeats installed by the page bootstrap script."""
    tags = list(button.tags or [])
    hold_tags = [
        "orographer-click-hold",
        f"orographer-click-hold-delay:{initial_delay}",
        f"orographer-click-hold-repeat:{repeat_delay}",
    ]
    button.tags = [*tags, *[tag for tag in hold_tags if tag not in tags]]


@dataclass(frozen=True)
class GeneTrackRecord:
    """Collapsed gene annotation used by the scalable gene track."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: tuple[tuple[int, int, int], ...]
    representative_transcript_id: str | None = None
    representative_transcript_name: str | None = None
    representative_selection_method: str = "gene annotation"


def _layout_mode_column(column_name: str, mode: str) -> str:
    """Return the alternate source-column name for a read filter layout mode."""
    return f"{column_name}_{mode}"


def _gene_strand_modal_label(strand_token: str) -> str:
    """Human-readable strand for the exon detail modal (not raw GTF +/-)."""
    if strand_token == "+":
        return "Forward (+)"
    if strand_token == "-":
        return "Reverse (-)"
    return strand_token


def _configure_genomic_x_axis(plot_figure) -> None:
    """Keep only labeled major tick marks on genomic x-axes."""
    axis = plot_figure.xaxis[0]
    axis.ticker = CustomJSTicker(
        major_code=load_javascript("genomic_x_major_ticker.js"),
        minor_code="return [];",
    )
    axis.formatter.use_scientific = False
    axis.formatter.power_limit_low = 0
    axis.formatter.power_limit_high = 0
    axis.minor_tick_line_color = None


def _annotation_interval(start: int, end: int) -> tuple[int, int]:
    """Return display interval boundaries for a 1-based inclusive annotation feature."""
    return start, end + 1


def _format_genomic_coordinates(
    chrom: str,
    start: int,
    end: int,
    *,
    include_size: bool = False,
) -> str:
    """Return a compact coordinate label for annotation metadata."""
    interval = f"{start:,}-{end:,}"
    coordinates = interval
    if chrom:
        coordinates = f"{chrom}:{interval}"
    if include_size:
        size_bp = end - start + 1
        return f"{coordinates} ({size_bp:,} bp)"
    return coordinates


def _initial_gene_track_mode(coordinate_start: int, coordinate_end: int) -> str:
    span_bp = max(1, coordinate_end - coordinate_start + 1)
    if span_bp <= GENE_DETAIL_SPAN_BP:
        return "detail"
    if span_bp <= GENE_MEDIUM_SPAN_BP:
        return "medium"
    return "overview"


def _mode_values(
    mode: str,
    *,
    overview: float,
    medium: float,
    detail: float,
) -> tuple[float, float, float, float]:
    active = {"overview": overview, "medium": medium, "detail": detail}
    return active[mode], overview, medium, detail


def _append_gene_arrows(
    arrow_data: dict[str, list],
    *,
    start: float,
    end: float,
    y_center: float,
    strand: str,
    mode: str,
) -> None:
    feature_length = end - start
    if feature_length <= 0:
        return
    arrow_count = max(1, int(feature_length / GENE_INTRON_ARROW_SPACING_BP))
    arrow_count = min(GENE_MAX_ARROWS_PER_FEATURE, arrow_count)
    arrow_alpha = _mode_values(
        mode,
        overview=GENE_HIDDEN_ALPHA,
        medium=GENE_HIDDEN_ALPHA,
        detail=0.9,
    )
    angle = -math.pi / 2 if strand == "+" else math.pi / 2
    for index in range(arrow_count):
        arrow_x = start + (index + 0.5) * (feature_length / arrow_count)
        arrow_data["x"].append(arrow_x)
        arrow_data["y"].append(y_center)
        arrow_data["angle"].append(angle)
        arrow_data["fill_alpha"].append(arrow_alpha[0])
        arrow_data["fill_alpha_overview"].append(arrow_alpha[1])
        arrow_data["fill_alpha_medium"].append(arrow_alpha[2])
        arrow_data["fill_alpha_detail"].append(arrow_alpha[3])
        arrow_data["line_alpha"].append(1.0 if arrow_alpha[0] else 0.0)
        arrow_data["line_alpha_overview"].append(0.0)
        arrow_data["line_alpha_medium"].append(0.0)
        arrow_data["line_alpha_detail"].append(1.0)


def _append_gene_intron(
    intron_data: dict[str, list],
    *,
    start: float,
    end: float,
    y_center: float,
    mode: str,
) -> None:
    if start >= end:
        return
    line_alpha = _mode_values(
        mode,
        overview=GENE_HIDDEN_ALPHA,
        medium=0.6,
        detail=0.8,
    )
    intron_data["xs"].append([start, end])
    intron_data["ys"].append([y_center, y_center])
    intron_data["line_alpha"].append(line_alpha[0])
    intron_data["line_alpha_overview"].append(line_alpha[1])
    intron_data["line_alpha_medium"].append(line_alpha[2])
    intron_data["line_alpha_detail"].append(line_alpha[3])


def _deduplicate_gene_annotations(gene_annotations) -> list[GeneTrackRecord]:
    records_by_id = {}
    for gene in gene_annotations:
        gene_id = str(getattr(gene, "gene_id", "") or getattr(gene, "gene_name", ""))
        gene_name = str(getattr(gene, "gene_name", "") or gene_id)
        strand = str(getattr(gene, "strand", "") or "")
        transcript_id = getattr(gene, "representative_transcript_id", None)
        transcript_name = getattr(gene, "representative_transcript_name", None)
        selection_method = str(
            getattr(gene, "representative_selection_method", "gene annotation") or "gene annotation"
        )
        chrom = str(getattr(gene, "chrom", "") or "")
        start = int(gene.start)
        end = int(gene.end)
        exons = tuple(
            sorted(
                (
                    (int(exon_start), int(exon_end), int(exon_number))
                    for exon_start, exon_end, exon_number in getattr(gene, "exons", [])
                ),
                key=lambda exon: (exon[0], exon[1], exon[2]),
            )
        )
        if gene_id not in records_by_id:
            records_by_id[gene_id] = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "exons": exons,
                "representative_transcript_id": transcript_id,
                "representative_transcript_name": transcript_name,
                "representative_selection_method": selection_method,
            }
            continue

    records = [
        GeneTrackRecord(
            gene_id=record["gene_id"],
            gene_name=record["gene_name"],
            chrom=record["chrom"],
            start=record["start"],
            end=record["end"],
            strand=record["strand"],
            exons=record["exons"],
            representative_transcript_id=record["representative_transcript_id"],
            representative_transcript_name=record["representative_transcript_name"],
            representative_selection_method=record["representative_selection_method"],
        )
        for record in records_by_id.values()
    ]
    return sorted(records, key=lambda record: (record.start, record.end, record.gene_name))


def _gene_label_width_bp(
    gene_name: str,
    coordinate_start: int,
    coordinate_end: int,
    plot_width_px: int,
) -> float:
    span_bp = max(1, coordinate_end - coordinate_start + 1)
    label_width_px = len(gene_name) * GENE_LABEL_SCREEN_PX_PER_CHAR
    label_width_px += GENE_LABEL_PADDING_SCREEN_PX
    return span_bp * label_width_px / max(1, plot_width_px)


def _gene_label_interval(
    gene: GeneTrackRecord,
    coordinate_start: int,
    coordinate_end: int,
    plot_width_px: int,
) -> tuple[float, float]:
    gene_left, gene_right = _annotation_interval(gene.start, gene.end)
    visible_start = max(gene_left, coordinate_start)
    visible_end = min(gene_right, coordinate_end + 1)
    label_center = (visible_start + visible_end) / 2
    half_label_width = _gene_label_width_bp(
        gene.gene_name,
        coordinate_start,
        coordinate_end,
        plot_width_px,
    )
    half_label_width /= 2
    return label_center - half_label_width, label_center + half_label_width


def _assign_gene_rows(
    genes: list[GeneTrackRecord],
    coordinate_start: int,
    coordinate_end: int,
    plot_width_px: int,
    *,
    reserve_label_width: bool,
) -> tuple[dict[str, int], int]:
    row_ends: list[float] = []
    row_assignments = {}
    span_bp = max(1, coordinate_end - coordinate_start + 1)
    body_padding_bp = span_bp * GENE_LABEL_PADDING_SCREEN_PX / max(1, plot_width_px)
    for gene in genes:
        gene_left, gene_right = _annotation_interval(gene.start, gene.end)
        if reserve_label_width:
            label_left, label_right = _gene_label_interval(
                gene,
                coordinate_start,
                coordinate_end,
                plot_width_px,
            )
        else:
            label_left = gene_left - body_padding_bp
            label_right = gene_right + body_padding_bp
        interval_start = min(gene_left, label_left)
        interval_end = max(gene_right, label_right)
        row_idx = 0
        for idx, row_end in enumerate(row_ends):
            if interval_start > row_end:
                row_idx = idx
                break
        else:
            row_idx = len(row_ends)
            row_ends.append(float("-inf"))
        row_ends[row_idx] = interval_end
        row_assignments[gene.gene_id] = row_idx
    return row_assignments, max(1, len(row_ends))


def _add_gene_track_zoom_callback(
    plot_figure,
    source_specs: list[tuple[ColumnDataSource, str]],
) -> None:
    callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "sources": [source for source, _target_column in source_specs],
            "target_columns": [target_column for _source, target_column in source_specs],
            "detail_span_bp": GENE_DETAIL_SPAN_BP,
            "medium_span_bp": GENE_MEDIUM_SPAN_BP,
        },
        code=load_javascript("gene_track_zoom_callback.js"),
    )
    plot_figure.x_range.js_on_change("start", callback)
    plot_figure.x_range.js_on_change("end", callback)


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
    box_zoom_tool = BoxZoomTool()
    plot_figure.add_tools(box_zoom_tool)
    plot_figure.add_tools(PanTool(dimensions="height"))
    plot_figure.toolbar.active_drag = box_zoom_tool
    wheel_pan = WheelPanTool(dimension="height")
    wheel_pan.visible = False
    plot_figure.add_tools(wheel_pan)
    plot_figure.toolbar.active_scroll = wheel_pan
    tap_tool = TapTool()
    plot_figure.add_tools(tap_tool)
    tap_tool.visible = False
    plot_figure.toolbar.active_tap = tap_tool
    _configure_genomic_x_axis(plot_figure)
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
        tools="reset",
        toolbar_location=None,
        sizing_mode="stretch_both",
    )
    box_zoom_tool = BoxZoomTool()
    plot_figure.add_tools(box_zoom_tool)
    plot_figure.add_tools(PanTool(dimensions="height"))
    # Shared-x rows have no visible toolbar; force box zoom as active drag so drag
    # interactions keep working on non-top samples (e.g. trio layouts).
    plot_figure.toolbar.active_drag = box_zoom_tool
    wheel_pan = WheelPanTool(dimension="height")
    wheel_pan.visible = False
    plot_figure.add_tools(wheel_pan)
    plot_figure.toolbar.active_scroll = wheel_pan
    tap_tool = TapTool()
    plot_figure.add_tools(tap_tool)
    tap_tool.visible = False
    plot_figure.toolbar.active_tap = tap_tool
    _configure_genomic_x_axis(plot_figure)
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
    _configure_genomic_x_axis(vcf_figure)
    vcf_figure.xgrid.grid_line_color = None
    vcf_figure.ygrid.grid_line_color = None
    vcf_figure.yaxis.major_tick_line_color = None
    vcf_figure.yaxis.minor_tick_line_color = None
    vcf_figure.yaxis.major_label_text_color = None
    vcf_figure.yaxis.axis_line_color = None
    return vcf_figure


def create_coverage_figure(main_figure):
    """Create a coverage depth track linked to *main_figure*'s x_range."""
    cov_figure = figure(
        width=(18 * 80),
        height=60,
        x_range=main_figure.x_range,
        y_range=(0, 1),
        tools="reset",
        toolbar_location=None,
        sizing_mode="stretch_width",
    )
    _configure_genomic_x_axis(cov_figure)
    cov_figure.xaxis.major_tick_line_color = None
    cov_figure.xaxis.major_label_text_color = None
    cov_figure.xgrid.grid_line_color = None
    cov_figure.ygrid.grid_line_color = "#eeeeee"
    cov_figure.yaxis.minor_tick_line_color = None
    cov_figure.yaxis.axis_label = "Depth"
    cov_figure.yaxis.axis_label_text_font_size = "10px"
    cov_figure.yaxis[0].ticker = BasicTicker(desired_num_ticks=4)
    cov_figure.yaxis[0].formatter = CustomJSTickFormatter(
        code="""
            if (!Number.isFinite(tick)) return "";
            const absTick = Math.abs(tick);
            if (Math.sign(absTick - 1000) === -1) {
                return String(Math.round(tick));
            }
            const exponent = Math.floor(Math.log10(absTick));
            const mantissa = tick / Math.pow(10, exponent);
            let label = mantissa.toFixed(1);
            if (label.endsWith(".0")) {
                label = label.slice(0, -2);
            }
            return label + "e" + String(exponent);
        """,
    )
    return cov_figure


def add_coverage_hover_to_figure(cov_figure, source, hp_source=None):
    """Add a dedicated invisible coverage line renderer for reliable hover hit testing."""
    hover_renderer = cov_figure.line(
        x="x",
        y="y",
        source=source,
        line_alpha=0,
        line_width=8,
    )
    cov_figure.add_tools(
        HoverTool(
            renderers=[hover_renderer],
            tooltips=[("Position", "$x{0,0}"), ("Depth", "@y{custom}")],
            formatters={
                "@y": CustomJSHover(
                    args={"source": source, "hp_source": hp_source},
                    code=load_javascript("coverage_depth_hover_callback.js"),
                )
            },
            mode="vline",
        )
    )
    return hover_renderer


def update_coverage_y_range(cov_figure, y_values):
    """Set a finite y-range for coverage tracks, including empty tracks."""
    max_depth = max([0, *list(y_values)])
    cov_figure.y_range.start = 0
    cov_figure.y_range.end = max(1, max_depth * 1.05)


def _coverage_haplotype_sort_key(haplotype):
    """Return a stable ordering for numeric and named haplotype coverage tracks."""
    if isinstance(haplotype, int):
        return (0, haplotype, "")
    haplotype_text = str(haplotype)
    try:
        return (0, int(haplotype_text), "")
    except ValueError:
        return (1, 0, haplotype_text.lower())


def _coverage_haplotype_label(haplotype) -> str:
    """Return a compact display label for a coverage haplotype key."""
    haplotype_text = str(haplotype)
    if haplotype == 0 or haplotype_text.lower() in {"0", "unassigned", ""}:
        return "Unassigned"
    try:
        return f"HP{int(haplotype_text)}"
    except ValueError:
        return haplotype_text


def add_cursor_guide_to_figures(plot_figures):
    """Add synchronized cursor-guide spans to figures and return the span models."""
    cursor_spans = []
    for plot_figure in plot_figures:
        cursor_span = Span(
            dimension="height",
            location=0,
            line_color=PLOT_CONFIG["cursor_guide_line_color"],
            line_alpha=PLOT_CONFIG["cursor_guide_line_alpha"],
            line_width=PLOT_CONFIG["cursor_guide_line_width"],
            level="overlay",
            visible=False,
        )
        plot_figure.add_layout(cursor_span)
        cursor_spans.append(cursor_span)
    return cursor_spans


def add_cursor_guide_callbacks(plot_figures, cursor_spans, cursor_guide_checkbox) -> None:
    """Synchronize cursor-guide spans across every plotted track."""
    if not plot_figures or not cursor_spans or cursor_guide_checkbox is None:
        return

    move_callback = CustomJS(
        args={
            "cursor_guide_spans": cursor_spans,
            "cursor_guide_checkbox": cursor_guide_checkbox,
        },
        code=load_javascript("cursor_guide_mousemove_callback.js"),
    )
    leave_callback = CustomJS(
        args={"cursor_guide_spans": cursor_spans},
        code=load_javascript("cursor_guide_mouseleave_callback.js"),
    )
    toggle_callback = CustomJS(
        args={"cursor_guide_spans": cursor_spans},
        code=load_javascript("cursor_guide_toggle_callback.js"),
    )
    cursor_guide_checkbox.js_on_change("active", toggle_callback)
    for plot_figure in plot_figures:
        plot_figure.js_on_event(MouseMove, move_callback)
        plot_figure.js_on_event(MouseLeave, leave_callback)


def build_coverage_source_data(
    coverage_tracks: dict,
    haplotype_color_map: dict | None = None,
) -> dict[str, dict | None]:
    """Compute raw data dicts for coverage sources without creating Bokeh models.

    Returns {"total": dict|None, "hp": dict|None} where each dict is ready to
    pass directly to ColumnDataSource(data=...).
    """
    if not coverage_tracks:
        return {"total": None, "hp": None}

    total = coverage_tracks.get(-1)
    total_data = None
    if total and any(y > 0 for y in total[1]):
        x_vals, y_vals = total
        total_data = {"x": list(x_vals), "y": list(y_vals)}

    hp_xs, hp_ys, hp_colors, hp_names = [], [], [], []
    hp_label_xs, hp_label_ys = [], []
    coverage_items = sorted(
        coverage_tracks.items(),
        key=lambda item: _coverage_haplotype_sort_key(item[0]),
    )
    coverage_haplotypes = [hp for hp, _track in coverage_items if hp != -1]
    color_map = haplotype_color_map or build_haplotype_color_map(coverage_haplotypes)
    for hp, (x_vals, y_vals) in coverage_items:
        if hp == -1:
            continue
        hp_xs.append(list(x_vals))
        hp_ys.append(list(y_vals))
        hp_colors.append(color_map[hp])
        hp_names.append(_coverage_haplotype_label(hp))
        peak_idx = max(range(len(y_vals)), key=lambda i: y_vals[i])
        hp_label_xs.append(x_vals[peak_idx])
        hp_label_ys.append(y_vals[peak_idx])

    hp_data = None
    if hp_xs:
        hp_data = {
            "xs": hp_xs,
            "ys": hp_ys,
            "colors": hp_colors,
            "names": hp_names,
            "label_x": hp_label_xs,
            "label_y": hp_label_ys,
        }
    return {"total": total_data, "hp": hp_data}


def add_coverage_to_figure(
    cov_figure, coverage_tracks, haplotype_color_map=None
) -> dict[str, object]:
    """Render coverage depth series into *cov_figure*.

    coverage_tracks: dict mapping haplotype value to (x_positions, y_depths).
    Key -1 is the total series; other keys are HP tag values.

    Returns {"total": ColumnDataSource|None, "hp": ColumnDataSource|None}.
    """
    if not coverage_tracks:
        return {"total": None, "hp": None}

    cov_data = build_coverage_source_data(coverage_tracks, haplotype_color_map)

    total_source = None
    if cov_data["total"] is not None:
        update_coverage_y_range(cov_figure, cov_data["total"]["y"])
        total_source = ColumnDataSource(cov_data["total"])
        cov_figure.varea(
            x="x",
            y1=0,
            y2="y",
            source=total_source,
            color="#cccccc",
            alpha=0.7,
        )

    hp_source = None
    if cov_data["hp"] is not None:
        hp_source = ColumnDataSource(cov_data["hp"])
        cov_figure.multi_line(
            xs="xs",
            ys="ys",
            line_color="colors",
            line_width=2,
            line_alpha=0.9,
            selection_line_width=3,
            selection_line_alpha=1.0,
            nonselection_line_alpha=0.15,
            source=hp_source,
        )
        cov_figure.text(
            x="label_x",
            y="label_y",
            text="names",
            source=hp_source,
            text_color="colors",
            text_font_size="9pt",
            text_font_style="bold",
            text_align="center",
            text_baseline="bottom",
            text_alpha=0,
            selection_text_alpha=1.0,
            nonselection_text_alpha=0,
            text_outline_color="white",
        )

    if total_source is not None:
        add_coverage_hover_to_figure(cov_figure, total_source, hp_source=hp_source)
    if hp_source is not None and not any(isinstance(t, TapTool) for t in cov_figure.tools):
        cov_figure.add_tools(TapTool())
    return {"total": total_source, "hp": hp_source}


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
    _configure_genomic_x_axis(gene_figure)
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
    _configure_genomic_x_axis(fig)
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


_REPEAT_DENSITY_HEIGHT = 20


def _build_repeat_density_source_data(
    density: "np.ndarray",
    region_start: int,
    region_end: int,
) -> dict:
    """Return ColumnDataSource data dict for the repeat density vbar track."""
    n = len(density)
    span = region_end - region_start
    bin_width = span / n
    return {
        "x": [region_start + (i + 0.5) * bin_width for i in range(n)],
        "top": density.tolist(),
        "width": [bin_width] * n,
    }


def create_repeat_density_figure(
    shared_x_range,
    density: "np.ndarray",
    region_start: int,
    region_end: int,
):
    """Thin bar-chart track showing per-bin off-diagonal repeat density.

    Placed below the gene track and above the genomic x-axis strip. Each bar
    represents one dotplot bin; height is the number of other bins that bin matched
    in the reference self-identity dotplot.

    Returns (figure, ColumnDataSource) so callers can update data on region swap.
    """
    y_max = float(max(density.max(), 1.0))
    source = ColumnDataSource(_build_repeat_density_source_data(density, region_start, region_end))
    y_range = Range1d(0, y_max * 1.1)
    fig = figure(
        height=_REPEAT_DENSITY_HEIGHT,
        x_range=shared_x_range,
        y_range=y_range,
        toolbar_location=None,
        tools=[],
        sizing_mode="stretch_width",
        outline_line_color=None,
    )
    fig.vbar(
        x="x",
        top="top",
        width="width",
        source=source,
        color="#5F249F",
        alpha=0.7,
        line_color=None,
    )
    fig.xgrid.grid_line_color = None
    fig.ygrid.grid_line_color = None
    fig.xaxis.visible = False
    fig.yaxis.minor_tick_line_color = None
    fig.yaxis.major_tick_line_color = None
    fig.yaxis.major_label_text_color = None
    fig.yaxis.axis_line_color = None
    fig.yaxis.axis_label = "Ident"
    fig.yaxis.axis_label_text_font_size = "8px"
    fig.min_border_top = 0
    fig.min_border_bottom = 0
    return fig, source


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
    phase_set_marker_renderers=None,
    default_hide_alignment_numbers=False,
    default_hide_1bp_indels=False,
    read_search_button=None,
    dotplot_thumbnail=None,
    enable_read_filter_checkboxes=False,
    show_checkbox_controls=True,
    model_name_suffix="",
    original_region_widget=None,
    nav_chrom_div=None,
    nav_orig_start_div=None,
    nav_orig_end_div=None,
):
    """Create coordinate displays: static full region at top + editable view row.
    If one_bp_renderers is provided, add a checkbox to hide/show 1bp indels.
    one_bp_markers, one_bp_texts, one_bp_segments restore LOD when re-showing.
    If alignment_label_renderers is provided, add checkbox to hide alignment numbers.
    If phase_set_marker_renderers is provided, add checkbox to hide phase-set markers.
    default_hide_alignment_numbers: if True (e.g. paraphase), numbers hidden initially.
    nav_chrom_div / nav_orig_start_div / nav_orig_end_div: when provided, view and
    go callbacks read from these Div models instead of hardcoded values, allowing the
    slot-swap callback to update navigation bounds after swapping regions.

    Returns (coord_controls, hide_1bp_checkbox, coord_input).
    """
    start_str = f"{coordinate_start:,}"
    end_str = f"{coordinate_end:,}"
    region_size = coordinate_end - coordinate_start + 1
    size_str = format_region_size(region_size)
    original_region_name = "orographer_original_region_wide"
    center_offset_name = "orographer_center_offset_spacer"
    if model_name_suffix:
        original_region_name = f"{original_region_name}_{model_name_suffix}"
        center_offset_name = f"{center_offset_name}_{model_name_suffix}"

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
        name=original_region_name,
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
    # Use caller-supplied Div models when available so swap callback can update nav bounds.
    if nav_chrom_div is None:
        nav_chrom_div = Div(text=chrom, visible=False)
    if nav_orig_start_div is None:
        nav_orig_start_div = Div(text=str(coordinate_start), visible=False)
    if nav_orig_end_div is None:
        nav_orig_end_div = Div(text=str(coordinate_end), visible=False)

    go_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "coord_input": coord_input,
            "error_div": error_div,
            "orig_start_div": nav_orig_start_div,
            "orig_end_div": nav_orig_end_div,
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
            "chrom_div": nav_chrom_div,
            "view_size_div": view_size_div,
        },
        code=load_javascript("view_callback.js"),
        module=True,
    )
    plot_figure.x_range.js_on_change("start", view_callback)
    plot_figure.x_range.js_on_change("end", view_callback)
    plot_figure.x_range.js_on_change("change", view_callback)

    nav_button_styles = {
        "font-size": "16px",
        "font-weight": "bold",
        "line-height": "1",
        "min-width": "0",
        "padding": "0",
        "text-align": "center",
    }
    zoom_out_btn = Button(
        label="\u2212",
        button_type="default",
        width=28,
        height=22,
        margin=(0, 1, 0, 0),
        styles=nav_button_styles.copy(),
    )
    zoom_in_btn = Button(
        label="+",
        button_type="default",
        width=28,
        height=22,
        margin=(0, 0, 0, 1),
        styles=nav_button_styles.copy(),
    )
    pan_left_btn = Button(
        label="\u2039",
        button_type="default",
        width=28,
        height=22,
        margin=(0, 1, 0, 0),
        styles=nav_button_styles.copy(),
    )
    pan_right_btn = Button(
        label="\u203a",
        button_type="default",
        width=28,
        height=22,
        margin=(0, 0, 0, 1),
        styles=nav_button_styles.copy(),
    )
    pan_left_btn.disabled = True
    pan_right_btn.disabled = True
    zoom_out_btn.disabled = True
    zoom_in_btn.disabled = coordinate_end - coordinate_start <= 10
    zoom_in_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "factor": 0.5,
            "orig_start_div": nav_orig_start_div,
            "orig_end_div": nav_orig_end_div,
        },
        code=load_javascript("zoom_buttons_callback.js"),
    )
    zoom_out_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "factor": 2,
            "orig_start_div": nav_orig_start_div,
            "orig_end_div": nav_orig_end_div,
        },
        code=load_javascript("zoom_buttons_callback.js"),
    )
    pan_left_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "direction": -1,
            "fraction": 0.05,
            "orig_start_div": nav_orig_start_div,
            "orig_end_div": nav_orig_end_div,
        },
        code=load_javascript("pan_buttons_callback.js"),
    )
    pan_right_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "direction": 1,
            "fraction": 0.05,
            "orig_start_div": nav_orig_start_div,
            "orig_end_div": nav_orig_end_div,
        },
        code=load_javascript("pan_buttons_callback.js"),
    )
    pan_state_callback = CustomJS(
        args={
            "x_range": plot_figure.x_range,
            "left_button": pan_left_btn,
            "right_button": pan_right_btn,
            "zoom_out_button": zoom_out_btn,
            "zoom_in_button": zoom_in_btn,
            "orig_start_div": nav_orig_start_div,
            "orig_end_div": nav_orig_end_div,
        },
        code=load_javascript("pan_buttons_state_callback.js"),
    )
    zoom_in_btn.js_on_click(zoom_in_callback)
    zoom_out_btn.js_on_click(zoom_out_callback)
    pan_left_btn.js_on_click(pan_left_callback)
    pan_right_btn.js_on_click(pan_right_callback)
    enable_button_click_hold(pan_left_btn)
    enable_button_click_hold(pan_right_btn)
    plot_figure.x_range.js_on_change("start", pan_state_callback)
    plot_figure.x_range.js_on_change("end", pan_state_callback)
    plot_figure.x_range.js_on_change("change", pan_state_callback)
    view_size_with_zoom = column(
        row(zoom_out_btn, zoom_in_btn, spacing=0),
        row(pan_left_btn, pan_right_btn, spacing=0),
        view_size_div,
        spacing=0,
        align="start",
    )

    checkbox_groups, hide_1bp_checkbox, _cursor_guide_checkbox = _plot_checkbox_items(
        plot_figure,
        one_bp_renderers=one_bp_renderers,
        one_bp_markers=one_bp_markers,
        one_bp_texts=one_bp_texts,
        one_bp_segments=one_bp_segments,
        alignment_label_renderers=alignment_label_renderers,
        phase_set_marker_renderers=phase_set_marker_renderers,
        default_hide_alignment_numbers=default_hide_alignment_numbers,
        default_hide_1bp_indels=default_hide_1bp_indels,
        enable_read_filter_checkboxes=enable_read_filter_checkboxes,
        enable_cursor_guide=False,
        checkbox_margin=(0, 8, 0, 0),
    )
    checkbox_items = _flatten_checkbox_groups(checkbox_groups)
    checkbox_controls = Spacer(width=0, height=1)
    if show_checkbox_controls and checkbox_items:
        checkbox_controls = column(*checkbox_items, spacing=0, align="start")

    # Coordinate navigation row. Stretch spacers keep the editable coordinate
    # controls centered while optional metadata blocks sit on the same row.
    input_center_offset = COORD_INPUT_CENTER_OFFSET_PX
    center_offset_spacer = Spacer(
        width=input_center_offset,
        height=1,
        sizing_mode="fixed",
        name=center_offset_name,
    )
    nav_row = row(
        checkbox_controls,
        center_offset_spacer,
        Spacer(sizing_mode="stretch_width"),
        dotplot_thumbnail if dotplot_thumbnail is not None else Spacer(width=0, height=1),
        coord_input,
        go_button,
        view_size_with_zoom,
        read_search_button if read_search_button is not None else Spacer(width=0, height=1),
        Spacer(sizing_mode="stretch_width"),
        original_region_widget if original_region_widget is not None else original_region_block_wide,
        sizing_mode="stretch_width",
        align="start",
        spacing=8,
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

    coord_controls = column(
        nav_row,
        error_row,
        sizing_mode="stretch_width",
        align="start",
    )
    return (
        coord_controls,
        hide_1bp_checkbox,
        coord_input,
    )


def _plot_checkbox_items(
    plot_figure,
    *,
    one_bp_renderers=None,
    one_bp_markers=None,
    one_bp_texts=None,
    one_bp_segments=None,
    alignment_label_renderers=None,
    phase_set_marker_renderers=None,
    enable_cursor_guide=True,
    default_hide_alignment_numbers=False,
    default_hide_1bp_indels=False,
    enable_read_filter_checkboxes=False,
    read_filter_sources=None,
    read_filter_y_ranges=None,
    read_filter_y_bounds=None,
    checkbox_margin=(0, 8, 0, 0),
):
    """Create plot-wide checkbox controls and return them with the 1bp checkbox."""
    display_items = []
    evidence_items = []
    hide_1bp_checkbox = None
    cursor_guide_checkbox = None
    if enable_cursor_guide:
        cursor_guide_checkbox = Checkbox(
            label="Show cursor guide",
            active=True,
            width=145,
            margin=checkbox_margin,
        )
        display_items.append(cursor_guide_checkbox)

    if one_bp_renderers is not None and len(one_bp_renderers) > 0:
        hide_1bp_checkbox = Checkbox(
            label="Hide 1bp INDELs",
            active=default_hide_1bp_indels,
            width=125,
            margin=checkbox_margin,
        )
        if default_hide_1bp_indels:
            for r in one_bp_markers or []:
                r.visible = True
            for r in one_bp_texts or []:
                r.visible = False
            for r in one_bp_segments or []:
                r.visible = False
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
        display_items.append(hide_1bp_checkbox)

    if alignment_label_renderers is not None and len(alignment_label_renderers) > 0:
        hide_alignment_numbers_checkbox = Checkbox(
            label="Hide algn numbers",
            active=default_hide_alignment_numbers,
            width=135,
            margin=checkbox_margin,
        )
        if default_hide_alignment_numbers:
            for r in alignment_label_renderers:
                r.visible = False
        hide_alignment_numbers_callback = CustomJS(
            args={"alignment_label_renderers": alignment_label_renderers},
            code=load_javascript("hide_alignment_numbers_callback.js"),
        )
        hide_alignment_numbers_checkbox.js_on_change("active", hide_alignment_numbers_callback)
        display_items.append(hide_alignment_numbers_checkbox)

    if phase_set_marker_renderers is not None and len(phase_set_marker_renderers) > 0:
        hide_phase_set_markers_checkbox = Checkbox(
            label="Hide phaseset markers",
            active=False,
            width=160,
            margin=checkbox_margin,
        )
        hide_phase_set_markers_callback = CustomJS(
            args={"phase_set_marker_renderers": phase_set_marker_renderers},
            code=load_javascript("hide_phase_set_markers_callback.js"),
        )
        hide_phase_set_markers_checkbox.js_on_change("active", hide_phase_set_markers_callback)
        display_items.append(hide_phase_set_markers_checkbox)

    if enable_read_filter_checkboxes:
        hide_non_split_checkbox = Checkbox(
            label="Show only split algns",
            active=False,
            width=155,
            margin=checkbox_margin,
        )
        hide_non_multiregion_checkbox = Checkbox(
            label="Show only multiregion algns",
            active=False,
            width=200,
            margin=checkbox_margin,
        )
        read_filter_callback = CustomJS(
            args={
                "hide_non_split_checkbox": hide_non_split_checkbox,
                "hide_non_multiregion_checkbox": hide_non_multiregion_checkbox,
                "read_filter_sources": read_filter_sources or [],
                "read_filter_y_ranges": read_filter_y_ranges or [],
                "read_filter_y_bounds": read_filter_y_bounds or [],
            },
            code=load_javascript("read_filter_callback.js"),
        )
        hide_non_split_checkbox.js_on_change("active", read_filter_callback)
        hide_non_multiregion_checkbox.js_on_change("active", read_filter_callback)
        evidence_items.extend([hide_non_split_checkbox, hide_non_multiregion_checkbox])

    return (
        {"Evidence": evidence_items, "Display": display_items},
        hide_1bp_checkbox,
        cursor_guide_checkbox,
    )


def _flatten_checkbox_groups(checkbox_groups):
    """Return checkbox controls in visual group order."""
    return [
        checkbox
        for group_name in ("Evidence", "Display")
        for checkbox in checkbox_groups.get(group_name, [])
    ]


def create_global_checkbox_controls(
    plot_figure,
    dotplot_thumbnail=None,
    read_search_controls=None,
    **kwargs,
):
    """Create a compact global checkbox strip for multi-region controls."""
    checkbox_groups, hide_1bp_checkbox, cursor_guide_checkbox = _plot_checkbox_items(
        plot_figure,
        checkbox_margin=(6, 6, 0, 0),
        **kwargs,
    )
    checkbox_items = _flatten_checkbox_groups(checkbox_groups)

    trailing_controls = []
    if read_search_controls is not None or dotplot_thumbnail is not None:
        trailing_children = [
            control for control in (read_search_controls, dotplot_thumbnail) if control is not None
        ]
        trailing_width = sum(getattr(control, "width", 0) or 0 for control in trailing_children)
        trailing_width += 12 * max(0, len(trailing_children) - 1)
        trailing_controls.append(
            row(
                *trailing_children,
                width=trailing_width,
                height=30,
                spacing=12,
                sizing_mode="fixed",
                align="center",
                margin=(6, 6, 0, 0),
            )
        )

    if not checkbox_items and not trailing_controls:
        return None, hide_1bp_checkbox, cursor_guide_checkbox
    controls = row(
        Spacer(sizing_mode="stretch_width"),
        *checkbox_items,
        *trailing_controls,
        Spacer(sizing_mode="stretch_width"),
        sizing_mode="stretch_width",
        align="center",
        spacing=12,
        margin=(0, 0, 3, 0),
        styles={
            "background": "#f8fafc",
            "border-bottom": "1px solid #e5e7eb",
            "box-sizing": "border-box",
            "column-gap": "12px",
            "flex-wrap": "wrap",
            "justify-content": "center",
            "row-gap": "2px",
            "min-height": "30px",
            "padding": "3px 8px",
        },
    )
    return controls, hide_1bp_checkbox, cursor_guide_checkbox


def _build_separator_raw_data(
    read_names: list,
    read_to_y_bottom: dict,
    read_heights: dict,
    coordinate_start: int,
    coordinate_end: int,
    layout_modes=None,
    read_filter_flags_by_read=None,
) -> dict:
    """Compute separator-line source data without creating Bokeh models.

    Used by add_separator_lines (rendering) and _serialize_region_for_swap (serialization).
    """
    read_filter_flags_by_read = read_filter_flags_by_read or {}
    xs = [[coordinate_start, coordinate_end]]
    ys = [[0, 0]]
    source_read_names = [""]
    has_split_alignment = [True]
    has_multiregion_connection = [True]
    for read_name in read_names:
        y_bottom = read_to_y_bottom[read_name] + read_heights[read_name]
        xs.append([coordinate_start, coordinate_end])
        ys.append([y_bottom, y_bottom])
        has_split, has_multiregion = read_filter_flags_by_read.get(read_name, (True, True))
        source_read_names.append(read_name)
        has_split_alignment.append(has_split)
        has_multiregion_connection.append(has_multiregion)
    source_data = {
        "xs": xs,
        "ys": ys,
        "read_name": source_read_names,
        "has_split_alignment": has_split_alignment,
        "has_multiregion_connection": has_multiregion_connection,
        "separator_line_alpha": [0.3] * len(xs),
    }
    if layout_modes:
        for mode in READ_FILTER_LAYOUT_MODES:
            mode_layout = layout_modes.get(mode, {})
            mode_bottoms = mode_layout.get("read_to_y_bottom", {})
            mode_heights = mode_layout.get("read_heights", {})
            mode_ys = [[0, 0]]
            for read_name in read_names:
                if read_name in mode_bottoms and read_name in mode_heights:
                    y_bottom = mode_bottoms[read_name] + mode_heights[read_name]
                    mode_ys.append([y_bottom, y_bottom])
                else:
                    mode_ys.append([0, 0])
            source_data[_layout_mode_column("ys", mode)] = mode_ys
    add_read_filter_visibility_columns(source_data)
    return source_data


def add_separator_lines(
    plot_figure,
    read_names,
    read_to_y_bottom,
    read_heights,
    coordinate_start,
    coordinate_end,
    layout_modes=None,
    read_filter_flags_by_read=None,
):
    """Add horizontal dotted lines to separate reads (single multi_line glyph)."""
    source_data = _build_separator_raw_data(
        read_names,
        read_to_y_bottom,
        read_heights,
        coordinate_start,
        coordinate_end,
        layout_modes,
        read_filter_flags_by_read,
    )
    source = ColumnDataSource(source_data)
    plot_figure.multi_line(
        xs="xs",
        ys="ys",
        source=source,
        line_color="grey",
        line_width=0.5,
        line_dash="dotted",
        line_alpha="separator_line_alpha",
        level="underlay",
    )
    return source


def get_haplotype_label(haplotype):
    """Generate a label for a haplotype value."""
    if haplotype == 0:
        return "Unassigned"
    return f"{haplotype}"


def _build_hp_label_raw_data(
    group_boundaries,
    haplotype_order,
    coordinate_start,
    coordinate_end,
    haplotype_color_fn=None,
    layout_modes=None,
) -> tuple[dict | None, dict | None]:
    """Compute HP label source data dicts without creating Bokeh models.

    Returns (separator_data, label_data) where each is a populated dict or None when empty.
    Used by add_haplotype_labels (rendering) and _serialize_region_for_swap (serialization).
    """
    label_x = []
    label_y = []
    label_texts = []
    label_colors = []
    label_haplotypes = []
    separator_xs = []
    separator_ys = []
    separator_haplotype_pairs = []
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
        label_colors.append(haplotype_color_fn(haplotype) if haplotype_color_fn else "black")
        label_haplotypes.append(haplotype)

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
                separator_haplotype_pairs.append((haplotype, next_haplotype))

    separator_data = None
    if separator_xs:
        separator_data = {
            "xs": separator_xs,
            "ys": separator_ys,
            "read_layout_mode": ["all"] * len(separator_xs),
        }
        if layout_modes:
            for mode in READ_FILTER_LAYOUT_MODES:
                mode_boundaries = layout_modes.get(mode, {}).get("group_boundaries", {})
                mode_ys = []
                for row_index, (haplotype, next_haplotype) in enumerate(separator_haplotype_pairs):
                    if haplotype in mode_boundaries and next_haplotype in mode_boundaries:
                        y_end = mode_boundaries[haplotype][1]
                        next_y_start = mode_boundaries[next_haplotype][0]
                        separator_y = (y_end + next_y_start) / 2
                    else:
                        separator_y = separator_ys[row_index][0]
                    mode_ys.append([separator_y, separator_y])
                separator_data[_layout_mode_column("ys", mode)] = mode_ys
        add_read_filter_visibility_columns(separator_data)

    label_data = None
    if label_x:
        swatch_offset = (coordinate_end - coordinate_start) * 0.006
        text_x = [x + swatch_offset for x in label_x]
        label_data = {
            "x": label_x,
            "y": label_y,
            "text_x": text_x,
            "text": label_texts,
            "color": label_colors,
            "read_layout_mode": ["all"] * len(label_x),
        }
        if layout_modes:
            for mode in READ_FILTER_LAYOUT_MODES:
                mode_boundaries = layout_modes.get(mode, {}).get("group_boundaries", {})
                mode_y = []
                for row_index, haplotype in enumerate(label_haplotypes):
                    if haplotype in mode_boundaries:
                        y_start, y_end = mode_boundaries[haplotype]
                        mode_y.append((y_start + y_end) / 2)
                    else:
                        mode_y.append(label_y[row_index])
                label_data[_layout_mode_column("y", mode)] = mode_y
        add_read_filter_visibility_columns(label_data)

    return separator_data, label_data


def add_haplotype_labels(
    plot_figure,
    group_boundaries,
    haplotype_order,
    coordinate_start,
    coordinate_end,
    haplotype_color_fn=None,
    layout_modes=None,
) -> dict:
    """Add text labels and separator lines for each haplotype group.

    haplotype_color_fn: optional callable(hp_int) -> color_str. When provided
    the label text is colored to match the per-HP coverage line, creating a
    direct visual link between the read-track label and the coverage track.

    Returns a dict {"separator": ColumnDataSource|None, "label": ColumnDataSource|None}.
    """
    separator_data, label_data = _build_hp_label_raw_data(
        group_boundaries, haplotype_order, coordinate_start, coordinate_end,
        haplotype_color_fn, layout_modes,
    )
    separator_source = None
    label_source = None

    if separator_data is not None:
        separator_source = ColumnDataSource(separator_data)
        plot_figure.multi_line(
            xs="xs",
            ys="ys",
            source=separator_source,
            line_color=PLOT_CONFIG["sample_label_color"],
            line_width=1.0,
            line_alpha=1,
            level="underlay",
        )
    if label_data is not None:
        label_source = ColumnDataSource(data=label_data)
        plot_figure.scatter(
            x="x",
            y="y",
            marker="square",
            size=12,
            fill_color="color",
            line_color="white",
            line_width=0.5,
            source=label_source,
        )
        plot_figure.text(
            x="text_x",
            y="y",
            text="text",
            source=label_source,
            text_font_size="11pt",
            text_color="black",
            text_font_style="bold",
            text_align="left",
            text_baseline="middle",
            text_alpha=1.0,
        )
    return {"separator": separator_source, "label": label_source}


def _build_gene_track_raw_data(
    gene_annotations, gene_track_y_start, coordinate_start, coordinate_end, plot_width_px=None
) -> tuple[float, dict, dict, dict, dict, dict]:
    """Compute gene track data dicts without creating Bokeh models.

    Returns (height, body_data, exon_data, intron_data, arrow_data, label_data).
    Used by add_gene_track (rendering) and _serialize_region_for_swap (serialization).
    """
    if plot_width_px is None:
        plot_width_px = GENE_TRACK_PLOT_WIDTH_PX
    mode = _initial_gene_track_mode(coordinate_start, coordinate_end)
    genes = _deduplicate_gene_annotations(gene_annotations)
    gene_row_assignments, num_rows = _assign_gene_rows(
        genes,
        coordinate_start,
        coordinate_end,
        plot_width_px,
        reserve_label_width=mode != "overview",
    )
    row_height = GENE_EXON_HEIGHT + GENE_ROW_SPACING

    body_data = {
        "x0": [],
        "x1": [],
        "y": [],
        "gene_id": [],
        "gene_name": [],
        "strand": [],
        "coordinates": [],
        "line_alpha": [],
        "line_alpha_overview": [],
        "line_alpha_medium": [],
        "line_alpha_detail": [],
    }
    exon_data = {
        "left": [],
        "right": [],
        "top": [],
        "bottom": [],
        "gene_name": [],
        "gene_strand": [],
        "exon_number": [],
        "exon_coordinates": [],
        "isoform_coordinates": [],
        "representative_transcript": [],
        "representative_selection_method": [],
        "fill_alpha": [],
        "fill_alpha_overview": [],
        "fill_alpha_medium": [],
        "fill_alpha_detail": [],
        "line_alpha": [],
        "line_alpha_overview": [],
        "line_alpha_medium": [],
        "line_alpha_detail": [],
    }
    intron_data = {
        "xs": [],
        "ys": [],
        "line_alpha": [],
        "line_alpha_overview": [],
        "line_alpha_medium": [],
        "line_alpha_detail": [],
    }
    arrow_data = {
        "x": [],
        "y": [],
        "angle": [],
        "fill_alpha": [],
        "fill_alpha_overview": [],
        "fill_alpha_medium": [],
        "fill_alpha_detail": [],
        "line_alpha": [],
        "line_alpha_overview": [],
        "line_alpha_medium": [],
        "line_alpha_detail": [],
    }
    label_data = {
        "x": [],
        "y": [],
        "text": [],
        "text_alpha": [],
        "text_alpha_overview": [],
        "text_alpha_medium": [],
        "text_alpha_detail": [],
    }

    for gene in genes:
        row_idx = gene_row_assignments.get(gene.gene_id, 0)
        gene_y_center = gene_track_y_start + (row_idx + 0.5) * row_height
        sorted_exons = sorted(gene.exons, key=lambda exon: (exon[0], exon[1], exon[2]))
        gene_left, gene_right = _annotation_interval(gene.start, gene.end)
        gene_vis_start = max(gene_left, coordinate_start)
        gene_vis_end = min(gene_right, coordinate_end + 1)
        if gene_vis_start < gene_vis_end:
            body_alpha = _mode_values(
                mode,
                overview=GENE_OVERVIEW_BODY_ALPHA,
                medium=0.45,
                detail=0.25,
            )
            body_data["x0"].append(gene_vis_start)
            body_data["x1"].append(gene_vis_end)
            body_data["y"].append(gene_y_center)
            body_data["gene_id"].append(gene.gene_id)
            body_data["gene_name"].append(gene.gene_name)
            body_data["strand"].append(gene.strand or "unknown")
            body_data["coordinates"].append(
                _format_genomic_coordinates(gene.chrom, gene.start, gene.end)
            )
            body_data["line_alpha"].append(body_alpha[0])
            body_data["line_alpha_overview"].append(body_alpha[1])
            body_data["line_alpha_medium"].append(body_alpha[2])
            body_data["line_alpha_detail"].append(body_alpha[3])

        for exon_start, exon_end, exon_number in sorted_exons:
            exon_left, exon_right = _annotation_interval(exon_start, exon_end)
            vis_start = max(exon_left, coordinate_start)
            vis_end = min(exon_right, coordinate_end + 1)
            if vis_start < vis_end:
                fill_alpha = _mode_values(
                    mode,
                    overview=GENE_HIDDEN_ALPHA,
                    medium=GENE_MEDIUM_EXON_ALPHA,
                    detail=GENE_DETAIL_EXON_ALPHA,
                )
                line_alpha = _mode_values(
                    mode,
                    overview=GENE_HIDDEN_ALPHA,
                    medium=0.8,
                    detail=1.0,
                )
                exon_data["left"].append(vis_start)
                exon_data["right"].append(vis_end)
                exon_data["top"].append(gene_y_center - GENE_EXON_HEIGHT / 2)
                exon_data["bottom"].append(gene_y_center + GENE_EXON_HEIGHT / 2)
                exon_data["gene_name"].append(gene.gene_name)
                strand_token = gene.strand or ""
                if strand_token in ("", "."):
                    exon_data["gene_strand"].append("unknown")
                else:
                    exon_data["gene_strand"].append(_gene_strand_modal_label(strand_token))
                exon_data["exon_number"].append(str(exon_number))
                exon_data["exon_coordinates"].append(
                    _format_genomic_coordinates(
                        gene.chrom,
                        exon_start,
                        exon_end,
                        include_size=True,
                    )
                )
                exon_data["isoform_coordinates"].append(
                    _format_genomic_coordinates(
                        gene.chrom,
                        gene.start,
                        gene.end,
                        include_size=True,
                    )
                )
                transcript_label = gene.representative_transcript_name or ""
                if gene.representative_transcript_id:
                    transcript_label = gene.representative_transcript_id
                exon_data["representative_transcript"].append(transcript_label)
                exon_data["representative_selection_method"].append(
                    gene.representative_selection_method
                )
                exon_data["fill_alpha"].append(fill_alpha[0])
                exon_data["fill_alpha_overview"].append(fill_alpha[1])
                exon_data["fill_alpha_medium"].append(fill_alpha[2])
                exon_data["fill_alpha_detail"].append(fill_alpha[3])
                exon_data["line_alpha"].append(line_alpha[0])
                exon_data["line_alpha_overview"].append(line_alpha[1])
                exon_data["line_alpha_medium"].append(line_alpha[2])
                exon_data["line_alpha_detail"].append(line_alpha[3])

        for i in range(len(sorted_exons) - 1):
            intron_start = sorted_exons[i][1] + 1
            intron_end = sorted_exons[i + 1][0]
            if intron_end > coordinate_start and intron_start < coordinate_end:
                vis_start = max(intron_start, coordinate_start)
                vis_end = min(intron_end, coordinate_end)
                if vis_start < vis_end:
                    _append_gene_intron(
                        intron_data,
                        start=vis_start,
                        end=vis_end,
                        y_center=gene_y_center,
                        mode=mode,
                    )
                    _append_gene_arrows(
                        arrow_data,
                        start=vis_start,
                        end=vis_end,
                        y_center=gene_y_center,
                        strand=gene.strand,
                        mode=mode,
                    )

        if sorted_exons:
            first_exon_left, _first_exon_right = _annotation_interval(
                sorted_exons[0][0], sorted_exons[0][1]
            )
            _last_exon_left, last_exon_right = _annotation_interval(
                sorted_exons[-1][0], sorted_exons[-1][1]
            )
            _append_gene_intron(
                intron_data,
                start=max(gene_left, coordinate_start),
                end=min(first_exon_left, coordinate_end),
                y_center=gene_y_center,
                mode=mode,
            )
            _append_gene_intron(
                intron_data,
                start=max(last_exon_right, coordinate_start),
                end=min(gene_right, coordinate_end),
                y_center=gene_y_center,
                mode=mode,
            )

        if len(sorted_exons) == 1:
            exon_start, exon_end, _ = sorted_exons[0]
            exon_left, exon_right = _annotation_interval(exon_start, exon_end)
            vis_start = max(exon_left, coordinate_start)
            vis_end = min(exon_right, coordinate_end + 1)
            exon_length = vis_end - vis_start
            if exon_length > GENE_MIN_EXON_ARROW_BP:
                _append_gene_arrows(
                    arrow_data,
                    start=vis_start,
                    end=vis_end,
                    y_center=gene_y_center,
                    strand=gene.strand,
                    mode=mode,
                )

        if gene_vis_start < gene_vis_end:
            label_alpha = _mode_values(
                mode,
                overview=GENE_HIDDEN_ALPHA,
                medium=GENE_HIDDEN_ALPHA,
                detail=GENE_LABEL_ALPHA,
            )
            label_data["x"].append((gene_vis_start + gene_vis_end) / 2)
            label_data["y"].append(gene_y_center - GENE_EXON_HEIGHT / 2)
            label_data["text"].append(gene.gene_name)
            label_data["text_alpha"].append(label_alpha[0])
            label_data["text_alpha_overview"].append(label_alpha[1])
            label_data["text_alpha_medium"].append(label_alpha[2])
            label_data["text_alpha_detail"].append(label_alpha[3])

    return num_rows * row_height + 1.0, body_data, exon_data, intron_data, arrow_data, label_data


def add_gene_track(
    plot_figure, gene_annotations, gene_track_y_start, coordinate_start, coordinate_end
) -> tuple[float, dict]:
    """Add a zoom-adaptive gene annotation track.

    Returns (height, sources) where sources maps "body"/"exon"/"intron"/"arrow"/"label"
    to their ColumnDataSource (or None when the data set was empty).
    """
    if not gene_annotations:
        return 0, {}

    gene_color = "#3366CC"
    plot_width_px = int(getattr(plot_figure, "width", GENE_TRACK_PLOT_WIDTH_PX))
    height, body_data, exon_data, intron_data, arrow_data, label_data = _build_gene_track_raw_data(
        gene_annotations, gene_track_y_start, coordinate_start, coordinate_end, plot_width_px
    )

    callback_sources = []
    body_source = None
    intron_source = None
    exon_source = None
    arrow_source = None
    label_source = None

    if body_data["x0"]:
        body_source = ColumnDataSource(data=body_data)
        body_renderer = plot_figure.segment(
            x0="x0",
            y0="y",
            x1="x1",
            y1="y",
            source=body_source,
            line_color=gene_color,
            line_width=3,
            line_alpha="line_alpha",
        )
        callback_sources.append((body_source, "line_alpha"))
        plot_figure.add_tools(
            HoverTool(
                renderers=[body_renderer],
                tooltips=[
                    ("Gene", "@gene_name"),
                    ("ID", "@gene_id"),
                    ("Strand", "@strand"),
                    ("Coordinates", "@coordinates"),
                ],
            )
        )

    if intron_data["xs"]:
        intron_source = ColumnDataSource(data=intron_data)
        plot_figure.multi_line(
            xs="xs",
            ys="ys",
            source=intron_source,
            line_color=gene_color,
            line_width=1.5,
            line_alpha="line_alpha",
        )
        callback_sources.append((intron_source, "line_alpha"))

    if exon_data["left"]:
        exon_source = ColumnDataSource(data=exon_data)
        exon_renderer = plot_figure.quad(
            left="left",
            right="right",
            top="top",
            bottom="bottom",
            source=exon_source,
            fill_color=gene_color,
            line_color=gene_color,
            fill_alpha="fill_alpha",
            line_alpha="line_alpha",
        )
        callback_sources.append((exon_source, "fill_alpha"))
        callback_sources.append((exon_source, "line_alpha"))
        tap_tool = TapTool()
        tap_tool.renderers = [exon_renderer]
        plot_figure.add_tools(tap_tool)
        plot_figure.toolbar.active_tap = tap_tool
        exon_click_callback = get_exon_click_callback(exon_source)
        exon_source.selected.js_on_change("indices", exon_click_callback)

    if arrow_data["x"]:
        arrow_source = ColumnDataSource(data=arrow_data)
        plot_figure.scatter(
            x="x",
            y="y",
            source=arrow_source,
            marker="triangle",
            size=8,
            angle="angle",
            fill_color="white",
            line_color=gene_color,
            fill_alpha="fill_alpha",
            line_alpha="line_alpha",
        )
        callback_sources.append((arrow_source, "fill_alpha"))
        callback_sources.append((arrow_source, "line_alpha"))

    if label_data["x"]:
        label_source = ColumnDataSource(data=label_data)
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
            text_alpha="text_alpha",
        )
        callback_sources.append((label_source, "text_alpha"))

    if callback_sources:
        _add_gene_track_zoom_callback(plot_figure, callback_sources)
    gene_sources = {
        "body": body_source,
        "intron": intron_source,
        "exon": exon_source,
        "arrow": arrow_source,
        "label": label_source,
    }
    return height, gene_sources


def _cluster_insertion_sites(sites: list[dict], cluster_bp: float) -> list[list[dict]]:
    """Cluster insertion sites by sorted reference position."""
    sorted_sites = sorted(sites, key=lambda site: site["pos"])
    clustered_sites = []
    for site in sorted_sites:
        if clustered_sites and site["pos"] - clustered_sites[-1][-1]["pos"] <= cluster_bp:
            clustered_sites[-1].append(site)
        else:
            clustered_sites.append([site])
    return clustered_sites


def _round_insertion_size(size: float) -> int:
    """Round positive insertion sizes to the nearest integer for display."""
    return math.floor(size + 0.5)


def _weighted_median_size_by_read_count(cluster: list[dict]) -> int:
    """Return the median insertion size weighted by per-site read count."""
    if not cluster:
        return 0
    values = sorted(
        (
            site["median_size"],
            max(1, site.get("total_count", site.get("count", 1))),
        )
        for site in cluster
    )
    total_weight = sum(weight for _size, weight in values)
    midpoint = total_weight / 2
    cumulative = 0
    for idx, (size, weight) in enumerate(values):
        previous_cumulative = cumulative
        cumulative += weight
        if cumulative == midpoint and idx + 1 < len(values):
            return _round_insertion_size((size + values[idx + 1][0]) / 2)
        if cumulative > midpoint or previous_cumulative <= midpoint <= cumulative:
            return _round_insertion_size(size)
    return _round_insertion_size(values[-1][0])


def _insertion_cluster_source_data(
    clustered_sites: list[list[dict]],
    hp: int,
    y_pos: float,
    y_by_layout_mode: dict[str, float] | None = None,
) -> dict:
    """Build ColumnDataSource data for displayed insertion marker clusters."""
    xs = [round((cluster[0]["pos"] + cluster[-1]["pos"]) / 2) for cluster in clustered_sites]
    weighted_median_sizes = [
        _weighted_median_size_by_read_count(cluster) for cluster in clustered_sites
    ]
    marker_labels = [
        f"{weighted_median_size:g}"
        if len(cluster) == 1
        else f"{weighted_median_size:g}({len(cluster)})"
        for cluster, weighted_median_size in zip(
            clustered_sites, weighted_median_sizes, strict=True
        )
    ]
    marker_widths = [
        max(
            PLOT_CONFIG["complex_sv_insertion_marker_min_width"],
            PLOT_CONFIG["complex_sv_insertion_marker_char_px"] * len(label) + 10,
        )
        for label in marker_labels
    ]
    marker_height = PLOT_CONFIG["complex_sv_insertion_marker_height"]
    source_data = {
        "x": xs,
        "y": [y_pos] * len(xs),
        "pos": [cluster[0]["pos"] for cluster in clustered_sites],
        "count": [sum(site["count"] for site in cluster) for cluster in clustered_sites],
        "median_size": weighted_median_sizes,
        "marker_label": marker_labels,
        "marker_width": marker_widths,
        "marker_height": [marker_height] * len(xs),
        "top_names": [
            "\n".join(name for site in cluster for name in site["read_names"][:5])
            for cluster in clustered_sites
        ],
        "total_count": [
            sum(site["total_count"] for site in cluster) for cluster in clustered_sites
        ],
        "hp_label": ["Unassigned" if hp == 0 else f"HP{hp}"] * len(xs),
        "chrom": [cluster[0]["chrom"] for cluster in clustered_sites],
        "cluster_pos": [[site["pos"] for site in cluster] for cluster in clustered_sites],
        "cluster_count": [[site["count"] for site in cluster] for cluster in clustered_sites],
        "cluster_median_size": [
            [site["median_size"] for site in cluster] for cluster in clustered_sites
        ],
        "cluster_top_names": [
            ["\n".join(site["read_names"]) for site in cluster] for cluster in clustered_sites
        ],
        "cluster_total_count": [
            [site["total_count"] for site in cluster] for cluster in clustered_sites
        ],
        "cluster_chrom": [[site["chrom"] for site in cluster] for cluster in clustered_sites],
        "read_layout_mode": ["all"] * len(xs),
    }
    y_by_layout_mode = y_by_layout_mode or {}
    for mode in READ_FILTER_LAYOUT_MODES:
        source_data[_layout_mode_column("y", mode)] = [y_by_layout_mode.get(mode, y_pos)] * len(xs)
    add_read_filter_visibility_columns(source_data)
    return source_data


def _insertion_raw_source_data(
    sites: list[dict],
    hp: int,
    y_pos: float,
    y_by_layout_mode: dict[str, float] | None = None,
) -> dict:
    """Build raw insertion-site data for browser-side reclustering."""
    hp_label = "Unassigned" if hp == 0 else f"HP{hp}"
    source_data = {
        "pos": [site["pos"] for site in sites],
        "count": [site["count"] for site in sites],
        "median_size": [site["median_size"] for site in sites],
        "top_names": ["\n".join(site["read_names"][:5]) for site in sites],
        "total_count": [site["total_count"] for site in sites],
        "hp_label": [hp_label] * len(sites),
        "chrom": [site["chrom"] for site in sites],
        "y": [y_pos] * len(sites),
        "read_layout_mode": ["all"] * len(sites),
    }
    y_by_layout_mode = y_by_layout_mode or {}
    for mode in READ_FILTER_LAYOUT_MODES:
        source_data[_layout_mode_column("y", mode)] = [y_by_layout_mode.get(mode, y_pos)] * len(
            sites
        )
    add_read_filter_visibility_columns(source_data)
    return source_data


def add_insertion_markers_to_figure(
    plot_figure,
    sites: list[dict],
    hp: int,
    y_pos: float,
    coordinate_start: int,
    coordinate_end: int,
    y_by_layout_mode: dict[str, float] | None = None,
    color: str | None = None,
):
    """Draw dynamically clustered insertion-site markers. Returns the marker renderer."""
    region_span = max(1, coordinate_end - coordinate_start + 1)
    target_px = PLOT_CONFIG["complex_sv_insertion_marker_target_px"]
    cluster_bp = region_span / plot_figure.width * target_px
    clustered_sites = _cluster_insertion_sites(sites, cluster_bp)
    raw_source = ColumnDataSource(_insertion_raw_source_data(sites, hp, y_pos, y_by_layout_mode))
    source = ColumnDataSource(
        _insertion_cluster_source_data(clustered_sites, hp, y_pos, y_by_layout_mode)
    )
    marker_color = color or get_haplotype_color(hp)
    renderer = plot_figure.rect(
        x="x",
        y="y",
        width="marker_width",
        height="marker_height",
        width_units="screen",
        height_units="screen",
        fill_color=marker_color,
        line_color=marker_color,
        fill_alpha=PLOT_CONFIG["insertion_fill_alpha"],
        line_alpha=PLOT_CONFIG["insertion_line_alpha"],
        source=source,
    )
    plot_figure.text(
        x="x",
        y="y",
        text="marker_label",
        source=source,
        text_font_size="9pt",
        text_font_style="bold",
        text_color="black",
        text_align="center",
        text_baseline="middle",
    )
    callback = CustomJS(
        args={"source": source},
        code=load_javascript("insertion_click_callback.js"),
    )
    source.selected.js_on_change("indices", callback)
    recluster_callback = CustomJS(
        args={
            "raw_source": raw_source,
            "display_source": source,
            "x_range": plot_figure.x_range,
            "plot_width": plot_figure.width,
            "target_px": target_px,
            "marker_height": PLOT_CONFIG["complex_sv_insertion_marker_height"],
            "marker_min_width": PLOT_CONFIG["complex_sv_insertion_marker_min_width"],
            "marker_char_px": PLOT_CONFIG["complex_sv_insertion_marker_char_px"],
        },
        code=load_javascript("insertion_marker_cluster_callback.js"),
    )
    plot_figure.x_range.js_on_change("start", recluster_callback)
    plot_figure.x_range.js_on_change("end", recluster_callback)
    plot_figure.x_range.js_on_change("change", recluster_callback)
    return renderer, [raw_source, source]


_DOTPLOT_CMAP = LinearSegmentedColormap.from_list(
    "orographer_dotplot",
    ["#ffffff", "#000000"],
).with_extremes(under="#d8d8d8")  # N-masked bins (sentinel -1.0) render as light gray


def _nice_tick_interval(span: int, target_count: int = 4) -> int:
    """Return a round tick interval that produces approximately target_count ticks."""
    if span <= 0:
        return 1
    raw = span / target_count
    magnitude = 10 ** math.floor(math.log10(max(raw, 1)))
    candidates = [1, 2, 5, 10]
    interval = min(candidates, key=lambda m: abs(m * magnitude - raw)) * magnitude
    return max(1, int(interval))


def _format_genomic_pos(pos: int) -> str:
    """Return a human-readable genomic position string."""
    if pos >= 1_000_000:
        return f"{pos / 1_000_000:.2f} Mb"
    if pos >= 1_000:
        return f"{pos / 1_000:.1f} kb"
    return f"{pos} bp"


def _build_block_ticks(blocks: list[dict]) -> tuple[list[float], list[str]]:
    """Return tick positions (in offset space) and labels for a multi-block dotplot."""
    positions: list[float] = []
    labels: list[str] = []
    for block in blocks:
        offset_span = block["offset_end"] - block["offset_start"]
        if offset_span <= 0:
            continue
        interval = _nice_tick_interval(offset_span)
        tick_genomic = math.ceil(block["start"] / interval) * interval
        while True:
            offset_pos = block["offset_start"] + (tick_genomic - block["start"])
            if offset_pos >= block["offset_end"]:
                break
            if offset_pos >= block["offset_start"]:
                positions.append(float(offset_pos))
                labels.append(f"{block['chromosome']}:{_format_genomic_pos(tick_genomic)}")
            tick_genomic += interval
    return positions, labels


def _render_dotplot_png_data_url(
    matrix: np.ndarray,
    region,
    *,
    with_axes: bool,
    figure_size: tuple[float, float],
    dpi: int,
    blocks: list[dict] | None = None,
) -> str:
    """Render a dotplot matrix PNG with Matplotlib."""
    fig = Figure(figsize=figure_size, dpi=dpi)
    canvas = FigureCanvasAgg(fig)
    axis = fig.add_subplot(111)
    if blocks:
        extent = (0, blocks[-1]["offset_end"], 0, blocks[-1]["offset_end"])
    else:
        extent = (region.start, region.end, region.start, region.end)
    axis.imshow(
        matrix,
        extent=extent,
        origin="lower",
        interpolation="nearest",
        aspect="equal",
        cmap=_DOTPLOT_CMAP,
        vmin=0.0,
        vmax=1.0,
    )
    if blocks:
        boundary_color = "#d36b9f"
        resolution = matrix.shape[0]
        total_len = blocks[-1]["offset_end"]
        for block in blocks[:-1]:
            # Align with the left edge of the first bin belonging to the next region.
            # The bin straddling the junction is occupied by the next region's k-mers
            # (not N-masked), so the grey area ends at this bin's left edge, not at
            # the exact sequence offset.
            first_next_bin = block["offset_end"] * resolution // total_len
            boundary = first_next_bin * total_len / resolution
            axis.axvline(boundary, color=boundary_color, linewidth=0.8, alpha=0.8)
            axis.axhline(boundary, color=boundary_color, linewidth=0.8, alpha=0.8)
    if not with_axes:
        axis.set_axis_off()
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        image_buffer = BytesIO()
        canvas.print_png(image_buffer)
        encoded = base64.b64encode(image_buffer.getvalue()).decode("ascii")
        return f"data:image/png;base64,{encoded}"

    if blocks:
        tick_positions, tick_labels = _build_block_ticks(blocks)
        axis.set_xticks(tick_positions)
        axis.set_yticks(tick_positions)
        axis.set_xticklabels(tick_labels, fontsize=7)
        axis.set_yticklabels(tick_labels, fontsize=7)
        axis.tick_params(axis="x", labelrotation=45)
    else:
        axis.set_xlabel("Reference position")
        axis.set_ylabel("Reference position")
        pos_formatter = FuncFormatter(lambda v, _p: _format_genomic_pos(int(v)))
        axis.xaxis.set_major_formatter(pos_formatter)
        axis.yaxis.set_major_formatter(pos_formatter)
        axis.tick_params(axis="x", labelrotation=45)
    axis.grid(False)
    fig.tight_layout()

    image_buffer = BytesIO()
    canvas.print_png(image_buffer)
    encoded = base64.b64encode(image_buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def create_dotplot_thumbnail(
    main_figure,
    matrix: np.ndarray,
    region,
    size: int = 42,
    blocks: list[dict] | None = None,
    region_label: str | None = None,
    modal_title: str = "Reference self-identity dotplot",
    individual_payloads: list[dict] | None = None,
):
    """Create a compact clickable dotplot thumbnail for the coordinate toolbar."""
    thumbnail_url = _render_dotplot_png_data_url(
        matrix,
        region,
        with_axes=False,
        figure_size=(1.0, 1.0),
        dpi=size * 2,
        blocks=blocks,
    )
    modal_url = _render_dotplot_png_data_url(
        matrix,
        region,
        with_axes=True,
        figure_size=(6.4, 6.4),
        dpi=120,
        blocks=blocks,
    )
    individual_images: list[dict] = []
    for ind in individual_payloads or []:
        ind_url = _render_dotplot_png_data_url(
            ind["matrix"],
            ind["region"],
            with_axes=True,
            figure_size=(6.4, 6.4),
            dpi=120,
            blocks=None,
        )
        ind_region = ind["region"]
        ind_span = ind_region.end - ind_region.start + 1
        ind_bin_size = round(ind_span / 512)
        individual_images.append(
            {
                "url": ind_url,
                "label": ind["label"],
                "title": ind["title"],
                "bin_size_bp": ind_bin_size,
            }
        )
    thumbnail_background = (
        f"url({thumbnail_url}), linear-gradient(#ffffff, #ffffff), "
        "linear-gradient(#6b7280, #6b7280)"
    )
    dot_button = Button(
        label="Show ref identity",
        width=156,
        height=30,
        margin=(0, 0, 0, 0),
        sizing_mode="fixed",
        styles={
            "cursor": "pointer",
            "min-width": "156px",
        },
        stylesheets=[
            f"""
            .bk-btn {{
                background-color: transparent;
                background-image: {thumbnail_background};
                background-origin: border-box;
                background-position: 9px center, 8px center, 7px center;
                background-repeat: no-repeat;
                background-size: 22px 22px, 24px 24px, 26px 26px;
                border: 0;
                border-radius: 3px;
                box-shadow: none;
                color: #1f1f1f;
                cursor: pointer;
                font-family: Arial, sans-serif;
                font-size: 13px;
                font-weight: 400;
                height: 30px;
                line-height: 30px;
                min-width: 156px;
                padding: 0 8px 0 42px;
                text-align: left;
            }}
            .bk-btn:hover {{
                background-color: #f3f4f6;
            }}
            """,
        ],
    )
    if blocks:
        total_span = blocks[-1]["offset_end"]
    else:
        total_span = region.end - region.start + 1
    bin_size_bp = round(total_span / 512)
    dot_button.js_on_click(
        CustomJS(
            args={
                "image_url": modal_url,
                "region_label": region_label
                or region.coordinate_str
                or f"{region.chromosome}:{region.start:,}-{region.end:,}",
                "title_text": modal_title,
                "individual_images": individual_images,
                "bin_size_bp": bin_size_bp,
            },
            code=load_javascript("dotplot_modal_callback.js"),
        ),
    )
    return dot_button
