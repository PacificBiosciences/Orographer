"""IsoSeq-specific Bokeh track construction and global controls."""

import logging
import math
from typing import Any

from bokeh.layouts import row
from bokeh.models import Button, Checkbox, ColumnDataSource, CustomJS, Div, Select, TapTool

from ..utils import ISOSEQ_REGION_TYPE
from . import isoseq as isoseq_data
from .figures import (
    add_coverage_hover_to_figure,
    create_bokeh_figure,
    create_bokeh_figure_shared_x,
    create_coverage_figure,
    update_coverage_y_range,
)
from .isoseq_store import prepare_lazy_chunks
from .renderers import add_arrows_to_plot, add_clickable_labels, add_variants_to_plot
from .sources import empty_arrow_source_data
from .utils import PLOT_CONFIG, load_javascript, tap_renderers

logger = logging.getLogger(__name__)


def _stable_gene_color(gene_id: str) -> str:
    return isoseq_data.stable_gene_color(gene_id)


def _format_span(start: int, end: int) -> str:
    return isoseq_data.format_span(start, end)


def _annotation_interval(start: int, end: int) -> tuple[int, int]:
    return start, end + 1


def _intron_direction_positions(
    intron_start: int,
    intron_end: int,
    region_span: int,
) -> list[float]:
    return isoseq_data.intron_direction_positions(intron_start, intron_end, region_span)


def _natural_sort_key(value: Any) -> list:
    return isoseq_data.natural_sort_key(value)


def _empty_arrow_source_data() -> dict:
    return empty_arrow_source_data()


def _build_isoseq_layout(groups: list[dict], segments_by_read: dict) -> dict:
    return isoseq_data.build_layout(groups, segments_by_read)


def _prepare_isoseq_lazy_chunks(
    groups: list[dict],
    segments_by_read: dict,
    layout: dict,
    coordinate_start: int,
    coordinate_end: int,
    region_idx: int,
    row_index: int,
    sample_label: str | None,
    _plot_dom_class: str,
    _plot_model_id: str,
    chunk_dir: str | None = None,
    chunk_url_prefix: str | None = None,
) -> dict:
    return prepare_lazy_chunks(
        groups,
        segments_by_read,
        layout,
        coordinate_start,
        coordinate_end,
        region_idx,
        row_index,
        sample_label,
        chunk_dir=chunk_dir,
        chunk_url_prefix=chunk_url_prefix,
    )


def add_isoseq_transcripts_to_plot(
    plot_figure,
    tap_tool,
    groups,
    layout,
    coordinate_start,
    coordinate_end,
    arrow_source=None,
    read_label_source=None,
    intron_source=None,
    intron_arrow_source=None,
    feature_label_source=None,
    coverage_source=None,
    total_coverage_source=None,
    coverage_figure=None,
):
    """Render clickable IsoSeq transcript rows and return their source."""
    rows = {
        "left": [],
        "right": [],
        "top": [],
        "bottom": [],
        "color": [],
        "alpha": [],
        "base_alpha": [],
        "line_alpha": [],
        "base_line_alpha": [],
        "gene_id": [],
        "gene_name": [],
        "transcript_id": [],
        "transcript_name": [],
        "coordinates": [],
        "span": [],
        "exon_count": [],
        "intron_count": [],
        "strand": [],
        "assigned_read_count": [],
        "gene_transcript_count": [],
        "source_kind": [],
        "chunk_url": [],
        "coverage_url": [],
        "chunk_key": [],
        "transcript_area_height": [],
        "selected_read_y_start": [],
        "selected_view_height": [],
        "transcript_plot_height": [],
        "selected_plot_height": [],
    }
    feature_data = {
        "x": [],
        "y": [],
        "text": [],
        "feature_type": [],
        "feature_number": [],
        "feature_label": [],
        "feature_length": [],
        "gene_id": [],
        "gene_name": [],
        "transcript_id": [],
        "transcript_name": [],
        "coordinates": [],
        "span": [],
        "exon_count": [],
        "intron_count": [],
        "strand": [],
        "assigned_read_count": [],
        "gene_transcript_count": [],
        "label_alpha": [],
    }
    intron_data = {
        "xs": [],
        "ys": [],
        "line_color": [],
        "line_alpha": [],
        "gene_id": [],
        "transcript_id": [],
    }
    intron_arrow_data = {
        "x": [],
        "y": [],
        "angle": [],
        "color": [],
        "alpha": [],
        "gene_id": [],
        "transcript_id": [],
    }
    label_x, label_y, label_text, label_color, label_alpha = [], [], [], [], []
    label_gene_id, label_transcript_id = [], []
    transcript_height = 0.36

    for group in groups:
        y = layout["transcript_y"].get(group["group_id"], 0)
        transcript = group.get("transcript")
        if transcript is None:
            rows["left"].append(coordinate_start)
            rows["right"].append(coordinate_end)
            rows["top"].append(y - transcript_height / 2)
            rows["bottom"].append(y + transcript_height / 2)
            rows["color"].append("#777777")
            rows["alpha"].append(0.18)
            rows["base_alpha"].append(0.18)
            rows["line_alpha"].append(0.6)
            rows["base_line_alpha"].append(0.6)
            rows["gene_id"].append("")
            rows["gene_name"].append("Unassigned")
            rows["transcript_id"].append("UNASSIGNED")
            rows["transcript_name"].append("Unassigned")
            rows["coordinates"].append(f"{coordinate_start:,}-{coordinate_end:,}")
            rows["span"].append(_format_span(coordinate_start, coordinate_end))
            rows["exon_count"].append(0)
            rows["intron_count"].append(0)
            rows["strand"].append("unknown")
            rows["assigned_read_count"].append(group["assigned_read_count"])
            rows["gene_transcript_count"].append(0)
            rows["source_kind"].append("unassigned")
            rows["chunk_url"].append(group.get("chunk_url", ""))
            rows["coverage_url"].append(group.get("coverage_url", ""))
            rows["chunk_key"].append(group.get("chunk_key", ""))
            rows["transcript_area_height"].append(layout["transcript_area_height"])
            rows["selected_read_y_start"].append(layout["selected_read_y_start"])
            rows["selected_view_height"].append(layout["height"])
            rows["transcript_plot_height"].append(
                int(max(3, layout["transcript_area_height"] * 0.6) * 80)
            )
            rows["selected_plot_height"].append(int(max(6, layout["height"] * 0.6) * 80))
            rows["selected_plot_height"][-1] = min(rows["selected_plot_height"][-1], 900)
            label_x.append(coordinate_start + (coordinate_end - coordinate_start) * 0.01)
            label_y.append(y)
            label_text.append(f"Unassigned ({group['assigned_read_count']:,} reads)")
            label_color.append("#666666")
            label_alpha.append(1.0)
            label_gene_id.append("")
            label_transcript_id.append("UNASSIGNED")
            continue
        color = _stable_gene_color(transcript.gene_id)
        exons = sorted(transcript.exons, key=lambda exon: exon[0])
        transcript_coordinates = f"{transcript.chrom}:{transcript.start:,}-{transcript.end:,}"
        transcript_span = _format_span(transcript.start, transcript.end)
        exon_count = len(exons)
        for exon_index, (exon_start, exon_end, _source_exon_number) in enumerate(
            exons,
            start=1,
        ):
            display_exon_number = (
                exon_count - exon_index + 1 if transcript.strand == "-" else exon_index
            )
            exon_left, exon_right = _annotation_interval(exon_start, exon_end)
            vis_start = max(exon_left, coordinate_start)
            vis_end = min(exon_right, coordinate_end + 1)
            if vis_start >= vis_end:
                continue
            rows["left"].append(vis_start)
            rows["right"].append(vis_end)
            rows["top"].append(y - transcript_height / 2)
            rows["bottom"].append(y + transcript_height / 2)
            rows["color"].append(color)
            rows["alpha"].append(0.86)
            rows["base_alpha"].append(0.86)
            rows["line_alpha"].append(1.0)
            rows["base_line_alpha"].append(1.0)
            rows["gene_id"].append(transcript.gene_id)
            rows["gene_name"].append(transcript.gene_name)
            rows["transcript_id"].append(transcript.transcript_id)
            rows["transcript_name"].append(transcript.transcript_name)
            rows["coordinates"].append(transcript_coordinates)
            rows["span"].append(transcript_span)
            rows["exon_count"].append(len(exons))
            rows["intron_count"].append(max(0, len(exons) - 1))
            rows["strand"].append(transcript.strand or "unknown")
            rows["assigned_read_count"].append(group["assigned_read_count"])
            rows["gene_transcript_count"].append(group["gene_transcript_count"])
            rows["source_kind"].append("transcript")
            rows["chunk_url"].append(group.get("chunk_url", ""))
            rows["coverage_url"].append(group.get("coverage_url", ""))
            rows["chunk_key"].append(group.get("chunk_key", ""))
            rows["transcript_area_height"].append(layout["transcript_area_height"])
            rows["selected_read_y_start"].append(layout["selected_read_y_start"])
            rows["selected_view_height"].append(layout["height"])
            rows["transcript_plot_height"].append(
                int(max(3, layout["transcript_area_height"] * 0.6) * 80)
            )
            rows["selected_plot_height"].append(int(max(6, layout["height"] * 0.6) * 80))
            rows["selected_plot_height"][-1] = min(rows["selected_plot_height"][-1], 900)
            feature_data["x"].append((vis_start + vis_end) / 2)
            feature_data["y"].append(y)
            feature_data["text"].append(str(display_exon_number))
            feature_data["feature_type"].append("exon")
            feature_data["feature_number"].append(display_exon_number)
            feature_data["feature_label"].append(f"Exon {display_exon_number}")
            feature_data["feature_length"].append(f"{max(0, exon_end - exon_start + 1):,} bp")
            feature_data["gene_id"].append(transcript.gene_id)
            feature_data["gene_name"].append(transcript.gene_name)
            feature_data["transcript_id"].append(transcript.transcript_id)
            feature_data["transcript_name"].append(transcript.transcript_name)
            feature_data["coordinates"].append(transcript_coordinates)
            feature_data["span"].append(transcript_span)
            feature_data["exon_count"].append(len(exons))
            feature_data["intron_count"].append(max(0, len(exons) - 1))
            feature_data["strand"].append(transcript.strand or "unknown")
            feature_data["assigned_read_count"].append(group["assigned_read_count"])
            feature_data["gene_transcript_count"].append(group["gene_transcript_count"])
            feature_data["label_alpha"].append(0.8)
        for i in range(len(exons) - 1):
            intron_start = max(exons[i][1] + 1, coordinate_start)
            intron_end = min(exons[i + 1][0], coordinate_end)
            if intron_start < intron_end:
                intron_count = len(exons) - 1
                intron_number = intron_count - i if transcript.strand == "-" else i + 1
                intron_data["xs"].append([intron_start, intron_end])
                intron_data["ys"].append([y, y])
                intron_data["line_color"].append(color)
                intron_data["line_alpha"].append(0.75)
                intron_data["gene_id"].append(transcript.gene_id)
                intron_data["transcript_id"].append(transcript.transcript_id)
                feature_data["x"].append((intron_start + intron_end) / 2)
                feature_data["y"].append(y)
                feature_data["text"].append(str(intron_number))
                feature_data["feature_type"].append("intron")
                feature_data["feature_number"].append(intron_number)
                feature_data["feature_label"].append(f"Intron {intron_number}")
                feature_data["feature_length"].append(f"{max(0, intron_end - intron_start):,} bp")
                feature_data["gene_id"].append(transcript.gene_id)
                feature_data["gene_name"].append(transcript.gene_name)
                feature_data["transcript_id"].append(transcript.transcript_id)
                feature_data["transcript_name"].append(transcript.transcript_name)
                feature_data["coordinates"].append(transcript_coordinates)
                feature_data["span"].append(transcript_span)
                feature_data["exon_count"].append(len(exons))
                feature_data["intron_count"].append(max(0, len(exons) - 1))
                feature_data["strand"].append(transcript.strand or "unknown")
                feature_data["assigned_read_count"].append(group["assigned_read_count"])
                feature_data["gene_transcript_count"].append(group["gene_transcript_count"])
                feature_data["label_alpha"].append(0.8)
                for marker_x in _intron_direction_positions(
                    intron_start,
                    intron_end,
                    coordinate_end - coordinate_start,
                ):
                    intron_arrow_data["x"].append(marker_x)
                    intron_arrow_data["y"].append(y)
                    intron_arrow_data["angle"].append(
                        -math.pi / 2 if transcript.strand != "-" else math.pi / 2
                    )
                    intron_arrow_data["color"].append(color)
                    intron_arrow_data["alpha"].append(0.88)
                    intron_arrow_data["gene_id"].append(transcript.gene_id)
                    intron_arrow_data["transcript_id"].append(transcript.transcript_id)
        transcript_left, _transcript_right = _annotation_interval(
            transcript.start,
            transcript.end,
        )
        label_x.append(max(transcript_left, coordinate_start))
        label_y.append(y - transcript_height / 2)
        label_text.append(
            f"{transcript.transcript_name or transcript.transcript_id} "
            f"({group['assigned_read_count']:,} reads)"
        )
        label_color.append(color)
        label_alpha.append(1.0)
        label_gene_id.append(transcript.gene_id)
        label_transcript_id.append(transcript.transcript_id)

    intron_source = None
    if intron_data["xs"]:
        intron_source = ColumnDataSource(data=intron_data)
        plot_figure.multi_line(
            xs="xs",
            ys="ys",
            source=intron_source,
            line_color="line_color",
            line_width=1.25,
            line_alpha="line_alpha",
        )
    if intron_arrow_data["x"]:
        intron_arrow_source = ColumnDataSource(data=intron_arrow_data)
        plot_figure.scatter(
            x="x",
            y="y",
            source=intron_arrow_source,
            marker="triangle",
            size=9,
            angle="angle",
            fill_color="color",
            fill_alpha="alpha",
            line_color="color",
            line_alpha="alpha",
            selection_fill_alpha="alpha",
            selection_line_alpha="alpha",
            nonselection_fill_alpha="alpha",
            nonselection_line_alpha="alpha",
        )

    transcript_source = None
    if rows["left"]:
        transcript_source = ColumnDataSource(data=rows)
        renderer = plot_figure.quad(
            left="left",
            right="right",
            top="top",
            bottom="bottom",
            source=transcript_source,
            fill_color="color",
            fill_alpha="alpha",
            line_color="color",
            line_alpha="line_alpha",
            selection_fill_alpha="alpha",
            selection_line_alpha="line_alpha",
            nonselection_fill_alpha="alpha",
            nonselection_line_alpha="line_alpha",
        )
        if tap_tool is None:
            tap_tool = TapTool(renderers=[renderer])
            plot_figure.add_tools(tap_tool)
        else:
            tap_tool.renderers = [*tap_renderers(tap_tool)]

    label_source = None
    if label_x:
        label_source = ColumnDataSource(
            data={
                "x": label_x,
                "y": label_y,
                "text": label_text,
                "color": label_color,
                "alpha": label_alpha,
                "base_alpha": label_alpha,
                "gene_id": label_gene_id,
                "transcript_id": label_transcript_id,
            }
        )
        plot_figure.text(
            x="x",
            y="y",
            text="text",
            source=label_source,
            text_font_size="8pt",
            text_color="color",
            text_alpha="alpha",
            text_font_style="italic",
            text_baseline="bottom",
        )
    if feature_data["x"]:
        feature_label_source = ColumnDataSource(data=feature_data)
        cfg = PLOT_CONFIG
        plot_figure.scatter(
            x="x",
            y="y",
            source=feature_label_source,
            size=cfg["alignment_label_visible_size"],
            marker="circle",
            fill_color=cfg["alignment_label_fill_color"],
            line_color=cfg["alignment_label_fill_color"],
            alpha="label_alpha",
        )
        feature_hit_renderer = plot_figure.scatter(
            x="x",
            y="y",
            source=feature_label_source,
            size=cfg["alignment_label_hit_size"],
            marker="circle",
            fill_color=cfg["alignment_label_fill_color"],
            line_color=cfg["alignment_label_fill_color"],
            alpha=0,
        )
        plot_figure.text(
            x="x",
            y="y",
            text="text",
            source=feature_label_source,
            text_font_size=cfg["alignment_label_text_font_size"],
            text_color=cfg["alignment_label_text_color"],
            text_alpha="label_alpha",
            text_align="center",
            text_baseline="middle",
            text_font_style="bold",
        )
        tap_tool.renderers = [*tap_renderers(tap_tool), feature_hit_renderer]
        feature_label_source.selected.js_on_change(
            "indices",
            CustomJS(
                args={"source": feature_label_source},
                code=load_javascript("transcript_feature_click_callback.js"),
                module=True,
            ),
        )
    if transcript_source is not None and arrow_source is not None and read_label_source is not None:
        transcript_source.selected.js_on_change(
            "indices",
            CustomJS(
                args={
                    "source": transcript_source,
                    "arrow_source": arrow_source,
                    "label_source": read_label_source,
                    "intron_source": intron_source,
                    "intron_arrow_source": intron_arrow_source,
                    "feature_label_source": feature_label_source,
                    "coverage_source": coverage_source,
                    "total_coverage_source": total_coverage_source,
                    "coverage_figure": coverage_figure,
                    "plot_figure": plot_figure,
                },
                code=load_javascript("transcript_click_callback.js"),
                module=True,
            ),
        )
        plot_figure.js_on_event(
            "tap",
            CustomJS(
                args={
                    "source": transcript_source,
                    "intron_source": intron_source,
                },
                code=load_javascript("transcript_select_tap_callback.js"),
                module=True,
            ),
        )
    return transcript_source, label_source, intron_source, intron_arrow_source, feature_label_source


def create_isoseq_controls(
    transcript_sources,
    read_sources,
    transcript_label_sources=None,
    intron_sources=None,
    intron_arrow_sources=None,
    feature_label_sources=None,
    read_arrow_sources=None,
    read_label_sources=None,
    plot_figures=None,
    dotplot_thumbnail=None,
):
    """Create IsoSeq transcript search and gene filter controls."""
    if not transcript_sources:
        return None
    transcript_label_sources = transcript_label_sources or []
    intron_sources = intron_sources or []
    intron_arrow_sources = intron_arrow_sources or []
    feature_label_sources = feature_label_sources or []
    read_arrow_sources = read_arrow_sources or []
    read_label_sources = read_label_sources or []
    plot_figures = plot_figures or []
    control_height = 30
    select_height = 22
    checkbox_margin = (6, 6, 0, 0)
    transcript_entries = []
    seen_transcripts = set()
    for source in transcript_sources:
        data = source.data
        transcript_rows = zip(
            data.get("transcript_id", []),
            data.get("transcript_name", []),
            data.get("source_kind", []),
            strict=False,
        )
        for transcript_id, transcript_name, source_kind in transcript_rows:
            seen_before = transcript_id in seen_transcripts
            if source_kind != "transcript" or not transcript_id or seen_before:
                continue
            seen_transcripts.add(transcript_id)
            transcript_label = transcript_name or transcript_id
            option_label = f"{transcript_label} ({transcript_id})"
            transcript_entries.append((transcript_id, option_label, transcript_label))
    transcript_entries.sort(
        key=lambda entry: (_natural_sort_key(entry[2]), _natural_sort_key(entry[0]))
    )
    transcript_options = [
        ("ALL", "All transcripts"),
        *[
            (transcript_id, option_label)
            for transcript_id, option_label, _label in transcript_entries
        ],
    ]
    transcript_select = Select(
        title="",
        value="ALL",
        options=transcript_options,
        width=245,
        height=select_height,
        margin=(0, 0, 0, 0),
    )
    transcript_select.js_on_change(
        "value",
        CustomJS(
            args={
                "transcript_select": transcript_select,
                "transcript_sources": transcript_sources,
            },
            code=load_javascript("isoseq_transcript_select_callback.js"),
            module=True,
        ),
    )

    gene_options = [("ALL", "All genes")]
    seen = set()
    for source in transcript_sources:
        data = source.data
        for gene_id, gene_name in zip(
            data.get("gene_id", []),
            data.get("gene_name", []),
            strict=False,
        ):
            if not gene_id or gene_id in seen:
                continue
            seen.add(gene_id)
            gene_options.append((gene_id, f"{gene_name} ({gene_id})"))
    gene_select = Select(
        title="",
        value="ALL",
        options=gene_options,
        width=220,
        height=select_height,
        margin=(0, 0, 0, 0),
    )
    gene_select.js_on_change(
        "value",
        CustomJS(
            args={
                "gene_select": gene_select,
                "transcript_sources": transcript_sources,
                "read_sources": read_sources,
            },
            code=load_javascript("isoseq_gene_filter_callback.js"),
            module=True,
        ),
    )
    hide_unselected_isoforms_checkbox = Checkbox(
        label="Hide unselected isoforms",
        active=False,
        width=165,
        height=control_height,
        margin=checkbox_margin,
    )
    hide_unselected_reads_checkbox = Checkbox(
        label="Hide unselected reads",
        active=False,
        width=150,
        height=control_height,
        margin=checkbox_margin,
    )
    visibility_callback = CustomJS(
        args={
            "transcript_sources": transcript_sources,
            "isoseq_label_sources": transcript_label_sources,
            "isoseq_intron_sources": intron_sources,
            "isoseq_intron_arrow_sources": intron_arrow_sources,
            "isoseq_feature_label_sources": feature_label_sources,
            "read_arrow_sources": read_arrow_sources,
            "read_label_sources": read_label_sources,
            "plot_figures": plot_figures,
            "hide_unselected_isoforms_checkbox": hide_unselected_isoforms_checkbox,
            "hide_unselected_reads_checkbox": hide_unselected_reads_checkbox,
        },
        code=load_javascript("isoseq_visibility_controls_callback.js"),
        module=True,
    )
    hide_unselected_isoforms_checkbox.js_on_change("active", visibility_callback)
    hide_unselected_reads_checkbox.js_on_change("active", visibility_callback)
    for source in [*read_arrow_sources, *read_label_sources]:
        source.selected.js_on_change("indices", visibility_callback)
    controls = [
        transcript_select,
        gene_select,
        *([dotplot_thumbnail] if dotplot_thumbnail is not None else []),
        hide_unselected_isoforms_checkbox,
        hide_unselected_reads_checkbox,
    ]
    control_width = sum(getattr(control, "width", 0) or 0 for control in controls)
    control_width += 8 * max(0, len(controls) - 1)
    return row(
        *controls,
        spacing=8,
        sizing_mode="fixed",
        width=control_width,
        height=control_height,
        align="center",
    )


def build_isoseq_bam_row_track(
    bam_row,
    region,
    coordinate_start,
    coordinate_end,
    region_state,
    region_idx,
    row_index,
    alignment_summaries=None,
    add_vcf_track_to_region_fn=None,
    isoseq_chunk_dir=None,
    isoseq_chunk_url_prefix=None,
):
    """Add one IsoSeq BAM row to the region using transcript-grouped layout."""
    segments_by_read = bam_row.get("segments_by_read", {})
    groups = bam_row.get("isoseq_groups", [])
    if not groups and not segments_by_read:
        logger.debug(f"No IsoSeq data to plot for region {region_idx + 1}, row {row_index + 1}.")
        return

    sample_label = bam_row.get("sample_label")
    layout = _build_isoseq_layout(groups, segments_by_read)
    alignments_height = layout["transcript_area_height"]

    if region_state.shared_x_range is None:
        plot_figure, tap_tool = create_bokeh_figure(
            coordinate_start, coordinate_end, alignments_height
        )
        region_state.shared_x_range = plot_figure.x_range
        region_state.first_plot_figure = plot_figure
        region_state.region_type = ISOSEQ_REGION_TYPE
    else:
        plot_figure, tap_tool = create_bokeh_figure_shared_x(
            region_state.shared_x_range, alignments_height
        )
    plot_dom_class = f"orographer-isoseq-plot-r{region_idx}-row{row_index}"
    plot_figure.css_classes = [*plot_figure.css_classes, plot_dom_class]
    total_coverage_data = _prepare_isoseq_lazy_chunks(
        groups,
        segments_by_read,
        layout,
        coordinate_start,
        coordinate_end,
        region_idx,
        row_index,
        sample_label,
        plot_dom_class,
        plot_figure.id,
        chunk_dir=isoseq_chunk_dir,
        chunk_url_prefix=isoseq_chunk_url_prefix,
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

    vcf_figure = add_vcf_track_to_region_fn(
        plot_figure,
        bam_row.get("vcf_variants"),
        coordinate_start,
        coordinate_end,
        sample_label=sample_label,
    )
    if vcf_figure:
        region_state.row_components.append(vcf_figure)
        region_state.cursor_guide_figures.append(vcf_figure)
        region_state.vcf_figures.append(vcf_figure)
    else:
        region_state.vcf_figures.append(None)

    coverage_figure = create_coverage_figure(plot_figure)
    coverage_source = ColumnDataSource(data=total_coverage_data)
    total_coverage_source = ColumnDataSource(data=total_coverage_data)
    update_coverage_y_range(coverage_figure, total_coverage_data["y"])
    coverage_figure.varea(
        x="x",
        y1=0,
        y2="y",
        source=coverage_source,
        color="#8fa8c8",
        alpha=0.72,
    )
    add_coverage_hover_to_figure(coverage_figure, coverage_source)
    region_state.row_components.append(coverage_figure)
    region_state.cursor_guide_figures.append(coverage_figure)

    page_button_styles = {
        "font-size": "14px",
        "font-weight": "bold",
        "line-height": "1",
        "min-width": "0",
        "padding": "0",
        "text-align": "center",
    }
    page_first_button = Button(
        label="\u00ab",
        button_type="default",
        width=32,
        height=22,
        disabled=True,
        margin=(0, 0, 0, 0),
        styles=page_button_styles.copy(),
    )
    page_prev_button = Button(
        label="\u2039",
        button_type="default",
        width=32,
        height=22,
        disabled=True,
        margin=(0, 0, 0, 0),
        styles=page_button_styles.copy(),
    )
    page_next_button = Button(
        label="\u203a",
        button_type="default",
        width=32,
        height=22,
        disabled=True,
        margin=(0, 0, 0, 0),
        styles=page_button_styles.copy(),
    )
    page_last_button = Button(
        label="\u00bb",
        button_type="default",
        width=32,
        height=22,
        disabled=True,
        margin=(0, 0, 0, 0),
        styles=page_button_styles.copy(),
    )
    page_status_div = Div(
        text="",
        width=180,
        height=24,
        margin=(3, 0, 0, 0),
        styles={"font-size": "11px", "color": "#555", "line-height": "18px"},
    )
    region_state.row_components.append(
        row(
            page_first_button,
            page_prev_button,
            page_next_button,
            page_last_button,
            page_status_div,
            spacing=6,
            height=26,
        )
    )

    arrow_result = add_arrows_to_plot(plot_figure, _empty_arrow_source_data(), allow_empty=True)
    arrow_source = arrow_result[0] if arrow_result else None
    arrow_renderer = arrow_result[1] if arrow_result else None
    if arrow_source:
        region_state.arrow_sources.append(arrow_source)
    if arrow_renderer:
        region_state.arrow_renderers.append(arrow_renderer)

    top_padding = max(alignments_height * 0.02, 0.5)
    region_state.y_bounds.append((plot_figure.y_range, alignments_height, -top_padding))

    (
        transcript_source,
        transcript_label_source,
        transcript_intron_source,
        transcript_intron_arrow_source,
        transcript_feature_label_source,
    ) = add_isoseq_transcripts_to_plot(
        plot_figure,
        tap_tool,
        groups,
        layout,
        coordinate_start,
        coordinate_end,
        arrow_source=arrow_source,
        read_label_source=None,
        coverage_source=coverage_source,
        total_coverage_source=total_coverage_source,
        coverage_figure=coverage_figure,
    )
    if transcript_source is not None:
        region_state.transcript_sources.append(transcript_source)
        region_state.isoseq_filter_sources.append(transcript_source)
        region_state.isoseq_transcript_label_sources.append(transcript_label_source)
        region_state.isoseq_intron_sources.append(transcript_intron_source)
        region_state.isoseq_intron_arrow_sources.append(transcript_intron_arrow_source)
        region_state.isoseq_feature_label_sources.append(transcript_feature_label_source)
    if transcript_label_source is not None:
        region_state.isoseq_filter_sources.append(transcript_label_source)
    if transcript_intron_source is not None:
        region_state.isoseq_filter_sources.append(transcript_intron_source)
    if transcript_intron_arrow_source is not None:
        region_state.isoseq_filter_sources.append(transcript_intron_arrow_source)
    if transcript_feature_label_source is not None:
        region_state.isoseq_filter_sources.append(transcript_feature_label_source)

    variant_renderers = add_variants_to_plot(
        plot_figure,
        {
            "mismatch": {
                "x": [],
                "y": [],
                "ref": [],
                "alt": [],
                "color": [],
                "has_split_alignment": [],
                "has_multiregion_connection": [],
                "read_filter_alpha": [],
                "marker_fill_alpha": [],
                "marker_line_alpha": [],
                "text_alpha": [],
            },
            "insertion": {
                "x": [],
                "y": [],
                "size": [],
                "count": [],
                "is_1bp": [],
                "has_split_alignment": [],
                "has_multiregion_connection": [],
                "read_filter_alpha": [],
                "marker_fill_alpha": [],
                "marker_line_alpha": [],
                "text_alpha": [],
            },
            "deletion": {
                "x0": [],
                "x1": [],
                "y": [],
                "is_1bp": [],
                "has_split_alignment": [],
                "has_multiregion_connection": [],
                "read_filter_alpha": [],
                "line_alpha": [],
            },
        },
    )
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
        {"x": [], "y": [], "customdata": []},
        arrow_source,
        arrow_renderer,
        arrow_tap_callback=CustomJS(code=""),
        allow_empty=True,
    )
    if label_result:
        label_source, label_renderers = label_result
        region_state.alignment_label_sources.append(label_source)
        region_state.alignment_label_renderers.extend(label_renderers)
        region_state.isoseq_filter_sources.append(label_source)
        if transcript_source is not None:
            transcript_source.selected.js_on_change(
                "indices",
                CustomJS(
                    args={
                        "source": transcript_source,
                        "arrow_source": arrow_source,
                        "label_source": label_source,
                        "transcript_label_source": transcript_label_source,
                        "intron_source": transcript_intron_source,
                        "intron_arrow_source": transcript_intron_arrow_source,
                        "feature_label_source": transcript_feature_label_source,
                        "coverage_source": coverage_source,
                        "total_coverage_source": total_coverage_source,
                        "coverage_figure": coverage_figure,
                        "plot_figure": plot_figure,
                        "page_status_div": page_status_div,
                        "page_first_button": page_first_button,
                        "page_prev_button": page_prev_button,
                        "page_next_button": page_next_button,
                        "page_last_button": page_last_button,
                    },
                    code=load_javascript("transcript_click_callback.js"),
                    module=True,
                ),
            )
            page_first_button.js_on_click(
                CustomJS(
                    args={
                        "source": transcript_source,
                        "arrow_source": arrow_source,
                        "label_source": label_source,
                        "transcript_label_source": transcript_label_source,
                        "intron_source": transcript_intron_source,
                        "intron_arrow_source": transcript_intron_arrow_source,
                        "feature_label_source": transcript_feature_label_source,
                        "coverage_source": coverage_source,
                        "total_coverage_source": total_coverage_source,
                        "coverage_figure": coverage_figure,
                        "plot_figure": plot_figure,
                        "page_status_div": page_status_div,
                        "page_first_button": page_first_button,
                        "page_prev_button": page_prev_button,
                        "page_next_button": page_next_button,
                        "page_last_button": page_last_button,
                        "page_target": "first",
                    },
                    code=load_javascript("transcript_click_callback.js"),
                    module=True,
                )
            )
            page_prev_button.js_on_click(
                CustomJS(
                    args={
                        "source": transcript_source,
                        "arrow_source": arrow_source,
                        "label_source": label_source,
                        "transcript_label_source": transcript_label_source,
                        "intron_source": transcript_intron_source,
                        "intron_arrow_source": transcript_intron_arrow_source,
                        "feature_label_source": transcript_feature_label_source,
                        "coverage_source": coverage_source,
                        "total_coverage_source": total_coverage_source,
                        "coverage_figure": coverage_figure,
                        "plot_figure": plot_figure,
                        "page_status_div": page_status_div,
                        "page_first_button": page_first_button,
                        "page_prev_button": page_prev_button,
                        "page_next_button": page_next_button,
                        "page_last_button": page_last_button,
                        "page_delta": -1,
                    },
                    code=load_javascript("transcript_click_callback.js"),
                    module=True,
                )
            )
            page_next_button.js_on_click(
                CustomJS(
                    args={
                        "source": transcript_source,
                        "arrow_source": arrow_source,
                        "label_source": label_source,
                        "transcript_label_source": transcript_label_source,
                        "intron_source": transcript_intron_source,
                        "intron_arrow_source": transcript_intron_arrow_source,
                        "feature_label_source": transcript_feature_label_source,
                        "coverage_source": coverage_source,
                        "total_coverage_source": total_coverage_source,
                        "coverage_figure": coverage_figure,
                        "plot_figure": plot_figure,
                        "page_status_div": page_status_div,
                        "page_first_button": page_first_button,
                        "page_prev_button": page_prev_button,
                        "page_next_button": page_next_button,
                        "page_last_button": page_last_button,
                        "page_delta": 1,
                    },
                    code=load_javascript("transcript_click_callback.js"),
                    module=True,
                )
            )
            page_last_button.js_on_click(
                CustomJS(
                    args={
                        "source": transcript_source,
                        "arrow_source": arrow_source,
                        "label_source": label_source,
                        "transcript_label_source": transcript_label_source,
                        "intron_source": transcript_intron_source,
                        "intron_arrow_source": transcript_intron_arrow_source,
                        "feature_label_source": transcript_feature_label_source,
                        "coverage_source": coverage_source,
                        "total_coverage_source": total_coverage_source,
                        "coverage_figure": coverage_figure,
                        "plot_figure": plot_figure,
                        "page_status_div": page_status_div,
                        "page_first_button": page_first_button,
                        "page_prev_button": page_prev_button,
                        "page_next_button": page_next_button,
                        "page_last_button": page_last_button,
                        "page_target": "last",
                    },
                    code=load_javascript("transcript_click_callback.js"),
                    module=True,
                )
            )
            plot_figure.js_on_event(
                "tap",
                CustomJS(
                    args={
                        "source": transcript_source,
                        "intron_source": transcript_intron_source,
                    },
                    code=load_javascript("transcript_select_tap_callback.js"),
                    module=True,
                ),
            )

    bam_row["segments_by_read"] = {}
    for group in groups:
        group["read_names"] = []

    region_state.row_components.append(plot_figure)
    region_state.plot_figures.append(plot_figure)
    region_state.cursor_guide_figures.append(plot_figure)
