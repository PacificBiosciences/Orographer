"""Reusable Bokeh glyph renderers for read, variant, and label tracks."""

import math
from typing import Any

from bokeh.models import ColumnDataSource, CustomJS

from .callbacks import get_arrow_tap_callback
from .data import add_read_filter_visibility_columns
from .sources import add_read_chevron_data, clickable_label_source_data
from .utils import PLOT_CONFIG, load_javascript, tap_renderers


def _initial_chevron_length(plot_figure: Any) -> float | None:
    x_range = getattr(plot_figure, "x_range", None)
    width = getattr(plot_figure, "width", None)
    range_start = getattr(x_range, "start", None)
    range_end = getattr(x_range, "end", None)
    if range_start is None or range_end is None or not width:
        return None
    if not math.isfinite(range_start) or not math.isfinite(range_end):
        return None
    return abs(range_end - range_start) / width * PLOT_CONFIG["read_chevron_target_px"]


def _initial_chevron_y_offset(plot_figure: Any) -> float | None:
    y_range = getattr(plot_figure, "y_range", None)
    height = getattr(plot_figure, "height", None)
    range_start = getattr(y_range, "start", None)
    range_end = getattr(y_range, "end", None)
    if range_start is None or range_end is None or not height:
        return None
    if not math.isfinite(range_start) or not math.isfinite(range_end):
        return None
    return abs(range_end - range_start) / height * PLOT_CONFIG["read_chevron_target_y_px"]


def add_arrows_to_plot(
    plot_figure: Any,
    arrow_data: dict,
    allow_empty: bool = False,
) -> tuple[Any, Any] | None:
    """Add arrow segments and arrowheads to the plot using batched glyphs."""
    if not arrow_data["x0"] and not allow_empty:
        return None
    arrow_data = dict(arrow_data)
    row_count = len(arrow_data["x0"])
    arrow_data.setdefault("y0", arrow_data["y"])
    arrow_data.setdefault("y1", arrow_data["y"])
    arrow_data.setdefault("read_filter_alpha", [1.0] * row_count)
    arrow_data.setdefault("has_split_alignment", [False] * row_count)
    arrow_data.setdefault("has_multiregion_connection", [False] * row_count)
    arrow_data.setdefault(
        "chevron_tip_fraction",
        [PLOT_CONFIG["read_chevron_tip_fraction"]] * row_count,
    )
    arrow_data.setdefault("arrow_line_alpha", [PLOT_CONFIG["arrow_line_alpha"]] * row_count)
    arrow_data.setdefault("arrowhead_alpha", [PLOT_CONFIG["arrowhead_fill_alpha"]] * row_count)
    arrow_data.setdefault("arrow_selected_alpha", [1.0] * row_count)
    arrow_data.setdefault(
        "arrow_nonselected_alpha",
        [PLOT_CONFIG["arrow_nonselection_line_alpha"]] * row_count,
    )
    arrow_data["angle"] = [
        -math.pi / 2 if x1 > x0 else math.pi / 2
        for x0, x1 in zip(arrow_data["x0"], arrow_data["x1"], strict=True)
    ]
    add_read_chevron_data(
        arrow_data,
        _initial_chevron_length(plot_figure),
        _initial_chevron_y_offset(plot_figure),
    )
    add_read_filter_visibility_columns(arrow_data)
    arrow_source = ColumnDataSource(data=arrow_data)
    arrow_renderer = plot_figure.segment(
        x0="x0",
        y0="y",
        x1="x1",
        y1="y",
        source=arrow_source,
        line_color="color",
        line_width=PLOT_CONFIG["arrow_line_width"],
        line_alpha="arrow_line_alpha",
        selection_line_color=PLOT_CONFIG["arrow_selection_line_color"],
        selection_line_alpha="arrow_selected_alpha",
        selection_line_width=PLOT_CONFIG["arrow_selection_line_width"],
        nonselection_line_alpha="arrow_nonselected_alpha",
        nonselection_line_width=PLOT_CONFIG["arrow_nonselection_line_width"],
    )
    plot_figure.multi_line(
        xs="chevron_xs",
        ys="chevron_ys",
        source=arrow_source,
        line_color="color",
        line_alpha="arrowhead_alpha",
        line_width=PLOT_CONFIG["arrow_line_width"],
        line_cap="butt",
        line_join="bevel",
        selection_line_color=PLOT_CONFIG["arrow_selection_line_color"],
        selection_line_alpha="arrow_selected_alpha",
        selection_line_width=PLOT_CONFIG["arrow_selection_line_width"],
        selection_line_cap="butt",
        selection_line_join="bevel",
        nonselection_line_alpha="arrow_nonselected_alpha",
        nonselection_line_width=PLOT_CONFIG["arrow_nonselection_line_width"],
        nonselection_line_cap="butt",
        nonselection_line_join="bevel",
    )
    chevron_callback = CustomJS(
        args={
            "source": arrow_source,
            "plot_figure": plot_figure,
            "x_range": plot_figure.x_range,
            "y_range": plot_figure.y_range,
            "target_px": PLOT_CONFIG["read_chevron_target_px"],
            "target_y_px": PLOT_CONFIG["read_chevron_target_y_px"],
            "tip_fraction": PLOT_CONFIG["read_chevron_tip_fraction"],
        },
        code=load_javascript("read_chevron_callback.js"),
    )
    plot_figure.x_range.js_on_change("start", chevron_callback)
    plot_figure.x_range.js_on_change("end", chevron_callback)
    plot_figure.y_range.js_on_change("start", chevron_callback)
    plot_figure.y_range.js_on_change("end", chevron_callback)
    return arrow_source, arrow_renderer


def _variant_defaults(source_data: dict, row_count: int, defaults: dict) -> dict:
    data = dict(source_data)
    data.setdefault("read_filter_alpha", [1.0] * row_count)
    data.setdefault("has_split_alignment", [False] * row_count)
    data.setdefault("has_multiregion_connection", [False] * row_count)
    for key, value in defaults.items():
        data.setdefault(key, [value] * row_count)
    add_read_filter_visibility_columns(data)
    return data


def _add_insertion_markers(plot_figure: Any, insertion_data: dict, renderers: dict) -> None:
    row_count = len(insertion_data["x"])
    insertion_data = _variant_defaults(
        insertion_data,
        row_count,
        {
            "marker_fill_alpha": PLOT_CONFIG["insertion_fill_alpha"],
            "marker_line_alpha": PLOT_CONFIG["insertion_line_alpha"],
            "text_alpha": 1.0,
        },
    )
    insertion_color = PLOT_CONFIG["variant_color_insertion"]
    idx_1bp = [idx for idx, is_one_bp in enumerate(insertion_data["is_1bp"]) if is_one_bp]
    idx_other = [idx for idx, is_one_bp in enumerate(insertion_data["is_1bp"]) if not is_one_bp]
    for idx_list, is_one_bp in [(idx_1bp, True), (idx_other, False)]:
        if not idx_list:
            continue
        sub = {key: [value[idx] for idx in idx_list] for key, value in insertion_data.items()}
        sub.pop("is_1bp", None)
        sub["text"] = [f"{count}I" for count in sub["count"]]
        add_read_filter_visibility_columns(sub)
        sub_source = ColumnDataSource(data=sub)
        renderers["sources"].append(sub_source)
        marker = plot_figure.scatter(
            x="x",
            y="y",
            source=sub_source,
            marker="square",
            size="size",
            fill_color=insertion_color,
            fill_alpha="marker_fill_alpha",
            line_color=insertion_color,
            line_alpha="marker_line_alpha",
            line_width=PLOT_CONFIG["insertion_line_width"],
        )
        renderers["marker"].append(marker)
        if is_one_bp:
            renderers["one_bp"].append(marker)
            renderers["one_bp_markers"].append(marker)
        text = plot_figure.text(
            x="x",
            y="y",
            text="text",
            source=sub_source,
            text_font_size=PLOT_CONFIG["insertion_text_font_size"],
            text_color=insertion_color,
            text_alpha="text_alpha",
            text_align="center",
            text_baseline="middle",
            text_font_style="bold",
        )
        text.visible = False
        renderers["text"].append(text)
        if is_one_bp:
            renderers["one_bp"].append(text)
            renderers["one_bp_texts"].append(text)


def _add_deletion_markers(plot_figure: Any, deletion_data: dict, renderers: dict) -> None:
    row_count = len(deletion_data["x0"])
    deletion_data = _variant_defaults(
        deletion_data,
        row_count,
        {"line_alpha": PLOT_CONFIG["deletion_line_alpha"]},
    )
    idx_1bp = [idx for idx, is_one_bp in enumerate(deletion_data["is_1bp"]) if is_one_bp]
    idx_other = [idx for idx, is_one_bp in enumerate(deletion_data["is_1bp"]) if not is_one_bp]
    for idx_list, is_one_bp in [(idx_1bp, True), (idx_other, False)]:
        if not idx_list:
            continue
        sub = {
            key: [deletion_data[key][idx] for idx in idx_list]
            for key in (
                "x0",
                "x1",
                "y",
                "line_alpha",
                "read_filter_alpha",
                "has_split_alignment",
                "has_multiregion_connection",
                "read_name",
                "layout_read_name",
            )
            if key in deletion_data
        }
        for key, value in deletion_data.items():
            if key.startswith("y_"):
                sub[key] = [value[idx] for idx in idx_list]
        add_read_filter_visibility_columns(sub)
        sub_source = ColumnDataSource(data=sub)
        renderers["sources"].append(sub_source)
        segment = plot_figure.segment(
            x0="x0",
            y0="y",
            x1="x1",
            y1="y",
            source=sub_source,
            line_color=PLOT_CONFIG["variant_color_deletion"],
            line_width=PLOT_CONFIG["deletion_line_width"],
            line_alpha="line_alpha",
        )
        if is_one_bp:
            renderers["one_bp"].append(segment)
            renderers["one_bp_segments"].append(segment)


def _build_variant_sources_data_list(variant_data: dict) -> list[dict]:
    """Return a list of source dicts parallel to the sources added by add_variants_to_plot.

    Used by _serialize_region_for_swap to produce variant source data for non-built regions.
    Order matches how add_variants_to_plot appends to renderers["sources"].
    """
    result: list[dict] = []
    mismatch_data = dict(variant_data.get("mismatch", {}))
    if mismatch_data.get("x"):
        row_count = len(mismatch_data["x"])
        result.append(
            _variant_defaults(
                mismatch_data,
                row_count,
                {
                    "marker_fill_alpha": PLOT_CONFIG["mismatch_fill_alpha"],
                    "marker_line_alpha": PLOT_CONFIG["mismatch_line_alpha"],
                    "text_alpha": 1.0,
                },
            )
        )
    insertion_data = dict(variant_data.get("insertion", {}))
    if insertion_data.get("x"):
        row_count = len(insertion_data["x"])
        ins_with_defaults = _variant_defaults(
            insertion_data,
            row_count,
            {
                "marker_fill_alpha": PLOT_CONFIG["insertion_fill_alpha"],
                "marker_line_alpha": PLOT_CONFIG["insertion_line_alpha"],
                "text_alpha": 1.0,
            },
        )
        idx_1bp = [idx for idx, flag in enumerate(ins_with_defaults["is_1bp"]) if flag]
        idx_other = [idx for idx, flag in enumerate(ins_with_defaults["is_1bp"]) if not flag]
        for idx_list in (idx_1bp, idx_other):
            if idx_list:
                sub = {k: [v[i] for i in idx_list] for k, v in ins_with_defaults.items()}
                sub.pop("is_1bp", None)
                sub["text"] = [f"{count}I" for count in sub["count"]]
                result.append(sub)
    deletion_data = dict(variant_data.get("deletion", {}))
    if deletion_data.get("x0"):
        row_count = len(deletion_data["x0"])
        del_with_defaults = _variant_defaults(
            deletion_data,
            row_count,
            {"line_alpha": PLOT_CONFIG["deletion_line_alpha"]},
        )
        idx_1bp = [idx for idx, flag in enumerate(del_with_defaults["is_1bp"]) if flag]
        idx_other = [idx for idx, flag in enumerate(del_with_defaults["is_1bp"]) if not flag]
        kept_keys = (
            "x0", "x1", "y", "line_alpha", "read_filter_alpha",
            "has_split_alignment", "has_multiregion_connection", "read_name", "layout_read_name",
        )
        for idx_list in (idx_1bp, idx_other):
            if idx_list:
                sub = {
                    k: [del_with_defaults[k][i] for i in idx_list]
                    for k in kept_keys
                    if k in del_with_defaults
                }
                for key, value in del_with_defaults.items():
                    if key.startswith("y_"):
                        sub[key] = [value[i] for i in idx_list]
                result.append(sub)
    return result


def add_variants_to_plot(plot_figure: Any, variant_data: dict) -> dict:
    """Add variant markers and return renderers grouped for global controls."""
    renderers = {
        "marker": [],
        "text": [],
        "sources": [],
        "one_bp": [],
        "one_bp_markers": [],
        "one_bp_texts": [],
        "one_bp_segments": [],
    }
    mismatch_data = variant_data["mismatch"]
    if mismatch_data["x"]:
        row_count = len(mismatch_data["x"])
        mismatch_data = _variant_defaults(
            mismatch_data,
            row_count,
            {
                "marker_fill_alpha": PLOT_CONFIG["mismatch_fill_alpha"],
                "marker_line_alpha": PLOT_CONFIG["mismatch_line_alpha"],
                "text_alpha": 1.0,
            },
        )
        mismatch_source = ColumnDataSource(data=mismatch_data)
        renderers["sources"].append(mismatch_source)
        mismatch_marker = plot_figure.scatter(
            x="x",
            y="y",
            source=mismatch_source,
            marker="square",
            size=PLOT_CONFIG["mismatch_size"],
            fill_color="color",
            fill_alpha="marker_fill_alpha",
            line_color="color",
            line_alpha="marker_line_alpha",
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
            text_alpha="text_alpha",
            text_align="center",
            text_baseline="middle",
            text_font_style="bold",
        )
        mismatch_text.visible = False
        renderers["text"].append(mismatch_text)
    if variant_data["insertion"]["x"]:
        _add_insertion_markers(plot_figure, variant_data["insertion"], renderers)
    if variant_data["deletion"]["x0"]:
        _add_deletion_markers(plot_figure, variant_data["deletion"], renderers)
    return renderers


def add_clickable_labels(
    plot_figure: Any,
    tap_tool: Any,
    clickable_data: dict,
    arrow_source: Any = None,
    arrow_renderer: Any = None,
    arrow_tap_callback: Any = None,
    reset_callback: Any = None,
    allow_empty: bool = False,
) -> tuple[Any, list[Any]] | None:
    """Add clickable number labels with modal callbacks."""
    if not clickable_data["x"] and not allow_empty:
        return None
    source = ColumnDataSource(data=clickable_label_source_data(clickable_data))
    cfg = PLOT_CONFIG
    visible_circles = plot_figure.scatter(
        "x",
        "y",
        source=source,
        size=cfg["alignment_label_visible_size"],
        marker="circle",
        fill_color=cfg["alignment_label_fill_color"],
        line_color=cfg["alignment_label_fill_color"],
        alpha="label_alpha",
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
    tap_tool.renderers = [*tap_renderers(tap_tool), circles]
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
        text_alpha="label_alpha",
        text_align="center",
        text_baseline="middle",
        text_font_style="bold",
    )
    return source, [visible_circles, circles, label_text]
