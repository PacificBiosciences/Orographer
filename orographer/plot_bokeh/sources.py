"""ColumnDataSource schema builders shared by plot tracks and lazy chunks."""

from .data import add_read_filter_visibility_columns
from .utils import PLOT_CONFIG


def read_chevron_vertices(
    x0: float,
    x1: float,
    y: float,
    chevron_length: float | None = None,
    chevron_y_offset: float | None = None,
    tip_fraction: float | None = None,
) -> tuple[list[float], list[float]]:
    """Return one read-direction chevron polyline kept inside the read span."""
    span = abs(x1 - x0)
    if span == 0:
        return [], []
    target_length = (
        chevron_length
        if chevron_length is not None
        else span * PLOT_CONFIG["read_chevron_fallback_length_fraction"]
    )
    length = min(max(1.0, target_length), span / 2)
    y_offset = (
        chevron_y_offset if chevron_y_offset is not None else PLOT_CONFIG["read_chevron_y_offset"]
    )
    direction = 1 if x1 > x0 else -1
    resolved_tip_fraction = tip_fraction or PLOT_CONFIG["read_chevron_tip_fraction"]
    end_tip = x0 + direction * span * resolved_tip_fraction
    end_base = end_tip - direction * length
    xs = [
        end_base,
        end_tip,
        end_base,
    ]
    ys = [
        y + y_offset,
        y,
        y - y_offset,
    ]
    return xs, ys


def add_read_chevron_data(
    arrow_data: dict,
    chevron_length: float | None = None,
    chevron_y_offset: float | None = None,
) -> None:
    """Populate chevron columns when absent or out of sync with read rows."""
    row_count = len(arrow_data["x0"])
    has_chevrons = (
        len(arrow_data.get("chevron_xs", [])) == row_count
        and len(arrow_data.get("chevron_ys", [])) == row_count
    )
    if has_chevrons and chevron_length is None and chevron_y_offset is None:
        return
    chevron_rows = [
        read_chevron_vertices(
            x0,
            x1,
            y,
            chevron_length,
            chevron_y_offset,
            tip_fraction,
        )
        for x0, x1, y, tip_fraction in zip(
            arrow_data["x0"],
            arrow_data["x1"],
            arrow_data["y"],
            arrow_data.get(
                "chevron_tip_fraction",
                [PLOT_CONFIG["read_chevron_tip_fraction"]] * row_count,
            ),
            strict=True,
        )
    ]
    arrow_data["chevron_xs"] = [xs for xs, _ys in chevron_rows]
    arrow_data["chevron_ys"] = [ys for _xs, ys in chevron_rows]


def clickable_label_source_data(clickable_data: dict) -> dict:
    """Return ColumnDataSource data for clickable alignment labels."""
    source = {
        "x": clickable_data["x"],
        "y": clickable_data["y"],
        "text": [str(info["alignment_number"]) for info in clickable_data["customdata"]],
        "read_name": [info["read_name"] for info in clickable_data["customdata"]],
        "layout_read_name": [
            info.get("layout_read_name", info["read_name"]) for info in clickable_data["customdata"]
        ],
        "alignment_number": [info["alignment_number"] for info in clickable_data["customdata"]],
        "strand": [info["strand"] for info in clickable_data["customdata"]],
        "coordinates": [info["coordinates"] for info in clickable_data["customdata"]],
        "haplotype": [info["haplotype"] for info in clickable_data["customdata"]],
        "sample_label": [info.get("sample_label", "") for info in clickable_data["customdata"]],
        "inclusion_reason": [
            info.get("inclusion_reason", "") for info in clickable_data["customdata"]
        ],
        "all_alignment_numbers": [
            info.get("all_alignment_numbers", []) for info in clickable_data["customdata"]
        ],
        "all_alignment_coordinates": [
            info.get("all_alignment_coordinates", []) for info in clickable_data["customdata"]
        ],
        "has_split_alignment": [
            info.get("has_split_alignment", False) for info in clickable_data["customdata"]
        ],
        "has_multiregion_connection": [
            info.get("has_multiregion_connection", False) for info in clickable_data["customdata"]
        ],
        "gene_id": [info.get("gene_id", "") for info in clickable_data["customdata"]],
        "gene_name": [info.get("gene_name", "") for info in clickable_data["customdata"]],
        "transcript_id": [info.get("transcript_id", "") for info in clickable_data["customdata"]],
        "annotation_label": [
            info.get("annotation_label", "") for info in clickable_data["customdata"]
        ],
        "read_filter_alpha": [1.0 for _info in clickable_data["customdata"]],
        "label_alpha": [0.8 for _info in clickable_data["customdata"]],
    }
    for key, value in clickable_data.items():
        if key.startswith("y_"):
            source[key] = value
    add_read_filter_visibility_columns(source)
    return source


def empty_clickable_label_source_data() -> dict:
    """Return empty label source columns matching clickable_label_source_data()."""
    return clickable_label_source_data({"x": [], "y": [], "customdata": []})


def empty_arrow_source_data() -> dict:
    """Return empty read-glyph columns compatible with add_arrows_to_plot()."""
    return {
        "x0": [],
        "x1": [],
        "y": [],
        "y0": [],
        "y1": [],
        "color": [],
        "read_name": [],
        "layout_read_name": [],
        "segment_id": [],
        "region_idx": [],
        "row_index": [],
        "alignment_order": [],
        "fwd_read_start": [],
        "fwd_read_end": [],
        "source_kind": [],
        "has_split_alignment": [],
        "has_multiregion_connection": [],
        "read_filter_alpha": [],
        "arrow_line_alpha": [],
        "arrowhead_alpha": [],
        "arrow_selected_alpha": [],
        "arrow_nonselected_alpha": [],
        "gene_id": [],
        "gene_name": [],
        "transcript_id": [],
        "group_id": [],
        "angle": [],
        "chevron_tip_fraction": [],
        "chevron_xs": [],
        "chevron_ys": [],
        "read_filter_visible_all": [],
        "read_filter_visible_split": [],
        "read_filter_visible_multiregion": [],
        "read_filter_visible_split_multiregion": [],
    }


def empty_isoseq_coverage_data() -> dict:
    """Return empty IsoSeq coverage source columns."""
    return {"x": [], "y": [], "coverage_alpha": []}
