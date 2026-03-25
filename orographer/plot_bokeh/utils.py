"""Utilities used only by the plot_bokeh submodule (variant colors, JS loader)."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class RegionBuildState:
    """Mutable state for building one region's layout (one or more BAM rows)."""

    row_components: list[Any] = field(default_factory=list)
    plot_figures: list[Any] = field(default_factory=list)
    vcf_figures: list[Any] = field(default_factory=list)
    arrow_sources: list[Any] = field(default_factory=list)
    arrow_renderers: list[Any] = field(default_factory=list)
    y_bounds: list[Any] = field(default_factory=list)
    one_bp_renderers: list[Any] = field(default_factory=list)
    one_bp_markers: list[Any] = field(default_factory=list)
    one_bp_texts: list[Any] = field(default_factory=list)
    one_bp_segments: list[Any] = field(default_factory=list)
    alignment_label_renderers: list[Any] = field(default_factory=list)
    variant_renderers: list[Any] = field(default_factory=list)
    one_bp_by_row: list[list[Any]] = field(default_factory=list)
    shared_x_range: Any | None = None
    first_plot_figure: Any | None = None
    region_type: str | None = None


# Variant plot colors (IGV-style)
VARIANT_COLOR_A = "#00CC00"  # Green
VARIANT_COLOR_T = "#CC0000"  # Dark red
VARIANT_COLOR_G = "#FFB300"  # Yellow-orange
VARIANT_COLOR_C = "#0000CC"  # Blue
VARIANT_COLOR_INSERTION = "#333333"  # Dark grey
VARIANT_COLOR_DELETION = "#E8E8E8"  # Lighter grey for deletions
VARIANT_COLOR_UNKNOWN = "#888888"  # Grey

# Plot configuration: colors (using constants above where applicable) and sizes.
PLOT_CONFIG = {
    "segment_paraphase_color": "grey",
    "segment_unassigned_color": "lightgrey",
    "segment_fwd_color": "red",
    "segment_rev_color": "blue",
    "base_colors": {
        "A": VARIANT_COLOR_A,
        "T": VARIANT_COLOR_T,
        "G": VARIANT_COLOR_G,
        "C": VARIANT_COLOR_C,
    },
    "variant_color_insertion": VARIANT_COLOR_INSERTION,
    "variant_color_deletion": VARIANT_COLOR_DELETION,
    "variant_color_unknown": VARIANT_COLOR_UNKNOWN,
    "arrow_line_width": 5,
    "arrow_line_alpha": 0.5,
    "arrow_selection_line_color": "black",
    "arrow_selection_line_width": 5,
    "arrow_nonselection_line_alpha": 0.5,
    "arrow_nonselection_line_width": 5,
    "arrowhead_size": 10,
    "arrowhead_fill_alpha": 0.35,
    "arrowhead_line_alpha": 0.35,
    "arrowhead_line_width": 0,
    "mismatch_size": 2.5,
    "mismatch_line_width": 0.25,
    "mismatch_fill_alpha": 0.8,
    "mismatch_line_alpha": 0.9,
    "mismatch_text_font_size": "8pt",
    "insertion_size_min": 2,
    "insertion_size_max": 6,
    "insertion_line_width": 0.25,
    "insertion_fill_alpha": 0.8,
    "insertion_line_alpha": 0.9,
    "insertion_text_font_size": "10pt",
    "deletion_line_width": 3,
    "deletion_line_alpha": 1.0,
    "vcf_triangle_size": 12,
    "vcf_triangle_line_width": 1.5,
    "vcf_triangle_fill_alpha": 1.0,
    "vcf_triangle_line_alpha": 1.0,
    "alignment_label_visible_size": 12,
    "alignment_label_hit_size": 20,
    "alignment_label_fill_color": "white",
    "alignment_label_text_font_size": "8pt",
    "alignment_label_text_color": "black",
    "sample_label_font_size": "16px",
    "sample_label_padding_bottom": "4px",
    "sample_label_color": "#df1995",
}

# Header layout: offset (px) used to center the Go + coordinate input. Injected into
# embed_replace.js at build time. Decrease to shift input left, increase to shift right.
COORD_INPUT_CENTER_OFFSET_PX = 20


def load_javascript(name: str, replacements: dict[str, str] | None = None) -> str:
    """Load a JS file from plot_bokeh/static/js/ and optionally apply replacements.
    Use for Bokeh CustomJS code and injected scripts.
    """
    base = Path(__file__).resolve().parent
    path = base / "static" / "js" / name
    text = path.read_text(encoding="utf-8")
    if replacements:
        for placeholder, value in replacements.items():
            text = text.replace(placeholder, value)
    return text
