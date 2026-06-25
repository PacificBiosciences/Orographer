"""Utilities used only by the plot_bokeh submodule (variant colors, JS loader)."""

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class RegionBuildState:
    """Mutable state for building one region's layout (one or more BAM rows)."""

    row_components: list[Any] = field(default_factory=list)
    plot_figures: list[Any] = field(default_factory=list)
    cursor_guide_figures: list[Any] = field(default_factory=list)
    vcf_figures: list[Any] = field(default_factory=list)
    arrow_sources: list[Any] = field(default_factory=list)
    arrow_renderers: list[Any] = field(default_factory=list)
    connector_sources: list[Any] = field(default_factory=list)
    read_filter_sources: list[Any] = field(default_factory=list)
    y_bounds: list[Any] = field(default_factory=list)
    read_filter_y_bounds: list[Any] = field(default_factory=list)
    one_bp_renderers: list[Any] = field(default_factory=list)
    one_bp_markers: list[Any] = field(default_factory=list)
    one_bp_texts: list[Any] = field(default_factory=list)
    one_bp_segments: list[Any] = field(default_factory=list)
    alignment_label_renderers: list[Any] = field(default_factory=list)
    alignment_label_sources: list[Any] = field(default_factory=list)
    phase_set_marker_renderers: list[Any] = field(default_factory=list)
    transcript_sources: list[Any] = field(default_factory=list)
    isoseq_filter_sources: list[Any] = field(default_factory=list)
    isoseq_transcript_label_sources: list[Any] = field(default_factory=list)
    isoseq_intron_sources: list[Any] = field(default_factory=list)
    isoseq_intron_arrow_sources: list[Any] = field(default_factory=list)
    isoseq_feature_label_sources: list[Any] = field(default_factory=list)
    isoseq_comparison_components: list[Any] = field(default_factory=list)
    variant_renderers: list[Any] = field(default_factory=list)
    one_bp_by_row: list[list[Any]] = field(default_factory=list)
    shared_x_range: Any | None = None
    first_plot_figure: Any | None = None
    region_type: str | None = None
    # Slot-swap tracking (populated only for complex-SV N>2 region plots)
    coverage_total_sources: list[Any] = field(default_factory=list)
    coverage_hp_sources: list[Any] = field(default_factory=list)
    vcf_sources: list[Any] = field(default_factory=list)
    marker_connector_sources: list[Any] = field(default_factory=list)
    arc_connector_sources: list[Any] = field(default_factory=list)
    clearable_sources: list[Any] = field(default_factory=list)
    slot_dom_classes: list[str] = field(default_factory=list)
    plot_model_ids: list[str] = field(default_factory=list)
    orig_size_div: Any | None = None
    orig_coord_div: Any | None = None
    # Coordinate-navigation Divs (updatable by swap callback so nav bounds track swapped region)
    coord_input: Any | None = None
    nav_chrom_div: Any | None = None
    nav_orig_start_div: Any | None = None
    nav_orig_end_div: Any | None = None
    # Gene track sources for the region (one dict per slot, for swap restoration)
    gene_track_sources: dict = field(default_factory=dict)
    # HP label sources per bam-row: list of {"separator": src|None, "label": src|None}
    hp_label_sources_per_row: list[dict] = field(default_factory=list)
    # Alignment (read number) label source per bam-row (None for rows without labels)
    alignment_label_source_per_row: list[Any] = field(default_factory=list)
    # Separator-line source per bam-row (horizontal dotted lines between reads)
    separator_source_per_row: list[Any] = field(default_factory=list)
    # Phase-set boundary source per bam-row (vertical lines; None if no phase data)
    phase_set_source_per_row: list[Any] = field(default_factory=list)
    # Variant sources per bam-row (variable-length list: mismatch + insertion + deletion)
    variant_sources_per_row: list[list[Any]] = field(default_factory=list)
    # Insertion-summary raw sources per bam-row (one per HP; display sources recluster on x_range)
    insertion_raw_sources_per_row: list[list[Any]] = field(default_factory=list)
    # Repeat-density (Ident) track source and figure (None when no density data)
    repeat_density_source: Any | None = None
    repeat_density_figure: Any | None = None
    # Gene track figure y_range (Range1d); None when region has no gene annotations
    gene_track_y_range: Any | None = None
    # Coverage figure y_range per bam-row (Range1d or None when row has no coverage)
    coverage_y_ranges_per_row: list[Any] = field(default_factory=list)


# Per-haplotype color palette.
# Colorblind-friendly; excludes red/blue reserved for strand direction.
# HP 0 (unassigned) is neutral gray; HP 1+ cycle through the palette.
_HAPLOTYPE_PALETTE = [
    "#E69F00",  # 1  orange
    "#009E73",  # 2  bluish-green
    "#CC79A7",  # 3  reddish purple
    "#D55E00",  # 4  vermillion
    "#56B4E9",  # 5  sky blue
    "#F0E442",  # 6  yellow
    "#0072B2",  # 7  dark blue
    "#8C564B",  # 8  brown
    "#E377C2",  # 9  pink
    "#17BECF",  # 10 teal
]


def get_haplotype_color(hp) -> str:
    """Return a stable color for haplotype *hp*; cycles through palette for HP > 10.

    Accepts int or str HP values; string HPs that are not plain integers use a
    hash-based index so non-numeric paraphase IDs still get consistent colors.
    """
    if hp is None or hp == 0 or str(hp).lower() in ("0", "unassigned", ""):
        return "#888888"
    try:
        idx = int(hp) - 1
    except (TypeError, ValueError):
        digest = hashlib.sha256(str(hp).encode("utf-8")).digest()
        idx = int.from_bytes(digest[:2], byteorder="big") % len(_HAPLOTYPE_PALETTE)
    return _HAPLOTYPE_PALETTE[idx % len(_HAPLOTYPE_PALETTE)]


def build_haplotype_color_map(haplotypes: list[Any]) -> dict[Any, str]:
    """Assign distinct colors for an ordered haplotype set until the palette is exhausted."""
    color_map: dict[Any, str] = {}
    used_colors: set[str] = set()
    named_haplotypes = []

    for haplotype in haplotypes:
        if haplotype in color_map:
            continue
        haplotype_text = str(haplotype).lower()
        if haplotype is None or haplotype == 0 or haplotype_text in {"0", "unassigned", ""}:
            color_map[haplotype] = "#888888"
            used_colors.add(color_map[haplotype])
            continue
        try:
            int(haplotype)
        except (TypeError, ValueError):
            named_haplotypes.append(haplotype)
            continue
        color_map[haplotype] = get_haplotype_color(haplotype)
        used_colors.add(color_map[haplotype])

    available_colors = [color for color in _HAPLOTYPE_PALETTE if color not in used_colors]
    fallback_colors = available_colors or _HAPLOTYPE_PALETTE
    for index, haplotype in enumerate(named_haplotypes):
        color_map[haplotype] = fallback_colors[index % len(fallback_colors)]

    return color_map


def tap_renderers(tap_tool: Any) -> list[Any]:
    """Return concrete renderers from a TapTool, treating Bokeh auto as empty."""
    renderers = getattr(tap_tool, "renderers", None)
    if renderers == "auto" or renderers is None:
        return []
    return list(renderers)


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
    "arrow_line_width": 4,
    "arrow_line_alpha": 0.5,
    "arrow_selection_line_color": "black",
    "arrow_selection_line_width": 4,
    "arrow_nonselection_line_alpha": 0.12,
    "arrow_nonselection_line_width": 4,
    "arrowhead_size": 12,
    "arrowhead_fill_alpha": 0.65,
    "arrowhead_line_alpha": 0.65,
    "arrowhead_line_width": 0,
    "read_chevron_target_px": 14,
    "read_chevron_target_y_px": 5,
    "read_chevron_plot_width_px": 1440,
    "read_chevron_tip_fraction": 0.8,
    "read_chevron_fallback_length_fraction": 0.08,
    "read_chevron_y_offset": 0.035,
    "connector_line_color": "lightgrey",
    "connector_line_width": 1.5,
    "connector_hit_line_width": 14,
    "connector_line_alpha": 0.68,
    "connector_lane_offset": 0.025,
    "connector_selection_line_color": "#92400e",
    "connector_selection_line_alpha": 0.82,
    "connector_selection_line_width": 3.5,
    "connector_marker_size": 11,
    "connector_marker_fill_color": "#f8fafc",
    "connector_marker_line_color": "#4b5563",
    "connector_selection_color": "black",
    "connector_nonselection_alpha": 0.28,
    "connector_endpoint_y_offset": 0.03,
    "phase_block_line_width": 2.0,
    "phase_block_line_alpha": 0.75,
    "phase_block_line_dash": "dotted",
    "cursor_guide_line_color": "#222222",
    "cursor_guide_line_alpha": 0.45,
    "cursor_guide_line_width": 1,
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
    "complex_sv_insertion_marker_target_px": 36,
    "complex_sv_insertion_marker_height": 10,
    "complex_sv_insertion_marker_min_width": 18,
    "complex_sv_insertion_marker_char_px": 7,
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
