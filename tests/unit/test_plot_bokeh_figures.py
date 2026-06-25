from __future__ import annotations

import base64
import unittest.mock
from types import SimpleNamespace

import numpy as np
import pytest
from bokeh.models import (
    BoxZoomTool,
    Button,
    Checkbox,
    ColumnDataSource,
    CustomJSTicker,
    Div,
    GlyphRenderer,
    HoverTool,
    Span,
)
from bokeh.plotting import figure

from orographer.plot_bokeh.figures import (
    _cluster_insertion_sites,
    _insertion_cluster_source_data,
    _insertion_raw_source_data,
    _render_dotplot_png_data_url,
    add_coverage_hover_to_figure,
    add_coverage_to_figure,
    add_cursor_guide_callbacks,
    add_cursor_guide_to_figures,
    add_gene_track,
    create_bokeh_figure,
    create_bokeh_figure_shared_x,
    create_coordinate_display,
    create_coverage_figure,
    create_dotplot_thumbnail,
    create_genomic_x_axis_strip,
    create_global_checkbox_controls,
    create_repeat_density_figure,
    update_coverage_y_range,
)
from orographer.utils import Region

PNG_DATA_URL_PREFIX = "data:image/png;base64,"
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def test_create_coverage_figure_disables_drag_zoom() -> None:
    main_figure = figure(x_range=(100, 200))
    coverage_figure = create_coverage_figure(main_figure)

    assert not coverage_figure.select(type=BoxZoomTool)
    assert coverage_figure.y_range.start == 0
    assert coverage_figure.y_range.end == 1
    assert coverage_figure.yaxis[0].ticker.desired_num_ticks == 4
    assert isinstance(coverage_figure.xaxis[0].ticker, CustomJSTicker)
    assert coverage_figure.xaxis[0].minor_tick_line_color is None


def test_genomic_x_axes_hide_unlabeled_minor_ticks() -> None:
    read_figure, _tap_tool = create_bokeh_figure(100, 200, 4)
    axis_strip = create_genomic_x_axis_strip(read_figure.x_range)

    assert read_figure.xaxis[0].minor_tick_line_color is None
    assert axis_strip.xaxis[0].minor_tick_line_color is None
    assert isinstance(read_figure.xaxis[0].ticker, CustomJSTicker)
    assert read_figure.xaxis[0].ticker.minor_code == "return [];"


def test_shared_x_figure_sets_active_box_zoom_without_select_one() -> None:
    read_figure, _tap_tool = create_bokeh_figure(100, 200, 4)
    shared_figure, _shared_tap_tool = create_bokeh_figure_shared_x(read_figure.x_range, 4)
    box_zoom_tools = shared_figure.select(type=BoxZoomTool)

    assert len(box_zoom_tools) == 1
    assert shared_figure.toolbar.active_drag is box_zoom_tools[0]


def test_gene_track_uses_gtf_inclusive_end_as_interval_boundary() -> None:
    gene = SimpleNamespace(
        gene_id="G1",
        gene_name="Gene1",
        chrom="chr1",
        start=100,
        end=250,
        strand="+",
        exons=[(100, 150, 1), (200, 250, 2)],
    )
    plot_figure = figure(x_range=(90, 260), y_range=(2, 0))

    add_gene_track(plot_figure, [gene], 0, 90, 260)

    quad_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Quad"
    ]
    assert len(quad_sources) == 1
    assert quad_sources[0].data["left"] == [100, 200]
    assert quad_sources[0].data["right"] == [151, 251]
    assert quad_sources[0].data["exon_coordinates"] == [
        "chr1:100-150 (51 bp)",
        "chr1:200-250 (51 bp)",
    ]
    assert quad_sources[0].data["isoform_coordinates"] == [
        "chr1:100-250 (151 bp)",
        "chr1:100-250 (151 bp)",
    ]


def test_gene_track_keeps_first_duplicate_gene_annotation() -> None:
    genes = [
        SimpleNamespace(
            gene_id="G1",
            gene_name="Gene1",
            start=100,
            end=250,
            strand="+",
            exons=[(100, 150, 1)],
        ),
        SimpleNamespace(
            gene_id="G1",
            gene_name="Gene1",
            start=120,
            end=300,
            strand="+",
            exons=[(200, 300, 2)],
        ),
    ]
    plot_figure = figure(x_range=(90, 320), y_range=(2, 0))

    add_gene_track(plot_figure, genes, 0, 90, 320)

    segment_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Segment"
    ]
    quad_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Quad"
    ]
    assert len(segment_sources) == 1
    assert segment_sources[0].data["gene_id"] == ["G1"]
    assert segment_sources[0].data["x0"] == [100]
    assert segment_sources[0].data["x1"] == [251]
    assert quad_sources[0].data["left"] == [100]
    assert quad_sources[0].data["right"] == [151]
    assert quad_sources[0].data["exon_number"] == ["1"]


def test_gene_track_exon_source_includes_representative_selection_method() -> None:
    gene = SimpleNamespace(
        gene_id="G1",
        gene_name="Gene1",
        start=100,
        end=250,
        strand="+",
        exons=[(100, 150, 1), (200, 250, 2)],
        representative_transcript_id="T1",
        representative_transcript_name="Transcript1",
        representative_selection_method="MANE Select",
    )
    plot_figure = figure(x_range=(90, 260), y_range=(2, 0))

    add_gene_track(plot_figure, [gene], 0, 90, 260)

    quad_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Quad"
    ]
    assert quad_sources[0].data["representative_transcript"] == ["T1", "T1"]
    assert quad_sources[0].data["representative_selection_method"] == [
        "MANE Select",
        "MANE Select",
    ]


def test_gene_track_reverse_strand_preserves_representative_exon_numbers() -> None:
    gene = SimpleNamespace(
        gene_id="G1",
        gene_name="Gene1",
        start=100,
        end=350,
        strand="-",
        exons=[(100, 150, 1), (200, 250, 2), (300, 350, 3)],
    )
    plot_figure = figure(x_range=(90, 360), y_range=(2, 0))

    add_gene_track(plot_figure, [gene], 0, 90, 360)

    quad_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Quad"
    ]
    assert quad_sources[0].data["left"] == [100, 200, 300]
    assert quad_sources[0].data["exon_number"] == ["1", "2", "3"]


def test_gene_track_draws_terminal_intron_stub_for_truncated_transcript() -> None:
    gene = SimpleNamespace(
        gene_id="G1",
        gene_name="Gene1",
        start=100,
        end=300,
        strand="+",
        exons=[(100, 150, 1)],
    )
    plot_figure = figure(x_range=(90, 220), y_range=(2, 0))

    add_gene_track(plot_figure, [gene], 0, 90, 220)

    multi_line_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "MultiLine"
    ]
    assert len(multi_line_sources) == 1
    assert multi_line_sources[0].data["xs"] == [[151, 220]]


def test_gene_track_large_regions_start_in_overview_mode() -> None:
    gene = SimpleNamespace(
        gene_id="G1",
        gene_name="Gene1",
        start=100_000,
        end=300_000,
        strand="+",
        exons=[(100_000, 100_100, 1), (299_900, 300_000, 2)],
    )
    plot_figure = figure(x_range=(1, 2_000_000), y_range=(2, 0))

    add_gene_track(plot_figure, [gene], 0, 1, 2_000_000)

    quad_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Quad"
    ]
    segment_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Segment"
    ]
    text_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Text"
    ]
    scatter_sources = [
        renderer.data_source
        for renderer in plot_figure.renderers
        if isinstance(renderer, GlyphRenderer) and renderer.glyph.__class__.__name__ == "Scatter"
    ]
    assert quad_sources[0].data["fill_alpha"] == [0.0, 0.0]
    assert quad_sources[0].data["fill_alpha_detail"] == [0.8, 0.8]
    assert segment_sources[0].data["line_alpha"] == [0.75]
    assert text_sources[0].data["text_alpha"] == [0.0]
    assert text_sources[0].data["text_alpha_medium"] == [0.0]
    assert text_sources[0].data["text_alpha_detail"] == [0.9]
    assert len(scatter_sources[0].data["x"]) == 12
    assert scatter_sources[0].data["fill_alpha"] == [0.0] * 12
    assert scatter_sources[0].data["fill_alpha_detail"] == [0.9] * 12


def test_gene_track_overview_rows_do_not_reserve_all_label_widths() -> None:
    genes = [
        SimpleNamespace(
            gene_id="G1",
            gene_name="LongGeneNameThatWouldOtherwiseReserveTooMuchSpaceOne",
            start=100_000,
            end=101_000,
            strand="+",
            exons=[(100_000, 101_000, 1)],
        ),
        SimpleNamespace(
            gene_id="G2",
            gene_name="LongGeneNameThatWouldOtherwiseReserveTooMuchSpaceTwo",
            start=400_000,
            end=401_000,
            strand="+",
            exons=[(400_000, 401_000, 1)],
        ),
    ]
    plot_figure = figure(x_range=(1, 2_000_000), y_range=(4, 0))

    gene_track_height, _gene_sources = add_gene_track(plot_figure, genes, 0, 1, 2_000_000)

    assert gene_track_height == 2.2


def test_update_coverage_y_range_uses_finite_minimum() -> None:
    main_figure = figure(x_range=(100, 200))
    coverage_figure = create_coverage_figure(main_figure)

    update_coverage_y_range(coverage_figure, [0, 10, 5])

    assert coverage_figure.y_range.start == 0
    assert coverage_figure.y_range.end == 10.5


def test_add_coverage_hover_to_figure_targets_invisible_line_renderer() -> None:
    main_figure = figure(x_range=(100, 200))
    coverage_figure = create_coverage_figure(main_figure)
    source = ColumnDataSource({"x": [100, 150, 200], "y": [0, 12, 0]})

    renderer = add_coverage_hover_to_figure(coverage_figure, source)

    hover_tools = coverage_figure.select(type=HoverTool)
    assert len(hover_tools) == 1
    assert hover_tools[0].renderers == [renderer]
    assert renderer.glyph.line_alpha == 0


def test_add_coverage_to_figure_adds_total_coverage_hover_renderer() -> None:
    main_figure = figure(x_range=(100, 200))
    coverage_figure = create_coverage_figure(main_figure)

    add_coverage_to_figure(coverage_figure, {-1: ([100, 150, 200], [0, 12, 0])})

    hover_tools = coverage_figure.select(type=HoverTool)
    assert len(hover_tools) == 1
    assert len(hover_tools[0].renderers) == 1


def test_add_coverage_to_figure_labels_named_haplotypes() -> None:
    main_figure = figure(x_range=(100, 200))
    coverage_figure = create_coverage_figure(main_figure)

    add_coverage_to_figure(
        coverage_figure,
        {
            -1: ([100, 150, 200], [0, 4, 0]),
            2: ([100, 150, 200], [0, 2, 0]),
            "hba_hba1hap1": ([100, 150, 200], [0, 3, 0]),
            "smn1_smn2hap1": ([100, 150, 200], [0, 2, 0]),
            "smn1_smn2hap2": ([100, 150, 200], [0, 1, 0]),
            "Unassigned": ([100, 150, 200], [0, 1, 0]),
        },
    )

    sources_with_names = [
        renderer.data_source.data
        for renderer in coverage_figure.renderers
        if hasattr(renderer, "data_source") and "names" in renderer.data_source.data
    ]
    hover_tools = coverage_figure.select(type=HoverTool)
    hover_formatter = hover_tools[0].formatters["@y"]

    assert sources_with_names
    assert sources_with_names[0]["names"] == [
        "HP2",
        "hba_hba1hap1",
        "smn1_smn2hap1",
        "smn1_smn2hap2",
        "Unassigned",
    ]
    assert sources_with_names[0]["colors"][2] != sources_with_names[0]["colors"][3]
    assert hover_formatter.args["hp_source"].data["names"] == [
        "HP2",
        "hba_hba1hap1",
        "smn1_smn2hap1",
        "smn1_smn2hap2",
        "Unassigned",
    ]


def _site(
    pos: int,
    *,
    size: int = 313,
    count: int = 1,
    chrom: str = "chr1",
    read_names: list[str] | None = None,
) -> dict:
    names = read_names if read_names is not None else [f"read-{pos}"]
    return {
        "pos": pos,
        "count": count,
        "median_size": size,
        "read_names": names,
        "total_count": count,
        "chrom": chrom,
    }


def _decode_png_url(data_url: str) -> bytes:
    assert data_url.startswith(PNG_DATA_URL_PREFIX)
    return base64.b64decode(data_url.removeprefix(PNG_DATA_URL_PREFIX))


def _collect_checkboxes(model) -> list[Checkbox]:
    checkboxes = []
    children = getattr(model, "children", [])
    for child in children:
        if isinstance(child, Checkbox):
            checkboxes.append(child)
        checkboxes.extend(_collect_checkboxes(child))
    return checkboxes


def test_cluster_insertion_sites_sorts_and_groups_by_distance() -> None:
    clustered = _cluster_insertion_sites(
        [_site(220), _site(100), _site(114), _site(180)],
        cluster_bp=20,
    )

    assert [[site["pos"] for site in cluster] for cluster in clustered] == [
        [100, 114],
        [180],
        [220],
    ]


def test_insertion_cluster_source_data_labels_merged_sites_as_size_and_count() -> None:
    clustered = [
        [_site(100, size=313, count=2), _site(108, size=500, count=3)],
        [_site(220, size=144)],
    ]

    data = _insertion_cluster_source_data(clustered, hp=2, y_pos=1.5)

    assert data["x"] == [104, 220]
    assert data["marker_label"] == ["500(2)", "144"]
    assert data["count"] == [5, 1]
    assert data["median_size"] == [500, 144]
    assert data["hp_label"] == ["HP2", "HP2"]
    assert data["cluster_pos"] == [[100, 108], [220]]
    assert data["cluster_median_size"] == [[313, 500], [144]]


def test_insertion_cluster_source_data_rounds_fractional_weighted_median() -> None:
    clustered = [[_site(100, size=313), _site(108, size=500)]]

    data = _insertion_cluster_source_data(clustered, hp=1, y_pos=1.0)

    assert data["marker_label"] == ["407(2)"]
    assert data["median_size"] == [407]


def test_insertion_raw_source_data_preserves_sites_for_browser_reclustering() -> None:
    sites = [
        _site(100, count=2, read_names=["r1", "r2", "r3", "r4", "r5", "r6"]),
        _site(160, size=42),
    ]

    data = _insertion_raw_source_data(sites, hp=0, y_pos=2.0)

    assert data["pos"] == [100, 160]
    assert data["median_size"] == [313, 42]
    assert data["hp_label"] == ["Unassigned", "Unassigned"]
    assert data["y"] == [2.0, 2.0]
    assert data["top_names"][0] == "r1\nr2\nr3\nr4\nr5"


def test_render_dotplot_png_data_url_returns_valid_png_with_and_without_axes() -> None:
    matrix = np.eye(6, dtype=np.float32)
    region = Region("chr1", 189_100_000, 189_106_000, "chr1:189100000-189106000")

    thumbnail_url = _render_dotplot_png_data_url(
        matrix,
        region,
        with_axes=False,
        figure_size=(1.0, 1.0),
        dpi=64,
    )
    modal_url = _render_dotplot_png_data_url(
        matrix,
        region,
        with_axes=True,
        figure_size=(2.4, 2.4),
        dpi=80,
    )

    assert _decode_png_url(thumbnail_url).startswith(PNG_SIGNATURE)
    assert _decode_png_url(modal_url).startswith(PNG_SIGNATURE)


def test_render_combined_dotplot_png_data_url_returns_valid_png() -> None:
    matrix = np.eye(8, dtype=np.float32)
    region = Region("chr8", 100, 180, "chr8:100-180")
    blocks = [
        {
            "label": "chr8:0.0-0.0 Mb",
            "chromosome": "chr8",
            "start": 0,
            "end": 40,
            "offset_start": 0,
            "offset_end": 40,
        },
        {
            "label": "chr14:0.0-0.0 Mb",
            "chromosome": "chr14",
            "start": 0,
            "end": 40,
            "offset_start": 40,
            "offset_end": 80,
        },
    ]

    modal_url = _render_dotplot_png_data_url(
        matrix,
        region,
        with_axes=True,
        figure_size=(2.4, 2.4),
        dpi=80,
        blocks=blocks,
    )

    assert _decode_png_url(modal_url).startswith(PNG_SIGNATURE)


def test_render_dotplot_png_boundary_line_at_exact_offset_end() -> None:
    """Pink boundary line must be at the exact offset_end coordinate, not bin-snapped.

    With resolution=8 and an asymmetric 700/1000 split, the old _bin_boundary formula
    returned 750.0 instead of 700.0 — a 50-unit error visible as a misaligned pink line.
    The line must sit at float(block["offset_end"]) because the imshow extent maps
    offset 0→total_len directly to visual data coordinates.
    """
    from matplotlib.figure import Figure

    matrix = np.eye(8, dtype=np.float32)
    region = Region("chr8", 0, 1000, "chr8:0-1000")
    blocks = [
        {
            "label": "chr8:0.0-0.0 Mb",
            "chromosome": "chr8",
            "start": 0,
            "end": 700,
            "offset_start": 0,
            "offset_end": 700,
        },
        {
            "label": "chr14:0.0-0.0 Mb",
            "chromosome": "chr14",
            "start": 0,
            "end": 300,
            "offset_start": 700,
            "offset_end": 1000,
        },
    ]

    captured_vline_xs: list[float] = []
    captured_hline_ys: list[float] = []
    original_add_subplot = Figure.add_subplot

    def patched_add_subplot(self, *args, **kwargs):
        ax = original_add_subplot(self, *args, **kwargs)
        orig_axvline = ax.axvline
        orig_axhline = ax.axhline

        def capture_axvline(x, **kw):
            captured_vline_xs.append(float(x))
            return orig_axvline(x, **kw)

        def capture_axhline(y, **kw):
            captured_hline_ys.append(float(y))
            return orig_axhline(y, **kw)

        ax.axvline = capture_axvline
        ax.axhline = capture_axhline
        return ax

    with unittest.mock.patch.object(Figure, "add_subplot", patched_add_subplot):
        _render_dotplot_png_data_url(
            matrix,
            region,
            with_axes=False,
            figure_size=(1.0, 1.0),
            dpi=64,
            blocks=blocks,
        )

    # Left edge of the first chr20 bin: floor(700 * 8 / 1000) * 1000/8 = 5 * 125 = 625.0.
    # This aligns the pink line with where the N-masked grey area visually ends.
    assert captured_vline_xs == [pytest.approx(625.0)], (
        f"axvline called with {captured_vline_xs}, expected [625.0]"
    )
    assert captured_hline_ys == [pytest.approx(625.0)], (
        f"axhline called with {captured_hline_ys}, expected [625.0]"
    )


def test_create_dotplot_thumbnail_returns_single_clickable_thumbnail_button() -> None:
    matrix = np.eye(4, dtype=np.float32)
    region = Region("chr1", 100, 104, "chr1:100-104")

    thumbnail = create_dotplot_thumbnail(None, matrix, region, size=32)

    assert thumbnail.label == "Show ref identity"
    assert thumbnail.width == 156
    assert thumbnail.height == 30
    assert thumbnail.styles["cursor"] == "pointer"
    assert thumbnail.styles["min-width"] == "156px"
    assert "background-position: 9px center, 8px center, 7px center" in thumbnail.stylesheets[0]
    assert "background-size: 22px 22px, 24px 24px, 26px 26px" in thumbnail.stylesheets[0]
    assert "padding: 0 8px 0 42px" in thumbnail.stylesheets[0]
    assert "button_click" in thumbnail.js_event_callbacks


def test_create_combined_dotplot_thumbnail_uses_combined_modal_title() -> None:
    matrix = np.eye(4, dtype=np.float32)
    region = Region("chr1", 100, 104, "chr1:100-104")
    blocks = [
        {
            "label": "chr1:0.0-0.0 Mb",
            "chromosome": "chr1",
            "start": 100,
            "end": 104,
            "offset_start": 0,
            "offset_end": 4,
        }
    ]

    thumbnail = create_dotplot_thumbnail(
        None,
        matrix,
        region,
        size=32,
        blocks=blocks,
        region_label="Regions are concatenated: chr1:100-104",
        modal_title="Combined reference self-identity dotplot",
    )
    callback = thumbnail.js_event_callbacks["button_click"][0]

    assert callback.args["title_text"] == "Combined reference self-identity dotplot"
    assert callback.args["region_label"] == "Regions are concatenated: chr1:100-104"
    assert callback.args["individual_images"] == []


def test_create_dotplot_thumbnail_with_individual_payloads_passes_images_to_callback() -> None:
    matrix = np.eye(8, dtype=np.float32)
    region = Region("chr1", 100, 108, "chr1:100-108")
    blocks = [
        {
            "label": "chr1:0.0-0.0 Mb",
            "chromosome": "chr1",
            "start": 100,
            "end": 104,
            "offset_start": 0,
            "offset_end": 4,
        },
        {
            "label": "chr2:0.0-0.0 Mb",
            "chromosome": "chr2",
            "start": 200,
            "end": 204,
            "offset_start": 4,
            "offset_end": 8,
        },
    ]
    individual_payloads = [
        {
            "matrix": np.eye(4, dtype=np.float32),
            "region": Region("chr1", 100, 104, "chr1:100-104"),
            "label": "chr1:100-104",
            "title": "chr1 self-identity",
        },
        {
            "matrix": np.eye(4, dtype=np.float32),
            "region": Region("chr2", 200, 204, "chr2:200-204"),
            "label": "chr2:200-204",
            "title": "chr2 self-identity",
        },
    ]

    thumbnail = create_dotplot_thumbnail(
        None,
        matrix,
        region,
        size=32,
        blocks=blocks,
        individual_payloads=individual_payloads,
    )
    callback = thumbnail.js_event_callbacks["button_click"][0]
    indiv = callback.args["individual_images"]

    assert len(indiv) == 2
    assert indiv[0]["label"] == "chr1:100-104"
    assert indiv[0]["title"] == "chr1 self-identity"
    assert indiv[0]["url"].startswith("data:image/png;base64,")
    assert indiv[1]["label"] == "chr2:200-204"
    assert indiv[1]["url"].startswith("data:image/png;base64,")


def test_create_dotplot_thumbnail_without_individual_payloads_passes_empty_list() -> None:
    matrix = np.eye(4, dtype=np.float32)
    region = Region("chr1", 100, 104, "chr1:100-104")

    thumbnail = create_dotplot_thumbnail(None, matrix, region, size=32)
    callback = thumbnail.js_event_callbacks["button_click"][0]

    assert callback.args["individual_images"] == []


def test_create_coordinate_display_adds_phase_set_marker_checkbox() -> None:
    plot_figure = figure(x_range=(100, 200))
    phase_renderer = plot_figure.segment([120], [0], [120], [1])

    controls, _hide_1bp_checkbox, _coord_input = create_coordinate_display(
        plot_figure,
        "chr1",
        100,
        200,
        phase_set_marker_renderers=[phase_renderer],
    )

    checkbox_by_label = {checkbox.label: checkbox for checkbox in _collect_checkboxes(controls)}
    checkbox = checkbox_by_label["Hide phaseset markers"]
    callbacks = checkbox.js_property_callbacks["change:active"]

    assert callbacks[0].args["phase_set_marker_renderers"] == [phase_renderer]


def test_create_coordinate_display_enables_hold_for_pan_arrows() -> None:
    plot_figure = figure(x_range=(100, 200))

    controls, _hide_1bp_checkbox, _cursor_guide_checkbox = create_coordinate_display(
        plot_figure, "chr1", 100, 200
    )

    hold_buttons = [
        button.label
        for button in controls.select({"type": Button})
        if "orographer-click-hold" in button.tags
    ]
    assert hold_buttons == ["\u2039", "\u203a"]


def test_create_global_checkbox_controls_groups_evidence_and_display_controls() -> None:
    plot_figure = figure(x_range=(100, 200))
    phase_renderer = plot_figure.segment([120], [0], [120], [1])

    controls, hide_1bp_checkbox, cursor_guide_checkbox = create_global_checkbox_controls(
        plot_figure,
        one_bp_renderers=[phase_renderer],
        alignment_label_renderers=[phase_renderer],
        phase_set_marker_renderers=[phase_renderer],
        enable_read_filter_checkboxes=True,
    )

    assert controls is not None
    assert hide_1bp_checkbox is not None
    assert cursor_guide_checkbox is not None
    assert controls.styles["background"] == "#f8fafc"
    assert controls.styles["flex-wrap"] == "wrap"
    assert controls.styles["justify-content"] == "center"
    labels = [child.label for child in controls.children if hasattr(child, "label")]
    assert labels == [
        "Show only split algns",
        "Show only multiregion algns",
        "Show cursor guide",
        "Hide 1bp INDELs",
        "Hide algn numbers",
        "Hide phaseset markers",
    ]
    widths = [child.width for child in controls.children if hasattr(child, "label")]
    assert widths == [155, 200, 145, 125, 135, 160]


def test_create_global_checkbox_controls_accepts_read_search_controls_without_checkboxes() -> None:
    plot_figure = figure(x_range=(100, 200))
    read_search_controls = Div(text="Select reads")
    dotplot_thumbnail = Div(text="self-identity")

    controls, hide_1bp_checkbox, cursor_guide_checkbox = create_global_checkbox_controls(
        plot_figure,
        read_search_controls=read_search_controls,
        dotplot_thumbnail=dotplot_thumbnail,
    )

    assert controls is not None
    assert hide_1bp_checkbox is None
    assert cursor_guide_checkbox is not None
    trailing_group = controls.children[-2]
    assert read_search_controls in trailing_group.children
    assert dotplot_thumbnail in trailing_group.children
    assert trailing_group.spacing == 12
    assert trailing_group.margin == (6, 6, 0, 0)
    assert trailing_group.width == 12
    assert trailing_group.sizing_mode == "fixed"


def test_cursor_guide_adds_span_and_event_callbacks() -> None:
    plot_figure = figure(x_range=(100, 200))
    shared_figure = figure(x_range=plot_figure.x_range)
    spans = add_cursor_guide_to_figures([plot_figure, shared_figure])
    checkbox = Checkbox(label="Show cursor guide", active=True)

    add_cursor_guide_callbacks([plot_figure, shared_figure], spans, checkbox)

    assert len(spans) == 2
    assert all(isinstance(span, Span) for span in spans)
    assert all(not span.visible for span in spans)
    assert spans[0] in plot_figure.center
    assert spans[1] in shared_figure.center
    assert "mousemove" in plot_figure.js_event_callbacks
    assert "mouseleave" in plot_figure.js_event_callbacks
    assert "change:active" in checkbox.js_property_callbacks


def test_create_repeat_density_figure_returns_figure() -> None:
    main = figure(x_range=(1_000, 2_000))
    density = np.zeros(512, dtype=np.float32)
    density[100] = 50.0
    fig, _src = create_repeat_density_figure(main.x_range, density, 1_000, 2_000)
    assert hasattr(fig, "renderers") and hasattr(fig, "x_range")


def test_create_repeat_density_figure_shares_x_range() -> None:
    main = figure(x_range=(500, 1_500))
    density = np.zeros(512, dtype=np.float32)
    fig, _src = create_repeat_density_figure(main.x_range, density, 500, 1_500)
    assert fig.x_range is main.x_range


def test_create_repeat_density_figure_has_glyph_renderer() -> None:
    main = figure(x_range=(0, 10_000))
    density = np.ones(512, dtype=np.float32) * 3.0
    fig, _src = create_repeat_density_figure(main.x_range, density, 0, 10_000)
    assert len(fig.renderers) > 0


def test_create_repeat_density_figure_height() -> None:
    main = figure(x_range=(0, 5_000))
    density = np.zeros(512, dtype=np.float32)
    fig, _src = create_repeat_density_figure(main.x_range, density, 0, 5_000)
    assert fig.height < 100


def test_create_repeat_density_figure_has_ident_label() -> None:
    main = figure(x_range=(0, 5_000))
    density = np.zeros(512, dtype=np.float32)
    fig, _src = create_repeat_density_figure(main.x_range, density, 0, 5_000)
    assert fig.yaxis[0].axis_label == "Ident"
