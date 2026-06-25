from __future__ import annotations

import functools
import gzip
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np
import pytest

from orographer.plot_bokeh.plot_bokeh import plot_reads_bokeh
from orographer.utils import COMPLEX_SV_REGION_TYPE, ISOSEQ_REGION_TYPE, OutputConfig, Region

if TYPE_CHECKING:
    from playwright.sync_api import Browser, Locator, Page, Playwright

playwright_sync = pytest.importorskip("playwright.sync_api")
expect = playwright_sync.expect
PlaywrightError = playwright_sync.Error
sync_playwright = playwright_sync.sync_playwright

BROWSER_NAMES = ("chromium", "firefox", "webkit")
SERVER_HOST = "127.0.0.1"
PAGE_LOAD_TIMEOUT_MS = 30_000
EXPECTED_INITIAL_COORDINATES = "chr1:100-200"
ZOOMED_COORDINATES = "chr1:125-175"
OVERLAY_GEOMETRY_TOLERANCE_PX = 2.0


class QuietRequestHandler(SimpleHTTPRequestHandler):
    """Serve generated plot files without noisy per-request test output."""

    def log_message(self, format: str, *args: object) -> None:
        return


@contextmanager
def serve_directory(directory: Path) -> Iterator[str]:
    handler = functools.partial(QuietRequestHandler, directory=str(directory))
    server = ThreadingHTTPServer((SERVER_HOST, 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://{SERVER_HOST}:{server.server_port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def create_segment(
    *,
    pos: int = 100,
    end: int = 180,
    alignment_order: int = 1,
    fwd_read_start: int = 0,
    fwd_read_end: int = 80,
    haplotype_tag: int | None = 1,
) -> SimpleNamespace:
    return SimpleNamespace(
        pos=pos,
        end=end,
        is_fwd_strand=True,
        haplotype_tag=haplotype_tag,
        readname="readA",
        chrom="chr1",
        color_tag=None,
        alignment_order=alignment_order,
        fwd_read_start=fwd_read_start,
        fwd_read_end=fwd_read_end,
        mismatches=[],
        insertions=[],
        deletions=[],
    )


def create_plot_page(tmp_path: Path) -> Path:
    region = Region("chr1", 100, 200, "chr1:100-200")
    output_file = plot_reads_bokeh(
        [
            {
                "region": region,
                "gene_annotations": [],
                "bam_rows": [
                    {
                        "segments_by_read": {"readA": [create_segment()]},
                        "region_type": COMPLEX_SV_REGION_TYPE,
                        "vcf_variants": [],
                        "sample_label": "browser smoke",
                    }
                ],
            }
        ],
        OutputConfig(str(tmp_path), "browser_smoke"),
    )
    assert output_file is not None
    return Path(output_file)


def create_multi_region_plot_page(tmp_path: Path) -> Path:
    first_segment = create_segment(
        pos=120,
        end=180,
        alignment_order=1,
        fwd_read_start=0,
        fwd_read_end=60,
        haplotype_tag=2,
    )
    second_segment = create_segment(
        pos=320,
        end=380,
        alignment_order=2,
        fwd_read_start=61,
        fwd_read_end=121,
        haplotype_tag=None,
    )
    first_segment.chrom = "chr1"
    second_segment.chrom = "chr2"
    solo_segment = create_segment(
        pos=185,
        end=195,
        alignment_order=1,
        fwd_read_start=0,
        fwd_read_end=10,
        haplotype_tag=1,
    )
    solo_segment.chrom = "chr1"
    solo_segment.readname = "readSolo"
    first_region_data = {
        "region": Region("chr1", 100, 200, "chr1:100-200"),
        "gene_annotations": [],
        "bam_rows": [
            {
                "segments_by_read": {
                    "readA_HP2": [first_segment],
                    "readSolo_HP1": [solo_segment],
                },
                "region_type": COMPLEX_SV_REGION_TYPE,
                "vcf_variants": [],
                "sample_label": "first region",
            }
        ],
    }
    first_region_data["dotplot_payload"] = {
        "matrix": np.eye(8, dtype=np.float32),
        "blocks": [
            {
                "label": "chr1:0.0-0.0 Mb",
                "coordinate_label": "chr1:100-200",
                "chromosome": "chr1",
                "start": 100,
                "end": 200,
                "offset_start": 0,
                "offset_end": 100,
                "span": 100,
            },
            {
                "label": "chr2:0.0-0.0 Mb",
                "coordinate_label": "chr2:300-400",
                "chromosome": "chr2",
                "start": 300,
                "end": 400,
                "offset_start": 100,
                "offset_end": 200,
                "span": 100,
            },
        ],
        "total_span": 200,
        "label": (
            "Regions are concatenated in display order for reference comparison only: "
            "chr1:100-200; chr2:300-400"
        ),
        "title": "Combined reference self-identity dotplot",
    }
    output_file = plot_reads_bokeh(
        [
            first_region_data,
            {
                "region": Region("chr2", 300, 400, "chr2:300-400"),
                "gene_annotations": [],
                "bam_rows": [
                    {
                        "segments_by_read": {"readA": [second_segment]},
                        "region_type": COMPLEX_SV_REGION_TYPE,
                        "vcf_variants": [],
                        "sample_label": "second region",
                    }
                ],
            },
        ],
        OutputConfig(str(tmp_path), "browser_multi_region"),
    )
    assert output_file is not None
    return Path(output_file)


def create_insertion_marker_plot_page(tmp_path: Path) -> Path:
    region = Region("chr1", 100, 5100, "chr1:100-5100")
    insertion_sites = [
        {
            "pos": 1000,
            "count": 1,
            "median_size": 50,
            "read_names": ["readA"],
            "total_count": 1,
            "chrom": "chr1",
        },
        {
            "pos": 1080,
            "count": 1,
            "median_size": 80,
            "read_names": ["readB"],
            "total_count": 1,
            "chrom": "chr1",
        },
    ]
    output_file = plot_reads_bokeh(
        [
            {
                "region": region,
                "gene_annotations": [],
                "bam_rows": [
                    {
                        "segments_by_read": {"readA": [create_segment(pos=900, end=1200)]},
                        "region_type": COMPLEX_SV_REGION_TYPE,
                        "vcf_variants": [],
                        "sample_label": "insertions",
                        "insertion_summary": {1: insertion_sites},
                    }
                ],
            }
        ],
        OutputConfig(str(tmp_path), "browser_insertions"),
    )
    assert output_file is not None
    return Path(output_file)


def create_isoseq_plot_page(tmp_path: Path) -> Path:
    transcript_one = SimpleNamespace(
        transcript_id="TX1",
        transcript_name="Isoform 1",
        gene_id="GENE1",
        gene_name="Gene One",
        chrom="chr1",
        start=100,
        end=260,
        strand="+",
        exons=[(100, 150, 1), (210, 260, 2)],
    )
    transcript_two = SimpleNamespace(
        transcript_id="TX2",
        transcript_name="Isoform 2",
        gene_id="GENE1",
        gene_name="Gene One",
        chrom="chr1",
        start=110,
        end=250,
        strand="+",
        exons=[(110, 250, 1)],
    )
    iso_segment = create_segment(
        pos=100,
        end=260,
        alignment_order=1,
        fwd_read_start=0,
        fwd_read_end=100,
        haplotype_tag=1,
    )
    iso_segment.readname = "isoRead1"
    iso_segment.aligned_blocks = [(100, 150), (210, 260)]
    unassigned_segment = create_segment(
        pos=120,
        end=240,
        alignment_order=1,
        fwd_read_start=0,
        fwd_read_end=80,
        haplotype_tag=0,
    )
    unassigned_segment.readname = "unassignedRead1"
    unassigned_segment.aligned_blocks = [(120, 240)]
    region_data = {
        "region": Region("chr1", 90, 270, "chr1:90-270"),
        "gene_annotations": [],
        "bam_rows": [
            {
                "segments_by_read": {
                    "isoRead1": [iso_segment],
                    "unassignedRead1": [unassigned_segment],
                },
                "isoseq_groups": [
                    {
                        "group_id": "TX1",
                        "transcript": transcript_one,
                        "read_names": ["isoRead1"],
                        "assigned_read_count": 1,
                        "gene_transcript_count": 2,
                    },
                    {
                        "group_id": "TX2",
                        "transcript": transcript_two,
                        "read_names": [],
                        "assigned_read_count": 0,
                        "gene_transcript_count": 2,
                    },
                    {
                        "group_id": "UNASSIGNED",
                        "transcript": None,
                        "read_names": ["unassignedRead1"],
                        "assigned_read_count": 1,
                        "gene_transcript_count": 0,
                    },
                ],
                "region_type": ISOSEQ_REGION_TYPE,
                "vcf_variants": [],
                "sample_label": "isoseq browser smoke",
            }
        ],
    }
    region_data["dotplot_payload"] = {
        "matrix": np.eye(8, dtype=np.float32),
        "blocks": [],
        "label": "chr1:90-270",
        "title": "Reference self-identity dotplot",
    }
    output_file = plot_reads_bokeh(
        [region_data],
        OutputConfig(str(tmp_path), "browser_isoseq"),
    )
    assert output_file is not None
    return Path(output_file)


def launch_browser(playwright: Playwright, browser_name: str) -> Browser:
    browser_type = getattr(playwright, browser_name)
    try:
        return browser_type.launch()
    except PlaywrightError as exc:
        message = str(exc)
        missing_binary = "Executable doesn't exist" in message
        missing_host_dependency = "Host system is missing dependencies" in message
        if missing_binary:
            pytest.skip(
                "Playwright browser binaries are not installed. Run "
                "`pixi run install-browsers` before browser tests."
            )
        if missing_host_dependency:
            pytest.skip(
                "Playwright browser host dependencies are missing. Run "
                "`playwright install-deps` or install the packages listed by Playwright."
            )
        raise


def coordinate_input(page: Page) -> Locator:
    return page.locator("input[type='text']").first


def attach_error_collectors(page: Page) -> tuple[list[str], list[str]]:
    console_errors: list[str] = []
    page_errors: list[str] = []
    page.on(
        "console",
        lambda message: console_errors.append(message.text) if message.type == "error" else None,
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    return console_errors, page_errors


def assert_no_browser_errors(console_errors: list[str], page_errors: list[str]) -> None:
    assert page_errors == []
    assert console_errors == []


@contextmanager
def open_generated_plot_page(
    tmp_path: Path,
    browser_name: str,
) -> Iterator[tuple[Page, list[str], list[str]]]:
    plot_page = create_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            yield page, console_errors, page_errors
        finally:
            browser.close()


@pytest.fixture(params=BROWSER_NAMES)
def generated_plot_page(
    tmp_path: Path,
    request: pytest.FixtureRequest,
) -> Iterator[tuple[Page, list[str], list[str]]]:
    with open_generated_plot_page(tmp_path, str(request.param)) as page_data:
        yield page_data


@pytest.fixture(params=BROWSER_NAMES)
def browser_name(request: pytest.FixtureRequest) -> str:
    return str(request.param)


def connector_selection_state(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const states = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                if (model.data && model.data.source_kind) {
                    if (model.data.source_kind[0] === "connector") {
                        states.push({
                            count: model.data.read_name.length,
                            selected: model.selected.indices.slice(),
                        });
                    }
                }
            }
            return { count: states.length, states };
        }
        """
    )


def insertion_marker_state(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const states = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (data && data.cluster_pos && data.marker_label) {
                    states.push({
                        clusterPos: data.cluster_pos.map((item) => item.slice()),
                        markerLabel: data.marker_label.slice(),
                        sourceSize: data.marker_label.length,
                    });
                }
            }
            return states[0] || { clusterPos: [], markerLabel: [], sourceSize: 0 };
        }
        """
    )


def overlay_line_state(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const paths = Array.from(
                document.querySelectorAll(
                    "#orographerReadConnectionOverlay "
                    + ".orographer-read-connection-overlay-line"
                )
            );
            return {
                count: paths.length,
                lines: paths.map((path) => ({
                    x1: Number(path.getAttribute("data-x1")),
                    y1: Number(path.getAttribute("data-y1")),
                    x2: Number(path.getAttribute("data-x2")),
                    y2: Number(path.getAttribute("data-y2")),
                    opacity: path.getAttribute("stroke-opacity"),
                    stroke: path.getAttribute("stroke"),
                    d: path.getAttribute("d"),
                })),
            };
        }
        """
    )


def overlay_line_midpoint(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const path = document.querySelector(
                "#orographerReadConnectionOverlay "
                + ".orographer-read-connection-overlay-line"
            );
            const point = path.getPointAtLength(path.getTotalLength() / 2);
            return { x: point.x, y: point.y };
        }
        """
    )


def overlay_modal_z_indices(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const overlay = document.getElementById("orographerReadConnectionOverlay");
            const modal = document.getElementById("alignmentModal");
            return {
                overlay: Number(window.getComputedStyle(overlay).zIndex),
                modal: Number(window.getComputedStyle(modal).zIndex),
            };
        }
        """
    )


def connection_debug_state(page: Page) -> dict:
    return page.evaluate("() => window.orographerReadConnectionDebugState()")


def combined_dotplot_state(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const callbacks = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const events = model.js_event_callbacks || {};
                const clickCallbacks = (events.tap || []).concat(events.button_click || []);
                clickCallbacks.forEach((callback) => {
                    if (callback.args && callback.args.title_text) {
                        callbacks.push({
                            modelId: model.id,
                            title: callback.args.title_text,
                            label: callback.args.region_label,
                        });
                    }
                });
            }
            return { count: callbacks.length, callbacks };
        }
        """
    )


def read_filter_alpha_state(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const rows = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (data && data.source_kind && data.source_kind[0] === "arrow") {
                    data.read_name.forEach((readName, rowIndex) => {
                        rows.push({
                            readName,
                            layoutReadName: data.layout_read_name
                                ? data.layout_read_name[rowIndex]
                                : readName,
                            alpha: data.read_filter_alpha[rowIndex],
                            y: data.y[rowIndex],
                            yAll: data.y_all ? data.y_all[rowIndex] : null,
                            ySplit: data.y_split ? data.y_split[rowIndex] : null,
                            yMultiRegion: data.y_multiregion
                                ? data.y_multiregion[rowIndex]
                                : null,
                            ySplitMultiRegion: data.y_split_multiregion
                                ? data.y_split_multiregion[rowIndex]
                                : null,
                            hasSplit: data.has_split_alignment[rowIndex],
                            hasMultiRegion: data.has_multiregion_connection[rowIndex],
                        });
                    });
                }
            }
            return { rows };
        }
        """
    )


def set_bokeh_checkbox(page: Page, label: str, active: bool) -> None:
    page.evaluate(
        """
        ({ label, active }) => {
            for (const model of Bokeh.documents[0]._all_models.values()) {
                if (model.label === label) {
                    model.setv({ active });
                    return;
                }
            }
            throw new Error("Missing checkbox: " + label);
        }
        """,
        {"label": label, "active": active},
    )


def bokeh_checkbox_label_counts(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const counts = {};
            for (const model of Bokeh.documents[0]._all_models.values()) {
                if (typeof model.label === "string" && model.active !== undefined) {
                    counts[model.label] = counts[model.label] ? counts[model.label] + 1 : 1;
                }
            }
            return counts;
        }
        """
    )


def alignment_plot_tops(page: Page) -> list[float]:
    return page.evaluate(
        """
        () => {
            const tops = [];
            const seen = new Set();
            function visit(value) {
                if (!value || typeof value !== "object" || seen.has(value)) return;
                seen.add(value);
                const cssClasses = value.model ? value.model.css_classes || [] : [];
                const isAlignmentPlot = cssClasses.some((item) => {
                    return item.indexOf("orographer-alignment-plot-r") === 0;
                });
                if (isAlignmentPlot && value.canvas_view) {
                    const rect = value.canvas_view.events_el.getBoundingClientRect();
                    tops.push(rect.top);
                }
                if (value[Symbol.iterator]) {
                    for (const item of value) visit(item);
                }
                for (const key of Object.keys(value)) {
                    if (key === "parent" || key === "root" || key === "owner") continue;
                    visit(value[key]);
                }
            }
            for (const rootView of Object.values(Bokeh.index)) {
                visit(rootView);
            }
            return tops.sort((a, b) => a - b);
        }
        """
    )


def combined_dotplot_screen_position(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            let targetId = null;
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const events = model.js_event_callbacks || {};
                const clickCallbacks = (events.tap || []).concat(events.button_click || []);
                clickCallbacks.forEach((callback) => {
                    if (callback.args && callback.args.title_text) {
                        targetId = model.id;
                    }
                });
            }
            if (!targetId) return null;
            const seen = new Set();
            function visit(value) {
                if (!value || typeof value !== "object" || seen.has(value)) return null;
                seen.add(value);
                if (value.model && value.model.id === targetId && value.canvas_view) {
                    const rect = value.canvas_view.events_el.getBoundingClientRect();
                    return {
                        x: rect.left + rect.width / 2,
                        y: rect.top + rect.height / 2,
                    };
                }
                if (value.model && value.model.id === targetId && value.el) {
                    const rect = value.el.getBoundingClientRect();
                    return {
                        x: rect.left + rect.width / 2,
                        y: rect.top + rect.height / 2,
                    };
                }
                if (value[Symbol.iterator]) {
                    for (const item of value) {
                        const found = visit(item);
                        if (found) return found;
                    }
                }
                for (const key of Object.keys(value)) {
                    if (key === "parent" || key === "root" || key === "owner") continue;
                    const found = visit(value[key]);
                    if (found) return found;
                }
                return null;
            }
            for (const rootView of Object.values(Bokeh.index)) {
                const found = visit(rootView);
                if (found) return found;
            }
            return null;
        }
        """
    )


def combined_dotplot_screen_rect(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            let targetId = null;
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const events = model.js_event_callbacks || {};
                const clickCallbacks = (events.tap || []).concat(events.button_click || []);
                clickCallbacks.forEach((callback) => {
                    if (callback.args && callback.args.title_text) {
                        targetId = model.id;
                    }
                });
            }
            if (!targetId) return null;
            const seen = new Set();
            function visit(value) {
                if (!value || typeof value !== "object" || seen.has(value)) return null;
                seen.add(value);
                if (value.model && value.model.id === targetId && value.canvas_view) {
                    const rect = value.canvas_view.events_el.getBoundingClientRect();
                    return {
                        left: rect.left,
                        right: rect.right,
                        top: rect.top,
                        width: rect.width,
                        height: rect.height,
                    };
                }
                if (value.model && value.model.id === targetId && value.el) {
                    const rect = value.el.getBoundingClientRect();
                    return {
                        left: rect.left,
                        right: rect.right,
                        top: rect.top,
                        width: rect.width,
                        height: rect.height,
                    };
                }
                if (value[Symbol.iterator]) {
                    for (const item of value) {
                        const found = visit(item);
                        if (found) return found;
                    }
                }
                for (const key of Object.keys(value)) {
                    if (key === "parent" || key === "root" || key === "owner") continue;
                    const found = visit(value[key]);
                    if (found) return found;
                }
                return null;
            }
            for (const rootView of Object.values(Bokeh.index)) {
                const found = visit(rootView);
                if (found) return found;
            }
            return null;
        }
        """
    )


def alignment_label_screen_position(page: Page, read_name: str) -> dict:
    return page.evaluate(
        """
        ({ readName }) => {
            let dataX = null;
            let dataY = null;
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (!data || !data.alignment_number || !data.read_name) continue;
                const rowIndex = data.read_name.indexOf(readName);
                if (rowIndex === -1) continue;
                dataX = data.x[rowIndex];
                dataY = data.y[rowIndex];
                break;
            }
            if (dataX === null || dataY === null) return null;

            const seen = new Set();
            function visit(value) {
                if (!value || typeof value !== "object" || seen.has(value)) return null;
                seen.add(value);
                const model = value.model ? value.model : null;
                const xRange = model ? model.x_range : null;
                const yRange = model ? model.y_range : null;
                if (value.frame && xRange && yRange) {
                    const insideX = dataX >= xRange.start && dataX <= xRange.end;
                    const insideForwardY = dataY >= yRange.start && dataY <= yRange.end;
                    const insideReverseY = dataY >= yRange.end && dataY <= yRange.start;
                    if (insideX && (insideForwardY || insideReverseY)) {
                        const rect = value.canvas_view.events_el.getBoundingClientRect();
                        return {
                            x: rect.left + value.frame.x_scale.compute(dataX),
                            y: rect.top + value.frame.y_scale.compute(dataY),
                        };
                    }
                }
                if (value[Symbol.iterator]) {
                    for (const item of value) {
                        const found = visit(item);
                        if (found) return found;
                    }
                }
                for (const key of Object.keys(value)) {
                    if (key === "parent" || key === "root" || key === "owner") continue;
                    const found = visit(value[key]);
                    if (found) return found;
                }
                return null;
            }
            for (const rootView of Object.values(Bokeh.index)) {
                const found = visit(rootView);
                if (found) return found;
            }
            return null;
        }
        """,
        {"readName": read_name},
    )


def isoseq_row_screen_position(page: Page, transcript_id: str, source_kind: str) -> dict:
    return page.evaluate(
        """
        ({ transcriptId, sourceKind }) => {
            let dataX = null;
            let dataY = null;
            let targetSourceId = null;
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (!data || !data.transcript_id || !data.source_kind) continue;
                const rowIndex = data.transcript_id.findIndex((value, index) => {
                    return value === transcriptId && data.source_kind[index] === sourceKind;
                });
                if (rowIndex === -1) continue;
                dataX = (data.left[rowIndex] + data.right[rowIndex]) / 2;
                dataY = (data.top[rowIndex] + data.bottom[rowIndex]) / 2;
                targetSourceId = model.id;
                break;
            }
            if (dataX === null || dataY === null) return null;

            const seen = new Set();
            function rendersTargetSource(model) {
                const renderers = model && model.renderers ? model.renderers : [];
                return renderers.some((renderer) => {
                    const source = renderer.data_source;
                    return source && source.id === targetSourceId;
                });
            }

            function visit(value) {
                if (!value || typeof value !== "object" || seen.has(value)) return null;
                seen.add(value);
                const model = value.model ? value.model : null;
                if (model && model.renderers) {
                    if (!rendersTargetSource(model)) return null;
                }
                const xRange = model ? model.x_range : null;
                const yRange = model ? model.y_range : null;
                if (value.frame && xRange && yRange) {
                    const insideX = dataX >= xRange.start && dataX <= xRange.end;
                    const insideForwardY = dataY >= yRange.start && dataY <= yRange.end;
                    const insideReverseY = dataY >= yRange.end && dataY <= yRange.start;
                    if (insideX && (insideForwardY || insideReverseY)) {
                        const rect = value.canvas_view.events_el.getBoundingClientRect();
                        return {
                            x: rect.left + value.frame.x_scale.compute(dataX),
                            y: rect.top + value.frame.y_scale.compute(dataY),
                        };
                    }
                }
                if (value[Symbol.iterator]) {
                    for (const item of value) {
                        const found = visit(item);
                        if (found) return found;
                    }
                }
                for (const key of Object.keys(value)) {
                    if (key === "parent" || key === "root" || key === "owner") continue;
                    const found = visit(value[key]);
                    if (found) return found;
                }
                return null;
            }
            for (const rootView of Object.values(Bokeh.index)) {
                const found = visit(rootView);
                if (found) return found;
            }
            return null;
        }
        """,
        {"transcriptId": transcript_id, "sourceKind": source_kind},
    )


def isoseq_transcript_screen_position(page: Page, transcript_id: str) -> dict:
    return isoseq_row_screen_position(page, transcript_id, "transcript")


def isoseq_unassigned_screen_position(page: Page) -> dict:
    return isoseq_row_screen_position(page, "UNASSIGNED", "unassigned")


def isoseq_transcript_alpha_state(page: Page) -> dict:
    return page.evaluate(
        """
        () => {
            const result = {};
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (!data || !data.transcript_id || !data.source_kind || !data.alpha) continue;
                data.transcript_id.forEach((transcriptId, index) => {
                    if (data.source_kind[index] !== "transcript") return;
                    result[transcriptId] = data.alpha[index];
                });
            }
            return result;
        }
        """
    )


def select_isoseq_transcript_option(page: Page, transcript_id: str) -> None:
    page.evaluate(
        """
        ({ transcriptId }) => {
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const hasOptions = Array.isArray(model.options);
                const hasTranscript = hasOptions
                    ? model.options.some((option) => option[0] === transcriptId)
                    : false;
                if (hasTranscript) {
                    model.setv({ value: transcriptId });
                    return;
                }
            }
            throw new Error("Missing IsoSeq transcript selector");
        }
        """,
        {"transcriptId": transcript_id},
    )


def isoseq_loaded_read_names(page: Page) -> list[str]:
    return page.evaluate(
        """
        () => {
            const names = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (!data || !data.source_kind || !data.read_name) continue;
                if (data.source_kind[0] !== "arrow") continue;
                for (const readName of data.read_name) {
                    if (readName && names.indexOf(readName) === -1) names.push(readName);
                }
            }
            return names.sort();
        }
        """
    )


def isoseq_read_source_column_mismatches(page: Page) -> list[dict]:
    return page.evaluate(
        """
        () => {
            const mismatches = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (!data || !data.read_name) continue;
                if (!data.x0 && !data.layout_read_name) continue;
                const lengths = {};
                Object.keys(data).forEach((key) => {
                    lengths[key] = Array.isArray(data[key]) ? data[key].length : -1;
                });
                const expected = Math.max(...Object.values(lengths));
                const badColumns = Object.keys(lengths).filter((key) => {
                    return lengths[key] !== expected;
                });
                if (badColumns.length) {
                    mismatches.push({ sourceId: model.id, expected, badColumns, lengths });
                }
            }
            return mismatches;
        }
        """
    )


def assert_close(actual: float, expected: float, tolerance: float) -> None:
    assert abs(actual - expected) <= tolerance


def assert_overlay_line_matches_connection_endpoints(debug_state: dict) -> None:
    assert len(debug_state["lines"]) == 1
    line = debug_state["lines"][0]
    endpoints = debug_state["endpoints"][line["connectionId"]]
    assert len(endpoints) == 2
    first, second = endpoints

    assert_close(line["x1"], first["x"], OVERLAY_GEOMETRY_TOLERANCE_PX)
    assert_close(line["y1"], first["y"], OVERLAY_GEOMETRY_TOLERANCE_PX)
    assert_close(line["x2"], second["x"], OVERLAY_GEOMETRY_TOLERANCE_PX)
    assert_close(line["y2"], second["y"], OVERLAY_GEOMETRY_TOLERANCE_PX)


def assert_endpoint_matches_independent_range_projection(endpoint: dict) -> None:
    frame = endpoint["frame"]
    x_range = endpoint["xRange"]
    y_range = endpoint["yRange"]
    x_span = x_range["end"] - x_range["start"]
    y_span = y_range["start"] - y_range["end"]
    expected_x = frame["left"] + frame["width"] * ((endpoint["dataX"] - x_range["start"]) / x_span)
    expected_y = frame["top"] + frame["height"] * ((endpoint["dataY"] - y_range["end"]) / y_span)

    assert_close(endpoint["x"], expected_x, OVERLAY_GEOMETRY_TOLERANCE_PX)
    assert_close(endpoint["y"], expected_y, OVERLAY_GEOMETRY_TOLERANCE_PX)


def assert_endpoint_matches_bokeh_scale(page: Page, endpoint: dict) -> None:
    expected = page.evaluate(
        """
        ({ plotModelId, dataX, dataY }) => {
            const seen = new Set();
            function visit(value) {
                if (!value || typeof value !== "object" || seen.has(value)) return null;
                seen.add(value);
                if (value.model && value.model.id === plotModelId && value.frame) {
                    const rect = value.canvas_view.events_el.getBoundingClientRect();
                    return {
                        x: rect.left + value.frame.x_scale.compute(dataX),
                        y: rect.top + value.frame.y_scale.compute(dataY),
                    };
                }
                if (value[Symbol.iterator]) {
                    for (const item of value) {
                        const found = visit(item);
                        if (found) return found;
                    }
                }
                for (const key of Object.keys(value)) {
                    if (key === "parent" || key === "root" || key === "owner") continue;
                    const found = visit(value[key]);
                    if (found) return found;
                }
                return null;
            }
            for (const rootView of Object.values(Bokeh.index)) {
                const found = visit(rootView);
                if (found) return found;
            }
            return null;
        }
        """,
        {
            "plotModelId": endpoint["plotModelId"],
            "dataX": endpoint["dataX"],
            "dataY": endpoint["dataY"],
        },
    )
    assert expected is not None
    assert_close(endpoint["x"], expected["x"], OVERLAY_GEOMETRY_TOLERANCE_PX)
    assert_close(endpoint["y"], expected["y"], OVERLAY_GEOMETRY_TOLERANCE_PX)


def assert_debug_state_matches_bokeh_scale(page: Page, debug_state: dict) -> None:
    for endpoints in debug_state["endpoints"].values():
        for endpoint in endpoints:
            assert_endpoint_matches_bokeh_scale(page, endpoint)


@pytest.mark.browser
def test_generated_plot_page_loads_without_console_errors(
    generated_plot_page: tuple[Page, list[str], list[str]],
) -> None:
    page, console_errors, page_errors = generated_plot_page

    expect(coordinate_input(page)).to_have_value(EXPECTED_INITIAL_COORDINATES)
    expect(page.get_by_role("button", name="Go")).to_be_visible()
    expect(page.get_by_role("button", name="Select reads")).to_be_visible()
    expect(page.get_by_role("button", name="Select reads")).to_have_count(1)
    expect(page.get_by_role("button", name="Clear selected")).to_have_count(1)
    expect(page.locator("canvas").first).to_be_visible()

    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_multi_region_plot_contains_connector_sources(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            page.wait_for_function(
                "() => document.querySelectorAll("
                "'#orographerReadConnectionOverlay .orographer-read-connection-overlay-line'"
                ").length === 1",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            state = connector_selection_state(page)
            overlay = overlay_line_state(page)
            debug_state = connection_debug_state(page)
            assert_debug_state_matches_bokeh_scale(page, debug_state)
        finally:
            browser.close()

    assert state["count"] == 2
    assert [item["count"] for item in state["states"]] == [1, 1]
    assert overlay["count"] == 1
    assert overlay["lines"][0]["x1"] < overlay["lines"][0]["x2"]
    assert overlay["lines"][0]["opacity"] == "0.32"
    assert overlay["lines"][0]["stroke"] == "lightgrey"
    assert overlay["lines"][0]["d"].startswith("M ")
    assert_overlay_line_matches_connection_endpoints(debug_state)
    for endpoints in debug_state["endpoints"].values():
        for endpoint in endpoints:
            assert_endpoint_matches_independent_range_projection(endpoint)
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_multi_region_plot_shows_one_combined_dotplot_thumbnail(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            state = combined_dotplot_state(page)
            position = combined_dotplot_screen_position(page)
            rect = combined_dotplot_screen_rect(page)
            select_box = page.get_by_role("button", name="Select reads").bounding_box()
            clear_box = page.get_by_role("button", name="Clear selected").bounding_box()
            assert position is not None
            assert rect is not None
            assert select_box is not None
            assert clear_box is not None
            expect(page.get_by_text("Show ref identity")).to_be_visible()
            page.mouse.click(rect["left"] + 20, rect["top"] + rect["height"] / 2)
            expect(page.locator("#alignmentModal")).to_be_visible()
            expect(page.locator("#alignmentModalTitle")).to_contain_text(
                "Combined reference self-identity dotplot"
            )
            expect(page.locator("#modalContent")).to_contain_text("Regions are concatenated")
        finally:
            browser.close()

    assert state["count"] == 1
    assert state["callbacks"][0]["title"] == "Combined reference self-identity dotplot"
    assert "chr1:100-200" in state["callbacks"][0]["label"]
    assert "chr2:300-400" in state["callbacks"][0]["label"]
    assert rect["width"] == 156
    assert rect["height"] == 30
    assert select_box["x"] + select_box["width"] <= clear_box["x"]
    assert clear_box["x"] + clear_box["width"] <= rect["left"]
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_read_filter_checkboxes_hide_reads_by_evidence_class(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )

            before = read_filter_alpha_state(page)
            checkbox_counts = bokeh_checkbox_label_counts(page)
            plot_tops = alignment_plot_tops(page)
            set_bokeh_checkbox(page, "Show only split algns", True)
            after_split = read_filter_alpha_state(page)
            set_bokeh_checkbox(page, "Show only multiregion algns", True)
            after_multiregion = read_filter_alpha_state(page)
        finally:
            browser.close()

    before_by_read = {row["readName"]: row for row in before["rows"]}
    split_by_read = {row["readName"]: row for row in after_split["rows"]}
    multiregion_by_read = {row["readName"]: row for row in after_multiregion["rows"]}
    before_by_layout = {row["layoutReadName"]: row for row in before["rows"]}
    split_by_layout = {row["layoutReadName"]: row for row in after_split["rows"]}
    multiregion_by_layout = {row["layoutReadName"]: row for row in after_multiregion["rows"]}

    assert before_by_read["readA"]["alpha"] == 1
    assert before_by_read["readSolo"]["alpha"] == 1
    assert before_by_read["readA"]["hasSplit"] is True
    assert before_by_read["readA"]["hasMultiRegion"] is True
    assert before_by_read["readSolo"]["hasSplit"] is False
    assert before_by_read["readSolo"]["hasMultiRegion"] is False
    assert split_by_read["readA"]["alpha"] == 1
    assert split_by_read["readSolo"]["alpha"] == 0
    assert multiregion_by_read["readA"]["alpha"] == 1
    assert multiregion_by_read["readSolo"]["alpha"] == 0
    assert before_by_layout["readA_HP2"]["y"] == before_by_layout["readA_HP2"]["yAll"]
    assert split_by_layout["readA_HP2"]["y"] == split_by_layout["readA_HP2"]["ySplit"]
    assert split_by_layout["readA_HP2"]["y"] != before_by_layout["readA_HP2"]["y"]
    assert (
        multiregion_by_layout["readA_HP2"]["y"]
        == multiregion_by_layout["readA_HP2"]["ySplitMultiRegion"]
    )
    assert checkbox_counts.get("Hide 1bp INDELs", 0) <= 1
    assert checkbox_counts.get("Hide algn numbers", 0) <= 1
    assert checkbox_counts.get("Hide phaseset markers", 0) <= 1
    assert checkbox_counts["Show cursor guide"] == 1
    assert checkbox_counts["Show only split algns"] == 1
    assert checkbox_counts["Show only multiregion algns"] == 1
    assert len(plot_tops) == 2
    assert abs(plot_tops[0] - plot_tops[1]) <= 1
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_isoseq_unassigned_selection_keeps_isoforms_visible_by_default(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_isoseq_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )

            click_position = isoseq_unassigned_screen_position(page)
            assert click_position is not None
            page.mouse.click(click_position["x"], click_position["y"])
            page.wait_for_function(
                """
                () => {
                    for (const model of Bokeh.documents[0]._all_models.values()) {
                        const selected = model.selected ? model.selected.indices || [] : [];
                        const data = model.data;
                        if (!data || !data.source_kind || !data.transcript_id) continue;
                        if (!selected.length) continue;
                        const index = selected[selected.length - 1];
                        if (data.source_kind[index] !== "unassigned") continue;
                        return data.transcript_id[index] === "UNASSIGNED";
                    }
                    return false;
                }
                """,
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            alpha_state = isoseq_transcript_alpha_state(page)
        finally:
            browser.close()

    assert alpha_state["TX1"] > 0
    assert alpha_state["TX2"] > 0
    assert_no_browser_errors(console_errors, page_errors)


def test_isoseq_rerender_uses_new_chunk_generation_directory(tmp_path: Path) -> None:
    plot_page = create_isoseq_plot_page(tmp_path)
    chunk_root = tmp_path / f"{plot_page.stem}_chunks"
    first_tokens = {path.name for path in chunk_root.iterdir() if path.is_dir()}
    assert len(first_tokens) == 1
    first_token = next(iter(first_tokens))
    first_coverage_files = list((chunk_root / first_token).glob("g*_UNASSIGNED_coverage.json.gz"))
    assert len(first_coverage_files) == 1
    first_coverage_file = first_coverage_files[0]

    plot_page = create_isoseq_plot_page(tmp_path)
    second_tokens = {path.name for path in chunk_root.iterdir() if path.is_dir()}
    new_tokens = second_tokens - first_tokens
    assert len(new_tokens) == 1
    second_token = next(iter(new_tokens))

    assert first_coverage_file.exists()
    with gzip.open(plot_page.with_suffix(".json.gz"), "rt", encoding="utf-8") as handle:
        json_text = handle.read()
    assert f"{chunk_root.name}/{second_token}/" in json_text
    assert f"{chunk_root.name}/{first_token}/" not in json_text


@pytest.mark.browser
def test_isoseq_transcript_selection_loads_reads_without_browser_errors(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_isoseq_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )

            gene_box = page.locator("select").nth(1).bounding_box()
            dotplot_box = page.get_by_role("button", name="Show ref identity").bounding_box()
            isoform_box = page.get_by_text("Hide unselected isoforms").bounding_box()
            assert gene_box is not None
            assert dotplot_box is not None
            assert isoform_box is not None
            assert gene_box["x"] + gene_box["width"] <= dotplot_box["x"]
            assert dotplot_box["x"] + dotplot_box["width"] <= isoform_box["x"]

            select_isoseq_transcript_option(page, "TX1")
            page.wait_for_function(
                """
                () => {
                    for (const model of Bokeh.documents[0]._all_models.values()) {
                        const data = model.data;
                        if (!data || !data.read_name) continue;
                        if (data.read_name.indexOf("isoRead1") !== -1) return true;
                    }
                    return false;
                }
                """,
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            names_after_dropdown = isoseq_loaded_read_names(page)
            mismatches_after_dropdown = isoseq_read_source_column_mismatches(page)

            select_isoseq_transcript_option(page, "ALL")
            page.wait_for_function(
                "() => !Bokeh.documents[0]._all_models.values().next().done",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            click_position = isoseq_transcript_screen_position(page, "TX1")
            assert click_position is not None
            page.mouse.click(click_position["x"], click_position["y"])
            page.wait_for_function(
                """
                () => {
                    for (const model of Bokeh.documents[0]._all_models.values()) {
                        const data = model.data;
                        if (!data || !data.read_name) continue;
                        if (data.read_name.indexOf("isoRead1") !== -1) return true;
                    }
                    return false;
                }
                """,
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            names_after_click = isoseq_loaded_read_names(page)
            mismatches_after_click = isoseq_read_source_column_mismatches(page)
        finally:
            browser.close()

    assert names_after_dropdown == ["isoRead1"]
    assert names_after_click == ["isoRead1"]
    assert mismatches_after_dropdown == []
    assert mismatches_after_click == []
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_read_search_selects_multi_region_connectors(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            page.get_by_role("button", name="Select reads").first.evaluate(
                "button => button.click()"
            )
            z_indices = overlay_modal_z_indices(page)
            page.locator("#modalContent textarea").fill("readA")
            page.locator("#modalContent button").filter(has_text="Highlight").click()
            state = connector_selection_state(page)
            overlay_path = page.locator(
                "#orographerReadConnectionOverlay .orographer-read-connection-overlay-line"
            )
            expect(overlay_path).to_have_count(1)
            overlay = overlay_line_state(page)
            debug_state = connection_debug_state(page)
            assert_debug_state_matches_bokeh_scale(page, debug_state)
        finally:
            browser.close()

    assert [item["selected"] for item in state["states"]] == [[0], [0]]
    assert z_indices["overlay"] < z_indices["modal"]
    assert overlay["lines"][0]["opacity"] == "0.84"
    assert overlay["lines"][0]["stroke"] == "#92400e"
    assert_overlay_line_matches_connection_endpoints(debug_state)
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_clicking_multi_region_connector_selects_read_without_popup_errors(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            overlay_path = page.locator(
                "#orographerReadConnectionOverlay .orographer-read-connection-overlay-hit-area"
            )
            expect(overlay_path).to_have_count(1, timeout=PAGE_LOAD_TIMEOUT_MS)
            click_point = overlay_line_midpoint(page)

            page.mouse.click(click_point["x"], click_point["y"])
            state = connector_selection_state(page)
            modal_display = page.locator("#alignmentModal").evaluate("modal => modal.style.display")
        finally:
            browser.close()

    assert [item["selected"] for item in state["states"]] == [[0], [0]]
    assert modal_display != "flex"
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_multi_region_connector_line_follows_zoomed_region(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.orographerReadConnectionDebugState"
                " && window.orographerReadConnectionDebugState().lines.length === 1",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            before = connection_debug_state(page)
            page.evaluate(
                """
                () => {
                    const state = window.orographerReadConnectionDebugState();
                    const connectionId = state.lines[0].connectionId;
                    const firstEndpoint = state.endpoints[connectionId][0];
                    const model = Bokeh.documents[0]._all_models.get(
                        firstEndpoint.plotModelId
                    );
                    model.x_range.start = 160;
                    model.x_range.end = 200;
                    model.x_range.change.emit();
                    window.orographerUpdateReadConnectionOverlay();
                }
                """
            )
            after = connection_debug_state(page)
            assert_debug_state_matches_bokeh_scale(page, after)
        finally:
            browser.close()

    before_line = before["lines"][0]
    after_line = after["lines"][0]
    assert after_line["x1"] < before_line["x1"]
    assert_close(after_line["x2"], before_line["x2"], OVERLAY_GEOMETRY_TOLERANCE_PX)
    assert_overlay_line_matches_connection_endpoints(after)
    for endpoints in after["endpoints"].values():
        for endpoint in endpoints:
            assert_endpoint_matches_independent_range_projection(endpoint)
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_insertion_markers_recluster_after_coordinate_zoom(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_insertion_marker_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            before = insertion_marker_state(page)

            coordinate_input(page).fill("chr1:950-1120")
            page.get_by_role("button", name="Go").click()
            page.wait_for_function(
                """
                () => {
                    for (const model of Bokeh.documents[0]._all_models.values()) {
                        const data = model.data;
                        if (data && data.cluster_pos && data.cluster_pos.length === 2) {
                            return true;
                        }
                    }
                    return false;
                }
                """,
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            after = insertion_marker_state(page)
        finally:
            browser.close()

    assert before["clusterPos"] == [[1000, 1080]]
    assert before["markerLabel"] == ["65(2)"]
    assert after["clusterPos"] == [[1000], [1080]]
    assert after["markerLabel"] == ["50", "80"]
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_coordinate_input_can_zoom_to_valid_region(
    generated_plot_page: tuple[Page, list[str], list[str]],
) -> None:
    page, console_errors, page_errors = generated_plot_page
    input_box = coordinate_input(page)

    input_box.fill("chr1:120-150")
    page.get_by_role("button", name="Go").click()

    expect(input_box).to_have_value("chr1:120-150")
    expect(page.get_by_text("Coordinate not found in target region")).not_to_be_visible()
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_coordinate_input_reports_out_of_region_coordinates(
    generated_plot_page: tuple[Page, list[str], list[str]],
) -> None:
    page, console_errors, page_errors = generated_plot_page

    coordinate_input(page).fill("chr1:500-550")
    page.get_by_role("button", name="Go").click()

    expect(page.get_by_text("Coordinate not found in target region")).to_be_visible()
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_zoom_button_updates_coordinate_display(
    generated_plot_page: tuple[Page, list[str], list[str]],
) -> None:
    page, console_errors, page_errors = generated_plot_page

    page.get_by_role("button", name="+").click()

    expect(coordinate_input(page)).to_have_value(ZOOMED_COORDINATES)
    expect(page.locator(".orographer-view-size")).to_contain_text("50 bp")
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_pan_button_moves_coordinate_display(
    generated_plot_page: tuple[Page, list[str], list[str]],
) -> None:
    page, console_errors, page_errors = generated_plot_page

    page.get_by_role("button", name="+").click()
    page.get_by_role("button", name="\u203a").click()

    expect(coordinate_input(page)).to_have_value("chr1:128-178")
    expect(page.locator(".orographer-view-size")).to_contain_text("50 bp")
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_read_search_modal_opens_without_clear_button(
    generated_plot_page: tuple[Page, list[str], list[str]],
) -> None:
    page, console_errors, page_errors = generated_plot_page

    page.get_by_role("button", name="Select reads").click()
    textarea = page.locator("#modalContent textarea")
    expect(textarea).to_be_visible()

    expect(page.locator("#modalContent button").filter(has_text="Highlight")).to_have_count(1)
    expect(page.locator("#modalContent button").filter(has_text="Clear")).to_have_count(0)
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_clear_selected_button_resets_read_selection(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            page.get_by_role("button", name="Select reads").first.evaluate(
                "button => button.click()"
            )
            page.locator("#modalContent textarea").fill("readA")
            page.locator("#modalContent button").filter(has_text="Highlight").click()
            before = connector_selection_state(page)
            page.get_by_role("button", name="Clear selected").first.click()
            after = connector_selection_state(page)
        finally:
            browser.close()

    assert [item["selected"] for item in before["states"]] == [[0], [0]]
    assert [item["selected"] for item in after["states"]] == [[], []]
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_alignment_number_click_opens_read_popup_without_browser_errors(
    generated_plot_page: tuple[Page, list[str], list[str]],
) -> None:
    page, console_errors, page_errors = generated_plot_page

    position = alignment_label_screen_position(page, "readA")
    assert position is not None
    page.mouse.click(position["x"], position["y"])

    modal = page.locator("#alignmentModal")
    expect(modal).to_be_visible()
    expect(page.locator("#modalContent")).to_contain_text("Read Name:")
    expect(page.locator("#modalContent")).to_contain_text("readA")
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_chimeric_read_popup_lists_all_alignment_coordinates(
    tmp_path: Path,
    browser_name: str,
) -> None:
    plot_page = create_multi_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            position = alignment_label_screen_position(page, "readA")
            assert position is not None
            page.mouse.click(position["x"], position["y"])

            modal_content = page.locator("#modalContent")
            expect(page.locator("#alignmentModal")).to_be_visible()
            expect(modal_content).to_contain_text("All alignment coordinates:")
            expect(modal_content).to_contain_text("Alignment 1: chr1:120-180")
            expect(modal_content).to_contain_text("Alignment 2: chr2:320-380")
        finally:
            browser.close()

    assert_no_browser_errors(console_errors, page_errors)


def create_swap_region_plot_page(tmp_path: Path) -> Path:
    """Three-region plot — triggers slot-swap Select UI (N > 2 regions)."""
    seg0 = create_segment(pos=100, end=180, alignment_order=1, haplotype_tag=1)
    seg0.chrom = "chr1"
    seg0.readname = "readRegion0"
    seg1 = create_segment(pos=300, end=380, alignment_order=1, haplotype_tag=1)
    seg1.chrom = "chr2"
    seg1.readname = "readRegion1"
    seg2 = create_segment(pos=500, end=580, alignment_order=1, haplotype_tag=1)
    seg2.chrom = "chr3"
    seg2.readname = "readRegion2"
    output_file = plot_reads_bokeh(
        [
            {
                "region": Region("chr1", 100, 200, "chr1:100-200"),
                "gene_annotations": [],
                "bam_rows": [
                    {
                        "segments_by_read": {"readRegion0": [seg0]},
                        "region_type": COMPLEX_SV_REGION_TYPE,
                        "vcf_variants": [],
                        "sample_label": "region 0",
                    }
                ],
            },
            {
                "region": Region("chr2", 300, 400, "chr2:300-400"),
                "gene_annotations": [],
                "bam_rows": [
                    {
                        "segments_by_read": {"readRegion1": [seg1]},
                        "region_type": COMPLEX_SV_REGION_TYPE,
                        "vcf_variants": [],
                        "sample_label": "region 1",
                    }
                ],
            },
            {
                "region": Region("chr3", 500, 600, "chr3:500-600"),
                "gene_annotations": [],
                "bam_rows": [
                    {
                        "segments_by_read": {"readRegion2": [seg2]},
                        "region_type": COMPLEX_SV_REGION_TYPE,
                        "vcf_variants": [],
                        "sample_label": "region 2",
                    }
                ],
            },
        ],
        OutputConfig(str(tmp_path), "browser_swap_regions"),
    )
    assert output_file is not None
    return Path(output_file)


def region_select_models(page: Page) -> dict:
    """Info about slot-swap Select widgets (options are [value, label] pairs)."""
    return page.evaluate(
        """
        () => {
            const selects = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const opts = model.options;
                const val = model.value;
                if (!Array.isArray(opts) || typeof val !== "string") continue;
                if (opts.length > 0 && Array.isArray(opts[0])) {
                    selects.push({ value: val, optionCount: opts.length });
                }
            }
            return { count: selects.length, selects };
        }
        """
    )


def swap_region(page: Page, panel: str, region_idx: int) -> None:
    """Trigger the slot-swap callback for the given panel ("left" or "right")."""
    page.evaluate(
        """
        ({ panel, regionIdx }) => {
            let leftId = null;
            let rightId = null;
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const cbs = model.js_property_callbacks
                    ? model.js_property_callbacks["change:value"]
                    : null;
                if (!cbs) continue;
                cbs.forEach(function(cb) {
                    if (cb.args && cb.args.left_select) {
                        leftId = cb.args.left_select.id;
                        rightId = cb.args.right_select.id;
                    }
                });
            }
            if (!leftId) throw new Error("Swap callback not found in document");
            const targetId = panel === "left" ? leftId : rightId;
            const target = Bokeh.documents[0]._all_models.get(targetId);
            if (!target) throw new Error("Select model not found: " + targetId);
            target.setv({ value: String(regionIdx) });
        }
        """,
        {"panel": panel, "regionIdx": region_idx},
    )


def arrow_source_read_names(page: Page) -> list:
    """Read-name arrays from each arrow (read-track) source in the document."""
    return page.evaluate(
        """
        () => {
            const result = [];
            for (const model of Bokeh.documents[0]._all_models.values()) {
                const data = model.data;
                if (data && data.source_kind && data.source_kind[0] === "arrow") {
                    result.push(Array.from(data.read_name));
                }
            }
            return result;
        }
        """
    )


@pytest.mark.browser
def test_swap_region_plot_loads_select_dropdowns(
    tmp_path: Path,
    browser_name: str,
) -> None:
    """Three-region plot renders two slot-swap Select widgets without browser errors."""
    plot_page = create_swap_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            state = region_select_models(page)
        finally:
            browser.close()

    assert state["count"] == 2
    assert all(s["optionCount"] >= 2 for s in state["selects"])
    assert_no_browser_errors(console_errors, page_errors)


@pytest.mark.browser
def test_swap_region_updates_coordinate_input_and_arrow_source(
    tmp_path: Path,
    browser_name: str,
) -> None:
    """Swapping left panel to region 2 updates the coord input and arrow source reads."""
    plot_page = create_swap_region_plot_page(tmp_path)
    with serve_directory(tmp_path) as base_url, sync_playwright() as playwright:
        browser = launch_browser(playwright, browser_name)
        try:
            page = browser.new_page()
            console_errors, page_errors = attach_error_collectors(page)
            page.goto(f"{base_url}/{plot_page.name}", wait_until="networkidle")
            page.wait_for_function(
                "() => window.Bokeh && Bokeh.documents.length > 0",
                timeout=PAGE_LOAD_TIMEOUT_MS,
            )
            before_reads = arrow_source_read_names(page)
            swap_region(page, "left", 2)
            page.wait_for_timeout(600)
            after_reads = arrow_source_read_names(page)
            coord_values = page.evaluate(
                """
                () => {
                    function collectTextInputs(root, result) {
                        var inputs = root.querySelectorAll("input[type='text']");
                        for (var i = 0; i < inputs.length; i++) {
                            result.push(inputs[i].value);
                        }
                        var els = root.querySelectorAll("*");
                        for (var j = 0; j < els.length; j++) {
                            if (els[j].shadowRoot) collectTextInputs(els[j].shadowRoot, result);
                        }
                        return result;
                    }
                    return collectTextInputs(document, []);
                }
                """
            )
        finally:
            browser.close()

    all_before = [name for row in before_reads for name in row]
    assert "readRegion0" in all_before

    all_after = [name for row in after_reads for name in row]
    assert "readRegion2" in all_after

    assert any("chr3" in v for v in coord_values)
    assert_no_browser_errors(console_errors, page_errors)
