import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath } from "node:url";

const SCRIPT_DIR = path.dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = path.resolve(SCRIPT_DIR, "../..");
const JS_DIR = path.join(REPO_ROOT, "orographer/plot_bokeh/static/js");

const tests = [];

function test(name, fn) {
    tests.push({ name, fn });
}

function runCallback(filename, context) {
    const callbackPath = path.join(JS_DIR, filename);
    const code = fs.readFileSync(callbackPath, "utf8");
    vm.runInNewContext(code, context, { filename: callbackPath });
}

function runWrappedCallback(filename, context) {
    const callbackPath = path.join(JS_DIR, filename);
    const code = fs.readFileSync(callbackPath, "utf8");
    vm.runInNewContext(`(function () {\n${code}\n}).call(globalThis);`, context, {
        filename: callbackPath,
    });
}

function runDefaultModuleCallback(filename, context, args, obj) {
    const callbackPath = path.join(JS_DIR, filename);
    const code = fs
        .readFileSync(callbackPath, "utf8")
        .replace("export default function", "globalThis.__callback = function");
    vm.runInNewContext(code, context, { filename: callbackPath });
    context.__callback(args, obj);
}

function renderer(visible = true) {
    return { visible };
}

function changeEmitter() {
    return {
        emitted: 0,
        emit() {
            this.emitted += 1;
        },
    };
}

function plain(value) {
    return JSON.parse(JSON.stringify(value));
}

function source(data, indices = []) {
    return {
        data,
        selected: { indices },
        change: changeEmitter(),
    };
}

function fakeElement(tagName, text = "") {
    const element = {
        tagName,
        textContent: text,
        innerHTML: "",
        value: "",
        placeholder: "",
        attributes: {},
        style: {},
        children: [],
        listeners: {},
        get firstChild() {
            return this.children.length === 0 ? null : this.children[0];
        },
        setAttribute(name, value) {
            this.attributes[name] = String(value);
            if (name === "id") {
                this.id = String(value);
            }
        },
        getAttribute(name) {
            return this.attributes[name];
        },
        appendChild(child) {
            this.children.push(child);
            child.parentNode = this;
        },
        contains(target) {
            let current = target;
            while (current) {
                if (current === this) return true;
                current = current.parentNode;
            }
            return false;
        },
        removeChild(child) {
            const index = this.children.indexOf(child);
            if (index !== -1) {
                this.children.splice(index, 1);
            }
            child.parentNode = null;
            return child;
        },
        addEventListener(name, callback) {
            this.listeners[name] = callback;
        },
        click() {
            if (this.listeners.click) {
                this.listeners.click({
                    preventDefault() {},
                    stopPropagation() {},
                });
            }
        },
        focus() {
            this.focused = true;
        },
    };
    return element;
}

function fakeDocument() {
    const modal = fakeElement("div");
    modal.setAttribute("id", "alignmentModal");
    const content = fakeElement("div");
    content.setAttribute("id", "modalContent");
    const body = fakeElement("body");
    const created = [];
    const listeners = {};
    function visit(element, callback) {
        callback(element);
        element.children.forEach(function (child) {
            visit(child, callback);
        });
    }
    return {
        created,
        listeners,
        modal,
        content,
        body,
        createElement(tagName) {
            const element = fakeElement(tagName);
            created.push(element);
            return element;
        },
        createElementNS(_namespace, tagName) {
            const element = fakeElement(tagName);
            created.push(element);
            return element;
        },
        getElementById(id) {
            if (id === "alignmentModal") return modal;
            if (id === "modalContent") return content;
            const match = created.find((element) => element.id === id);
            if (match) return match;
            return null;
        },
        querySelectorAll(selector) {
            const matches = [];
            if (selector.indexOf(".orographer-read-connection-overlay-line") === -1) {
                return matches;
            }
            created.forEach(function (element) {
                visit(element, function (candidate) {
                    if (
                        candidate.getAttribute
                        && candidate.getAttribute("class")
                        === "orographer-read-connection-overlay-line"
                    ) {
                        matches.push(candidate);
                    }
                });
            });
            return matches;
        },
        addEventListener(name, callback) {
            listeners[name] = callback;
        },
        querySelector() {
            return null;
        },
    };
}

test("insertion marker reclustering groups visible sites by current range", () => {
    const displaySource = {
        data: {},
        change: changeEmitter(),
    };
    const rawSource = {
        data: {
            pos: [300, 100, 108, 700],
            count: [1, 2, 3, 1],
            median_size: [42, 313, 500, 99],
            top_names: ["r4", "r1\nr2", "r3", "r5"],
            total_count: [1, 2, 3, 1],
            hp_label: ["HP1", "HP1", "HP1", "HP1"],
            chrom: ["chr1", "chr1", "chr1", "chr1"],
            y: [2, 2, 2, 2],
        },
    };

    runCallback("insertion_marker_cluster_callback.js", {
        raw_source: rawSource,
        display_source: displaySource,
        x_range: { start: 90, end: 320 },
        plot_width: 460,
        target_px: 40,
        marker_height: 18,
        marker_min_width: 24,
        marker_char_px: 7,
        Math,
        String,
    });

    assert.deepEqual(plain(displaySource.data.pos), [100, 300]);
    assert.deepEqual(plain(displaySource.data.x), [104, 300]);
    assert.deepEqual(plain(displaySource.data.median_size), [500, 42]);
    assert.deepEqual(plain(displaySource.data.marker_label), ["500(2)", "42"]);
    assert.deepEqual(plain(displaySource.data.count), [5, 1]);
    assert.deepEqual(plain(displaySource.data.cluster_pos), [[100, 108], [300]]);
    assert.equal(displaySource.change.emitted, 1);
});

test("insertion marker reclustering emits empty data when no sites are visible", () => {
    const displaySource = {
        data: {},
        change: changeEmitter(),
    };
    const rawSource = {
        data: {
            pos: [100, 200],
            count: [1, 1],
            median_size: [10, 20],
            top_names: ["a", "b"],
            total_count: [1, 1],
            hp_label: ["HP2", "HP2"],
            chrom: ["chr1", "chr1"],
            y: [1, 1],
        },
    };

    runCallback("insertion_marker_cluster_callback.js", {
        raw_source: rawSource,
        display_source: displaySource,
        x_range: { start: 500, end: 600 },
        plot_width: 400,
        target_px: 40,
        marker_height: 18,
        marker_min_width: 24,
        marker_char_px: 7,
        Math,
        String,
    });

    assert.deepEqual(plain(displaySource.data.marker_label), []);
    assert.deepEqual(plain(displaySource.data.cluster_pos), []);
    assert.equal(displaySource.change.emitted, 1);
});

test("insertion marker reclustering rounds fractional weighted medians", () => {
    const displaySource = {
        data: {},
        change: changeEmitter(),
    };
    const rawSource = {
        data: {
            pos: [100, 108],
            count: [1, 1],
            median_size: [313, 500],
            top_names: ["r1", "r2"],
            total_count: [1, 1],
            hp_label: ["HP1", "HP1"],
            chrom: ["chr1", "chr1"],
            y: [2, 2],
        },
    };

    runCallback("insertion_marker_cluster_callback.js", {
        raw_source: rawSource,
        display_source: displaySource,
        x_range: { start: 90, end: 120 },
        plot_width: 300,
        target_px: 120,
        marker_height: 18,
        marker_min_width: 24,
        marker_char_px: 7,
        Math,
        String,
    });

    assert.deepEqual(plain(displaySource.data.median_size), [407]);
    assert.deepEqual(plain(displaySource.data.marker_label), ["407(2)"]);
});

test("alignment number modal includes complex sv inclusion reason when present", () => {
    const document = fakeDocument();
    const alignmentSource = source(
        {
            read_name: ["readA"],
            alignment_number: [1],
            strand: ["Forward +"],
            coordinates: ["chr1:100-200"],
            haplotype: ["HP:1"],
            sample_label: ["sample"],
            inclusion_reason: ["10 bp INS at chr1:130"],
        },
        [0],
    );

    runWrappedCallback("number_click_callback.js", {
        source: alignmentSource,
        all_sources: [alignmentSource],
        document,
        window: {},
        String,
    });

    assert.match(document.content.innerHTML, /Primary inclusion reason:/);
    assert.match(document.content.innerHTML, /10 bp INS at chr1:130/);
    assert.equal(document.modal.style.display, "flex");
    assert.deepEqual(plain(alignmentSource.selected.indices), []);
});

test("alignment number modal lists chimeric read alignment coordinates", () => {
    const document = fakeDocument();
    const firstSource = source(
        {
            read_name: ["readA"],
            alignment_number: [2],
            strand: ["Forward +"],
            coordinates: ["chr14:105000000-105000200"],
            haplotype: ["HP:1"],
            sample_label: ["sample"],
            inclusion_reason: ["Chimeric"],
        },
        [0],
    );
    const secondSource = source(
        {
            read_name: ["readA", "readB"],
            alignment_number: [1, 1],
            strand: ["Reverse -", "Forward +"],
            coordinates: ["chr8:127000000-127000200", "chr1:1-2"],
            haplotype: ["HP:1", "HP:2"],
            sample_label: ["sample", "sample"],
            inclusion_reason: ["", ""],
        },
        [],
    );

    runWrappedCallback("number_click_callback.js", {
        source: firstSource,
        all_sources: [firstSource, secondSource],
        document,
        window: {},
        String,
        Number,
    });

    assert.match(document.content.innerHTML, /All alignment coordinates:/);
    assert.match(document.content.innerHTML, /Alignment 1: chr8:127000000-127000200/);
    assert.match(document.content.innerHTML, /Alignment 2: chr14:105000000-105000200/);
});

test("hide 1 bp callback hides all one-base renderers when active", () => {
    const oneBpRenderers = [renderer(true), renderer(true)];
    const oneBpMarkers = [renderer(true)];
    const oneBpTexts = [renderer(false)];
    const oneBpSegments = [renderer(true)];

    runCallback("hide_1bp_callback.js", {
        cb_obj: { active: true },
        one_bp_renderers: oneBpRenderers,
        one_bp_markers: oneBpMarkers,
        one_bp_texts: oneBpTexts,
        one_bp_segments: oneBpSegments,
        x_range: { start: 100, end: 500 },
        Math,
    });

    assert.deepEqual(oneBpRenderers.map((item) => item.visible), [false, false]);
    assert.equal(oneBpMarkers[0].visible, true);
    assert.equal(oneBpTexts[0].visible, false);
    assert.equal(oneBpSegments[0].visible, true);
});

test("hide 1 bp callback restores text mode for narrow ranges", () => {
    const oneBpRenderers = [renderer(false)];
    const oneBpMarkers = [renderer(true)];
    const oneBpTexts = [renderer(false)];
    const oneBpSegments = [renderer(false)];

    runCallback("hide_1bp_callback.js", {
        cb_obj: { active: false },
        one_bp_renderers: oneBpRenderers,
        one_bp_markers: oneBpMarkers,
        one_bp_texts: oneBpTexts,
        one_bp_segments: oneBpSegments,
        x_range: { start: 100, end: 500 },
        Math,
    });

    assert.equal(oneBpRenderers[0].visible, true);
    assert.equal(oneBpMarkers[0].visible, false);
    assert.equal(oneBpTexts[0].visible, true);
    assert.equal(oneBpSegments[0].visible, true);
});

test("hide alignment numbers callback follows checkbox state", () => {
    const labels = [renderer(true), renderer(true)];

    runCallback("hide_alignment_numbers_callback.js", {
        cb_obj: { active: true },
        alignment_label_renderers: labels,
    });

    assert.deepEqual(labels.map((item) => item.visible), [false, false]);
});

test("hide phaseset markers callback follows checkbox state", () => {
    const markers = [renderer(true), renderer(true)];

    runCallback("hide_phase_set_markers_callback.js", {
        cb_obj: { active: true },
        phase_set_marker_renderers: markers,
    });

    assert.deepEqual(markers.map((item) => item.visible), [false, false]);
});

test("alignment label selection callback hides unselected read labels", () => {
    const first = source({ read_name: ["readA", "readB"] }, [0]);
    const connector = source({ read_name: ["readA"] }, [0]);
    const labels = source(
        {
            read_name: ["readA", "readB", "readC"],
            label_alpha: [0.8, 0.8, 0.8],
        },
        [],
    );

    runCallback("alignment_label_selection_callback.js", {
        selection_sources: [first, connector],
        alignment_label_sources: [labels],
        Math,
    });

    assert.deepEqual(plain(labels.data.label_alpha), [0.8, 0, 0]);
    assert.equal(labels.change.emitted, 1);
});

test("alignment label selection callback restores labels with no selection", () => {
    const first = source({ read_name: ["readA", "readB"] }, []);
    const labels = source(
        {
            read_name: ["readA", "readB"],
            label_alpha: [0, 0],
        },
        [],
    );

    runCallback("alignment_label_selection_callback.js", {
        selection_sources: [first],
        alignment_label_sources: [labels],
        Math,
    });

    assert.deepEqual(plain(labels.data.label_alpha), [0.8, 0.8]);
    assert.equal(labels.change.emitted, 1);
});

test("read filter callback hides reads without requested evidence", () => {
    const arrows = source(
        {
            read_name: ["readA", "readSolo"],
            source_kind: ["arrow", "arrow"],
            has_split_alignment: [true, false],
            has_multiregion_connection: [true, false],
            read_filter_alpha: [1, 1],
            arrow_line_alpha: [0.5, 0.5],
            arrowhead_alpha: [0.65, 0.65],
            arrow_selected_alpha: [1, 1],
            arrow_nonselected_alpha: [0.12, 0.12],
        },
        [],
    );
    const labels = source(
        {
            read_name: ["readA", "readSolo"],
            has_split_alignment: [true, false],
            has_multiregion_connection: [true, false],
            read_filter_alpha: [1, 1],
            label_alpha: [0.8, 0.8],
        },
        [],
    );
    const window = {
        Bokeh: {
            documents: [
                {
                    _all_models: new Map([
                        ["arrows", arrows],
                        ["labels", labels],
                    ]),
                },
            ],
        },
        orographerUpdateReadConnectionOverlay() {
            this.overlayUpdated = true;
        },
    };

    runWrappedCallback("read_filter_callback.js", {
        hide_non_split_checkbox: { active: true },
        hide_non_multiregion_checkbox: { active: false },
        window,
        Set,
        Map,
    });

    assert.deepEqual(plain(arrows.data.read_filter_alpha), [1, 0]);
    assert.deepEqual(plain(arrows.data.arrow_line_alpha), [0.5, 0]);
    assert.deepEqual(plain(labels.data.label_alpha), [0.8, 0]);
    assert.equal(arrows.change.emitted, 1);
    assert.equal(labels.change.emitted, 1);
    assert.equal(window.overlayUpdated, true);
});

test("coordinate go callback accepts comma-formatted coordinates in range", () => {
    const xRange = { start: 0, end: 0, change: changeEmitter() };
    const errorDiv = { text: "old", styles: {} };

    runCallback("coord_go_callback.js", {
        coord_input: { value: "chr1:1,100-1,250" },
        x_range: xRange,
        error_div: errorDiv,
        orig_start: 1000,
        orig_end: 2000,
        Math,
        String,
        RegExp,
        parseInt,
        isNaN,
    });

    assert.equal(xRange.start, 1100);
    assert.equal(xRange.end, 1250);
    assert.equal(xRange.change.emitted, 1);
    assert.equal(errorDiv.text, "");
    assert.equal(errorDiv.styles.height, "0");
});

test("coordinate go callback expands tiny windows to at least ten bases", () => {
    const xRange = { start: 0, end: 0, change: changeEmitter() };
    const errorDiv = { text: "", styles: {} };

    runCallback("coord_go_callback.js", {
        coord_input: { value: "1005-1006" },
        x_range: xRange,
        error_div: errorDiv,
        orig_start: 1000,
        orig_end: 2000,
        Math,
        String,
        RegExp,
        parseInt,
        isNaN,
    });

    assert.equal(xRange.start, 1000);
    assert.equal(xRange.end, 1011);
    assert.equal(xRange.change.emitted, 1);
    assert.equal(errorDiv.text, "");
});

test("coordinate go callback rejects coordinates outside the original region", () => {
    const xRange = { start: 1100, end: 1200 };
    const errorDiv = { text: "", styles: {} };

    runCallback("coord_go_callback.js", {
        coord_input: { value: "3000-3010" },
        x_range: xRange,
        error_div: errorDiv,
        orig_start: 1000,
        orig_end: 2000,
        Math,
        String,
        RegExp,
        parseInt,
        isNaN,
    });

    assert.equal(xRange.start, 1100);
    assert.equal(xRange.end, 1200);
    assert.equal(errorDiv.text, "Coordinate not found in target region");
    assert.equal(errorDiv.styles.height, "auto");
});

test("range reset callback resets ranges and clears selections", () => {
    const sourceA = {
        selected: { indices: [1, 2] },
        _lastNonEmpty: [2],
        change: changeEmitter(),
    };
    const sourceB = {
        selected: { indices: [3] },
        change: changeEmitter(),
    };
    const yRanges = [
        { start: 0, end: 0 },
        { start: 0, end: 0 },
    ];

    runCallback("range_reset_callback.js", {
        x_range: { start: 0, end: 0 },
        x_start: 100,
        x_end: 200,
        y_ranges: yRanges,
        y_bounds: [[10, 1], [20, 2]],
        all_sources: [sourceA, sourceB],
    });

    assert.deepEqual(yRanges, [
        { start: 10, end: 1 },
        { start: 20, end: 2 },
    ]);
    assert.deepEqual(plain(sourceA.selected.indices), []);
    assert.deepEqual(plain(sourceA._lastNonEmpty), []);
    assert.deepEqual(plain(sourceB.selected.indices), []);
    assert.equal(sourceA.change.emitted, 1);
    assert.equal(sourceB.change.emitted, 1);
});

test("variant lod callback switches between marker and text renderers", () => {
    const markerRenderers = [renderer(true), renderer(true)];
    const textRenderers = [renderer(false), renderer(false)];

    runDefaultModuleCallback(
        "variant_lod_callback.js",
        { Math },
        {
            x_range: { start: 100, end: 500 },
            marker_renderers: markerRenderers,
            text_renderers: textRenderers,
            one_bp_renderers: [],
            hide_1bp_checkbox: null,
        },
    );

    assert.deepEqual(markerRenderers.map((item) => item.visible), [false, false]);
    assert.deepEqual(textRenderers.map((item) => item.visible), [true, true]);
});

test("variant lod callback keeps one-base renderers hidden when checkbox is active", () => {
    const markerRenderers = [renderer(false)];
    const textRenderers = [renderer(true)];
    const oneBpRenderers = [renderer(true), renderer(true)];

    runDefaultModuleCallback(
        "variant_lod_callback.js",
        { Math },
        {
            x_range: { start: 100, end: 5000 },
            marker_renderers: markerRenderers,
            text_renderers: textRenderers,
            one_bp_renderers: oneBpRenderers,
            hide_1bp_checkbox: { active: true },
        },
    );

    assert.equal(markerRenderers[0].visible, true);
    assert.equal(textRenderers[0].visible, false);
    assert.deepEqual(oneBpRenderers.map((item) => item.visible), [false, false]);
});

test("view callback formats coordinate and small view size", () => {
    const viewSizeDiv = { text: "", change: changeEmitter() };
    const coordInput = { value: "" };
    const errorDiv = { text: "old", styles: {} };

    runDefaultModuleCallback(
        "view_callback.js",
        { RegExp, Math },
        {
            chrom: "chr1",
            coord_input: coordInput,
            error_div: errorDiv,
            view_size_div: viewSizeDiv,
        },
        { start: 1100.2, end: 1200.7 },
    );

    assert.equal(coordInput.value, "chr1:1,100-1,201");
    assert.equal(errorDiv.text, "");
    assert.equal(errorDiv.styles.height, "0");
    assert.equal(viewSizeDiv.text, "101 bp");
    assert.equal(viewSizeDiv.change.emitted, 1);
});

test("view callback formats kilobase and megabase spans", () => {
    const kbDiv = { text: "", change: changeEmitter() };
    const mbDiv = { text: "", change: changeEmitter() };

    runDefaultModuleCallback(
        "view_callback.js",
        { RegExp, Math },
        {
            chrom: "chr1",
            coord_input: { value: "" },
            error_div: { text: "", styles: {} },
            view_size_div: kbDiv,
        },
        { start: 1000, end: 3500 },
    );
    runDefaultModuleCallback(
        "view_callback.js",
        { RegExp, Math },
        {
            chrom: "chr1",
            coord_input: { value: "" },
            error_div: { text: "", styles: {} },
            view_size_div: mbDiv,
        },
        { start: 1000, end: 1_501_000 },
    );

    assert.equal(kbDiv.text, "2.5 kb");
    assert.equal(mbDiv.text, "1.5 Mb");
});

test("zoom buttons callback zooms in around center and emits range change", () => {
    const xRange = { start: 100, end: 300, change: changeEmitter() };

    runCallback("zoom_buttons_callback.js", {
        x_range: xRange,
        factor: 0.5,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(xRange.start, 150);
    assert.equal(xRange.end, 250);
    assert.equal(xRange.change.emitted, 1);
});

test("zoom buttons callback clamps zoom out to original region bounds", () => {
    const xRange = { start: 100, end: 900, change: changeEmitter() };

    runCallback("zoom_buttons_callback.js", {
        x_range: xRange,
        factor: 2,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(xRange.start, 0);
    assert.equal(xRange.end, 1000);
    assert.equal(xRange.change.emitted, 1);
});

test("pan buttons callback shifts the view by five percent", () => {
    const xRange = { start: 100, end: 300, change: changeEmitter() };

    runCallback("pan_buttons_callback.js", {
        x_range: xRange,
        direction: 1,
        fraction: 0.05,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(xRange.start, 110);
    assert.equal(xRange.end, 310);
    assert.equal(xRange.change.emitted, 1);
});

test("pan buttons callback clamps at original region bounds", () => {
    const xRange = { start: 10, end: 210, change: changeEmitter() };

    runCallback("pan_buttons_callback.js", {
        x_range: xRange,
        direction: -1,
        fraction: 0.05,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(xRange.start, 0);
    assert.equal(xRange.end, 200);
    assert.equal(xRange.change.emitted, 1);
});

test("pan buttons state callback disables unavailable directions", () => {
    const leftButton = { disabled: false };
    const rightButton = { disabled: false };
    const zoomOutButton = { disabled: false };
    const zoomInButton = { disabled: false };

    runCallback("pan_buttons_state_callback.js", {
        x_range: { start: 0, end: 1000 },
        left_button: leftButton,
        right_button: rightButton,
        zoom_out_button: zoomOutButton,
        zoom_in_button: zoomInButton,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(leftButton.disabled, true);
    assert.equal(rightButton.disabled, true);
    assert.equal(zoomOutButton.disabled, true);
    assert.equal(zoomInButton.disabled, false);

    runCallback("pan_buttons_state_callback.js", {
        x_range: { start: 100, end: 900 },
        left_button: leftButton,
        right_button: rightButton,
        zoom_out_button: zoomOutButton,
        zoom_in_button: zoomInButton,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(leftButton.disabled, false);
    assert.equal(rightButton.disabled, false);
    assert.equal(zoomOutButton.disabled, false);
    assert.equal(zoomInButton.disabled, false);

    runCallback("pan_buttons_state_callback.js", {
        x_range: { start: 0, end: 500 },
        left_button: leftButton,
        right_button: rightButton,
        zoom_out_button: zoomOutButton,
        zoom_in_button: zoomInButton,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(leftButton.disabled, true);
    assert.equal(rightButton.disabled, false);
    assert.equal(zoomOutButton.disabled, false);
    assert.equal(zoomInButton.disabled, false);

    runCallback("pan_buttons_state_callback.js", {
        x_range: { start: 100, end: 110 },
        left_button: leftButton,
        right_button: rightButton,
        zoom_out_button: zoomOutButton,
        zoom_in_button: zoomInButton,
        orig_start: 0,
        orig_end: 1000,
        Math,
        String,
    });

    assert.equal(zoomOutButton.disabled, false);
    assert.equal(zoomInButton.disabled, true);
});

test("y scroll callback pans and clamps to the full y extent", () => {
    const yRange = { start: 20, end: 80, change: changeEmitter() };

    runCallback("y_scroll_callback.js", {
        cb_obj: { delta: 1 },
        y_range: yRange,
        y_min: 0,
        y_max: 100,
        Math,
    });

    assert.equal(yRange.start, 8);
    assert.equal(yRange.end, 68);
    assert.equal(yRange.change.emitted, 1);

    runCallback("y_scroll_callback.js", {
        cb_obj: { delta: -1 },
        y_range: yRange,
        y_min: 0,
        y_max: 100,
        Math,
    });
    runCallback("y_scroll_callback.js", {
        cb_obj: { delta: -1 },
        y_range: yRange,
        y_min: 0,
        y_max: 100,
        Math,
    });
    runCallback("y_scroll_callback.js", {
        cb_obj: { delta: -1 },
        y_range: yRange,
        y_min: 0,
        y_max: 100,
        Math,
    });

    assert.equal(yRange.start, 40);
    assert.equal(yRange.end, 100);
});

test("shared and arrow reset callbacks clear cached selections", () => {
    const first = source({}, [1]);
    const second = source({}, [2, 3]);
    first._lastNonEmpty = [1];
    second._lastNonEmpty = [3];

    runCallback("shared_reset_callback.js", {
        all_sources: [first, second],
    });

    assert.deepEqual(plain(first.selected.indices), []);
    assert.deepEqual(plain(second.selected.indices), []);
    assert.deepEqual(plain(first._lastNonEmpty), []);
    assert.deepEqual(plain(second._lastNonEmpty), []);
    assert.equal(first.change.emitted, 1);
    assert.equal(second.change.emitted, 1);

    first.selected.indices = [4];
    first._lastNonEmpty = [4];
    runCallback("arrow_reset_callback.js", { source: first });

    assert.deepEqual(plain(first.selected.indices), []);
    assert.deepEqual(plain(first._lastNonEmpty), []);
    assert.equal(first.change.emitted, 2);
});

test("arrow tap callback selects and deselects all segments for a read", () => {
    const arrowSource = source(
        {
            x0: [100, 200, 500],
            x1: [150, 250, 550],
            y: [1, 1, 2],
            read_name: ["readA", "readA", "readB"],
        },
        [],
    );

    runWrappedCallback("arrow_tap_callback.js", {
        source: arrowSource,
        cb_obj: { x: 125, y: 1.05 },
        Math,
        Infinity,
    });

    assert.deepEqual(plain(arrowSource.selected.indices), [0, 1]);
    assert.deepEqual(plain(arrowSource._lastNonEmpty), [0, 1]);
    assert.equal(arrowSource.change.emitted, 1);

    runWrappedCallback("arrow_tap_callback.js", {
        source: arrowSource,
        cb_obj: { x: 125, y: 1.05 },
        Math,
        Infinity,
    });

    assert.deepEqual(plain(arrowSource.selected.indices), []);
    assert.deepEqual(plain(arrowSource._lastNonEmpty), []);
    assert.equal(arrowSource.change.emitted, 2);
});

test("arrow tap callback preserves selection when tap misses every segment", () => {
    const arrowSource = source(
        {
            x0: [100],
            x1: [150],
            y: [1],
            read_name: ["readA"],
        },
        [0],
    );
    arrowSource._lastNonEmpty = [0];

    runWrappedCallback("arrow_tap_callback.js", {
        source: arrowSource,
        cb_obj: { x: 900, y: 10 },
        Math,
        Infinity,
    });

    assert.deepEqual(arrowSource.selected.indices, [0]);
    assert.deepEqual(arrowSource._lastNonEmpty, [0]);
    assert.equal(arrowSource.change.emitted, 0);
});

test("multi-region arrow tap selects matching read across sources", () => {
    const primary = source(
        {
            x0: [100, 200],
            x1: [150, 250],
            y: [1, 2],
            read_name: ["readA", "readB"],
        },
        [],
    );
    const secondary = source(
        {
            x0: [300, 400, 500],
            x1: [350, 450, 550],
            y: [1, 2, 3],
            read_name: ["readC", "readA", "readA"],
        },
        [0],
    );

    runWrappedCallback("arrow_tap_multi_region_callback.js", {
        source: primary,
        all_sources: [primary, secondary],
        cb_obj: { x: 125, y: 1 },
        window: {},
        Math,
        Infinity,
    });

    assert.deepEqual(plain(primary.selected.indices), [0]);
    assert.deepEqual(plain(secondary.selected.indices), [0, 1, 2]);
    assert.equal(primary.change.emitted, 1);
    assert.equal(secondary.change.emitted, 1);

    runWrappedCallback("arrow_tap_multi_region_callback.js", {
        source: primary,
        all_sources: [primary, secondary],
        cb_obj: { x: 125, y: 1 },
        window: {},
        Math,
        Infinity,
    });

    assert.deepEqual(plain(primary.selected.indices), []);
    assert.deepEqual(plain(secondary.selected.indices), [0]);
});

test("multi-region arrow tap selects connector sources for matching read", () => {
    const primary = source(
        {
            x0: [100],
            x1: [150],
            y: [1],
            read_name: ["readA"],
            source_kind: ["arrow"],
        },
        [],
    );
    const connector = source(
        {
            read_name: ["readA", "readB"],
            source_kind: ["connector", "connector"],
        },
        [],
    );

    runWrappedCallback("arrow_tap_multi_region_callback.js", {
        source: primary,
        all_sources: [primary, connector],
        cb_obj: { x: 125, y: 1 },
        window: {},
        Math,
        Infinity,
    });

    assert.deepEqual(plain(primary.selected.indices), [0]);
    assert.deepEqual(plain(connector.selected.indices), [0]);
    assert.equal(connector.change.emitted, 1);
});

test("multi-region arrow tap selects read when connector line is clicked", () => {
    const primary = source(
        {
            x0: [100, 300],
            x1: [150, 350],
            y: [1, 2],
            read_name: ["readA", "readB"],
            source_kind: ["arrow", "arrow"],
        },
        [],
    );
    const connector = source(
        {
            xs: [[150, 150, 300, 300]],
            ys: [[1, 1.3, 1.3, 2]],
            read_name: ["readB"],
            source_kind: ["connector"],
            connection_id: ["connection-1"],
        },
        [],
    );

    runWrappedCallback("arrow_tap_multi_region_callback.js", {
        source: primary,
        all_sources: [primary, connector],
        cb_obj: { x: 225, y: 1.3 },
        window: {},
        Math,
        Infinity,
    });

    assert.deepEqual(plain(primary.selected.indices), [1]);
    assert.deepEqual(plain(connector.selected.indices), [0]);
    assert.equal(primary.change.emitted, 1);
    assert.equal(connector.change.emitted, 1);
});

test("read search modal highlights matching reads and reports missing names", () => {
    const first = source({ read_name: ["readA", "readB", "readA"] }, []);
    const second = source({ read_name: ["readC", "readA"] }, []);
    const document = fakeDocument();

    runWrappedCallback("read_search_modal_callback.js", {
        all_sources: [first, second],
        document,
        window: {},
        Math,
        Object,
        String,
        setTimeout(fn) {
            fn();
        },
    });

    const input = document.created.find((element) => element.tagName === "textarea");
    const buttons = document.created.filter((element) => element.tagName === "button");
    const status = document.content.children[document.content.children.length - 1];
    input.value = "readA, missing";
    buttons[0].onclick();

    assert.deepEqual(plain(first.selected.indices), [0, 2]);
    assert.deepEqual(plain(second.selected.indices), [1]);
    assert.equal(first.change.emitted, 1);
    assert.equal(second.change.emitted, 1);
    assert.equal(status.textContent, "Matched 1 read across 3 segments. Missing: missing");
    assert.equal(document.modal.style.display, "none");
});

test("read search modal selects connectors without counting them as segments", () => {
    const arrowSource = source(
        {
            read_name: ["readA"],
            source_kind: ["arrow"],
        },
        [],
    );
    const connectorSource = source(
        {
            read_name: ["readA", "readA"],
            source_kind: ["connector", "connector"],
        },
        [],
    );
    const document = fakeDocument();

    runWrappedCallback("read_search_modal_callback.js", {
        all_sources: [arrowSource, connectorSource],
        document,
        window: {},
        Math,
        Object,
        String,
        setTimeout(fn) {
            fn();
        },
    });

    const input = document.created.find((element) => element.tagName === "textarea");
    const highlightButton = document.created.filter((element) => element.tagName === "button")[0];
    const status = document.content.children[document.content.children.length - 1];
    input.value = "readA";
    highlightButton.onclick();

    assert.deepEqual(plain(arrowSource.selected.indices), [0]);
    assert.deepEqual(plain(connectorSource.selected.indices), [0, 1]);
    assert.equal(status.textContent, "Matched 1 read across 1 segment");
});

test("read search modal omits clear button", () => {
    const first = source({ read_name: ["readA"] }, [0]);
    first._lastNonEmpty = [0];
    const document = fakeDocument();

    runWrappedCallback("read_search_modal_callback.js", {
        all_sources: [first],
        document,
        window: {},
        Math,
        Object,
        String,
        setTimeout(fn) {
            fn();
        },
    });

    const buttons = document.created.filter((element) => element.tagName === "button");

    assert.equal(buttons.length, 1);
    assert.equal(buttons[0].textContent, "Highlight");
});

test("clear read selection callback resets selections", () => {
    const first = source({ read_name: ["readA"] }, [0]);
    first._lastNonEmpty = [0];
    const clearButton = { disabled: false };
    const window = {
        orographerUpdateReadConnectionOverlay() {
            this.overlayUpdated = true;
        },
    };

    runCallback("clear_read_selection_callback.js", {
        all_sources: [first],
        clear_button: clearButton,
        window,
    });

    assert.deepEqual(plain(first.selected.indices), []);
    assert.deepEqual(plain(first._lastNonEmpty), []);
    assert.equal(first.change.emitted, 1);
    assert.equal(clearButton.disabled, true);
    assert.equal(window.overlayUpdated, true);
});

test("clear read selection state callback follows selected reads", () => {
    const first = source({ read_name: ["readA"] }, []);
    const second = source({ read_name: ["readB"] }, [0]);
    const clearButton = { disabled: true };

    runCallback("clear_read_selection_state_callback.js", {
        all_sources: [first, second],
        clear_button: clearButton,
    });

    assert.equal(clearButton.disabled, false);

    second.selected.indices = [];
    runCallback("clear_read_selection_state_callback.js", {
        all_sources: [first, second],
        clear_button: clearButton,
    });

    assert.equal(clearButton.disabled, true);
});

test("modal close button consumes click before it can reach plot handlers", () => {
    const document = fakeDocument();
    const closeButton = document.createElement("span");
    closeButton.setAttribute("id", "closeModal");
    const title = document.createElement("h2");
    title.setAttribute("id", "alignmentModalTitle");
    const dialog = document.createElement("div");
    dialog.setAttribute("id", "alignmentModalDialog");
    const wrapper = document.createElement("div");
    wrapper.setAttribute("id", "alignmentModalWrapper");
    document.modal.style.display = "flex";
    document.modal.appendChild(wrapper);
    wrapper.appendChild(dialog);
    dialog.appendChild(closeButton);
    const event = {
        prevented: false,
        stopped: false,
        immediateStopped: false,
        preventDefault() {
            this.prevented = true;
        },
        stopPropagation() {
            this.stopped = true;
        },
        stopImmediatePropagation() {
            this.immediateStopped = true;
        },
    };
    const window = {
        innerHeight: 800,
        addEventListener() {},
    };

    runCallback("modal.js", {
        document,
        window,
        Math,
    });
    closeButton.listeners.click(event);

    assert.equal(document.modal.style.display, "none");
    assert.equal(event.prevented, true);
    assert.equal(event.stopped, true);
    assert.equal(event.immediateStopped, true);
});

test("modal dialog does not prevent default text selection", () => {
    const document = fakeDocument();
    const dialog = document.createElement("div");
    dialog.setAttribute("id", "alignmentModalDialog");
    document.modal.style.display = "flex";
    document.modal.appendChild(dialog);
    const event = {
        prevented: false,
        stopped: false,
        immediateStopped: false,
        preventDefault() {
            this.prevented = true;
        },
        stopPropagation() {
            this.stopped = true;
        },
        stopImmediatePropagation() {
            this.immediateStopped = true;
        },
    };
    const window = {
        innerHeight: 800,
        addEventListener() {},
    };

    runCallback("modal.js", {
        document,
        window,
        Math,
    });
    dialog.listeners.mousedown(event);

    assert.equal(event.prevented, false);
    assert.equal(event.stopped, true);
    assert.equal(event.immediateStopped, true);
});

test("read connection overlay draws selected endpoint pairs", () => {
    const document = fakeDocument();
    const firstSource = source(
        {
            marker_x: [100],
            y: [2],
            read_name: ["readA"],
            source_kind: ["connector"],
            connection_id: ["connection-1"],
            haplotype_transition: [false],
        },
        [0],
    );
    const secondSource = source(
        {
            marker_x: [300],
            y: [4],
            read_name: ["readA"],
            connection_id: ["connection-1"],
            haplotype_transition: [true],
        },
        [0],
    );
    const markerView = function (dataSource, offsetLeft, offsetTop) {
        return {
            model: {
                name: "orographer_read_connector_marker",
                data_source: dataSource,
            },
            plot_view: {
                frame: {
                    x_scale: { compute: (value) => value },
                    y_scale: { compute: (value) => value },
                },
                canvas_view: {
                    events_el: {
                        getBoundingClientRect() {
                            return { left: offsetLeft, top: offsetTop };
                        },
                    },
                },
            },
        };
    };
    const window = {
        orographerSelectableSources: [firstSource, secondSource],
        Bokeh: {
            index: {
                root: {
                    child_views: {
                        first: markerView(firstSource, 10, 20),
                        second: markerView(secondSource, 30, 40),
                    },
                },
            },
        },
        addEventListener() {},
    };

    runCallback("read_connection_overlay.js", {
        document,
        window,
        Set,
        Object,
        String,
    });
    window.orographerUpdateReadConnectionOverlay();

    const overlay = document.getElementById("orographerReadConnectionOverlay");
    assert.equal(overlay.children.length, 1);
    assert.equal(overlay.style.zIndex, "900");
    const group = overlay.children[0];
    const path = group.children[0];
    const hitPath = group.children[1];
    assert.equal(group.getAttribute("data-connection-id"), "connection-1");
    assert.equal(group.getAttribute("data-read-name"), "readA");
    assert.equal(path.getAttribute("data-x1"), "110");
    assert.equal(path.getAttribute("data-y1"), "22");
    assert.equal(path.getAttribute("data-x2"), "330");
    assert.equal(path.getAttribute("data-y2"), "44");
    assert.equal(path.getAttribute("stroke-dasharray"), "6 4");
    assert.equal(path.getAttribute("stroke"), "#92400e");
    assert.match(path.getAttribute("d"), /^M 110 22 L /);
    assert.equal(path.style.pointerEvents, "none");
    assert.equal(hitPath.getAttribute("class"), "orographer-read-connection-overlay-hit-area");
    assert.equal(hitPath.getAttribute("stroke-width"), "14");
    assert.equal(hitPath.style.pointerEvents, "stroke");
});

test("read connection overlay path click toggles selected read across sources", () => {
    const document = fakeDocument();
    const arrowSource = source(
        {
            read_name: ["readA", "readB"],
            source_kind: ["arrow", "arrow"],
        },
        [],
    );
    const firstSource = source(
        {
            marker_x: [100],
            y: [2],
            read_name: ["readA"],
            source_kind: ["connector"],
            connection_id: ["connection-1"],
            haplotype_transition: [false],
        },
        [],
    );
    const secondSource = source(
        {
            marker_x: [300],
            y: [4],
            read_name: ["readA"],
            source_kind: ["connector"],
            connection_id: ["connection-1"],
            haplotype_transition: [false],
        },
        [],
    );
    const markerView = function (dataSource, offsetLeft, offsetTop) {
        return {
            model: {
                name: "orographer_read_connector_marker",
                data_source: dataSource,
            },
            plot_view: {
                frame: {
                    x_scale: { compute: (value) => value },
                    y_scale: { compute: (value) => value },
                },
                canvas_view: {
                    events_el: {
                        getBoundingClientRect() {
                            return { left: offsetLeft, top: offsetTop };
                        },
                    },
                },
            },
        };
    };
    const window = {
        orographerSelectableSources: [arrowSource, firstSource, secondSource],
        Bokeh: {
            index: {
                root: {
                    child_views: {
                        first: markerView(firstSource, 10, 20),
                        second: markerView(secondSource, 30, 40),
                    },
                },
            },
        },
        addEventListener() {},
    };

    runCallback("read_connection_overlay.js", {
        document,
        window,
        Set,
        Object,
        String,
    });
    window.orographerUpdateReadConnectionOverlay();
    document.getElementById("orographerReadConnectionOverlay").children[0].click();

    assert.deepEqual(plain(arrowSource.selected.indices), [0]);
    assert.deepEqual(plain(firstSource.selected.indices), [0]);
    assert.deepEqual(plain(secondSource.selected.indices), [0]);
    assert.equal(arrowSource.change.emitted, 1);
    assert.equal(firstSource.change.emitted, 1);
    assert.equal(secondSource.change.emitted, 1);
});

test("read connection overlay ignores clicks inside visible modal", () => {
    const document = fakeDocument();
    const closeButton = document.createElement("span");
    closeButton.setAttribute("id", "closeModal");
    document.modal.style.display = "flex";
    document.modal.appendChild(closeButton);
    const arrowSource = source(
        {
            read_name: ["readA"],
            source_kind: ["arrow"],
        },
        [],
    );
    const firstSource = source(
        {
            marker_x: [100],
            y: [2],
            read_name: ["readA"],
            source_kind: ["connector"],
            connection_id: ["connection-1"],
            haplotype_transition: [false],
        },
        [],
    );
    const secondSource = source(
        {
            marker_x: [300],
            y: [4],
            read_name: ["readA"],
            source_kind: ["connector"],
            connection_id: ["connection-1"],
            haplotype_transition: [false],
        },
        [],
    );
    const markerView = function (dataSource, offsetLeft, offsetTop) {
        return {
            model: {
                name: "orographer_read_connector_marker",
                data_source: dataSource,
            },
            plot_view: {
                frame: {
                    x_scale: { compute: (value) => value },
                    y_scale: { compute: (value) => value },
                },
                canvas_view: {
                    events_el: {
                        getBoundingClientRect() {
                            return { left: offsetLeft, top: offsetTop };
                        },
                    },
                },
            },
        };
    };
    const window = {
        listeners: {},
        orographerSelectableSources: [arrowSource, firstSource, secondSource],
        Bokeh: {
            index: {
                root: {
                    child_views: {
                        first: markerView(firstSource, 10, 20),
                        second: markerView(secondSource, 30, 40),
                    },
                },
            },
        },
        addEventListener(name, callback) {
            this.listeners[name] = callback;
        },
    };

    runCallback("read_connection_overlay.js", {
        document,
        window,
        Set,
        Object,
        String,
    });
    window.orographerUpdateReadConnectionOverlay();
    const event = {
        clientX: 220,
        clientY: 33,
        target: closeButton,
        preventDefault() {
            this.prevented = true;
        },
        stopPropagation() {
            this.stopped = true;
        },
        stopImmediatePropagation() {
            this.immediateStopped = true;
        },
    };
    window.listeners.click(event);

    assert.deepEqual(plain(arrowSource.selected.indices), []);
    assert.deepEqual(plain(firstSource.selected.indices), []);
    assert.deepEqual(plain(secondSource.selected.indices), []);
    assert.equal(event.prevented, undefined);
    assert.equal(event.stopped, undefined);
    assert.equal(event.immediateStopped, undefined);
});

let failures = 0;
for (const { name, fn } of tests) {
    try {
        fn();
        console.log(`ok - ${name}`);
    } catch (error) {
        failures += 1;
        console.error(`not ok - ${name}`);
        console.error(error);
    }
}

if (failures > 0) {
    process.exitCode = 1;
}

console.log(`${tests.length - failures}/${tests.length} JS unit tests passed`);
