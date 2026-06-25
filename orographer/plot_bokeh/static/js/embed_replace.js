// Load JSON data from external .json.gz file and decompress in browser
const jsonFile = "__JSON_FILENAME__";
const headerBreakpoint = 825;
const centerOffsetWide = __CENTER_OFFSET_WIDE__;
const centerOffsetNarrow = 0;
const centerOffsetFullWidthAt = 1500;
const compactRepeatKey = "__orog_repeat__";
const compactDictKey = "__orog_dict__";
const compactRangeKey = "__orog_range__";
const compactRleKey = "__orog_rle__";

function hasOwn(objectValue, key) {
    return Object.prototype.hasOwnProperty.call(objectValue, key);
}

function expandCompactRepeat(encoded) {
    const compactValue = encoded[0];
    const count = encoded[1];
    const expanded = [];
    for (let index = 0; index < count; index += 1) {
        expanded.push(expandOrographerCompact(compactValue));
    }
    return expanded;
}

function expandCompactDict(encoded) {
    const values = expandOrographerCompact(encoded.values);
    const indices = expandOrographerCompact(encoded.indices);
    const expanded = [];
    for (let index = 0; index < indices.length; index += 1) {
        expanded.push(values[indices[index]]);
    }
    return expanded;
}

function expandCompactRange(encoded) {
    const start = encoded[0];
    const step = encoded[1];
    const count = encoded[2];
    const expanded = [];
    for (let index = 0; index < count; index += 1) {
        expanded.push(start + step * index);
    }
    return expanded;
}

function expandCompactRle(encoded) {
    const expanded = [];
    for (let runIndex = 0; runIndex < encoded.length; runIndex += 1) {
        const run = encoded[runIndex];
        const value = expandOrographerCompact(run[0]);
        const count = run[1];
        for (let index = 0; index < count; index += 1) {
            expanded.push(value);
        }
    }
    return expanded;
}

function expandOrographerCompact(value) {
    if (Array.isArray(value)) {
        const expanded = [];
        for (let index = 0; index < value.length; index += 1) {
            expanded.push(expandOrographerCompact(value[index]));
        }
        return expanded;
    }
    if (!value || typeof value !== "object") return value;
    if (hasOwn(value, compactRepeatKey)) return expandCompactRepeat(value[compactRepeatKey]);
    if (hasOwn(value, compactDictKey)) return expandCompactDict(value[compactDictKey]);
    if (hasOwn(value, compactRangeKey)) return expandCompactRange(value[compactRangeKey]);
    if (hasOwn(value, compactRleKey)) return expandCompactRle(value[compactRleKey]);
    const expanded = {};
    const keys = Object.keys(value);
    for (let index = 0; index < keys.length; index += 1) {
        const key = keys[index];
        expanded[key] = expandOrographerCompact(value[key]);
    }
    return expanded;
}

function getScaledCenterOffset(windowWidth) {
    if (windowWidth <= headerBreakpoint) return centerOffsetNarrow;
    if (windowWidth >= centerOffsetFullWidthAt) return centerOffsetWide;
    const span = centerOffsetFullWidthAt - headerBreakpoint;
    const ratio = (windowWidth - headerBreakpoint) / span;
    return Math.round(centerOffsetNarrow + ratio * (centerOffsetWide - centerOffsetNarrow));
}

function updateResponsiveHeaderByModel() {
    if (!root || !root.Bokeh || !root.Bokeh.documents) return false;
    const docs = root.Bokeh.documents;
    const useWide = window.innerWidth > headerBreakpoint;
    for (let i = 0; i < docs.length; i++) {
        const doc = docs[i];
        const wide = doc.get_model_by_name("orographer_original_region_wide");
        const narrow = doc.get_model_by_name("orographer_original_region_narrow_row");
        const centerOffsetSpacer = doc.get_model_by_name("orographer_center_offset_spacer");
        if (wide && narrow && centerOffsetSpacer) {
            wide.visible = useWide;
            narrow.visible = !useWide;
            centerOffsetSpacer.width = getScaledCenterOffset(window.innerWidth);
            wide.change.emit();
            narrow.change.emit();
            centerOffsetSpacer.change.emit();
            return true;
        }
    }
    return false;
}

function orographerHasTag(model, tag) {
    const tags = model.tags || [];
    for (let index = 0; index < tags.length; index += 1) {
        if (tags[index] === tag) return true;
    }
    return false;
}

function orographerTagNumber(model, prefix, fallback) {
    const tags = model.tags || [];
    for (let index = 0; index < tags.length; index += 1) {
        const tag = String(tags[index]);
        if (tag.indexOf(prefix) === 0) {
            const value = Number(tag.slice(prefix.length));
            if (Number.isFinite(value)) return value;
        }
    }
    return fallback;
}

function orographerFindBokehView(modelId) {
    const targetId = String(modelId);
    const seen = new Set();
    const queue = Object.values(root.Bokeh.index || {}).slice();
    let qi = 0;
    while (qi < queue.length) {
        const view = queue[qi];
        qi += 1;
        if (!view || typeof view !== "object" || seen.has(view)) continue;
        seen.add(view);
        if (view.model && String(view.model.id) === targetId && view.el) return view;
        const cv = view.child_views;
        if (cv instanceof Map) {
            cv.forEach(function (child) { queue.push(child); });
        } else if (cv && Array.isArray(cv)) {
            for (let j = 0; j < cv.length; j += 1) queue.push(cv[j]);
        } else if (cv && typeof cv === "object") {
            const vals = Object.values(cv);
            for (let j = 0; j < vals.length; j += 1) queue.push(vals[j]);
        }
        const vs = view.views;
        if (vs instanceof Map) {
            vs.forEach(function (child) { queue.push(child); });
        }
    }
    return null;
}


function orographerInstallButtonHold(button) {
    if (!window.orographerButtonHoldRegistry) {
        window.orographerButtonHoldRegistry = {};
    }
    if (window.orographerButtonHoldRegistry[button.id]) return true;

    const view = orographerFindBokehView(button.id);
    if (!view || !view.el) return false;
    // Attach to view.el (stable wrapper) not button_el — Bokeh may re-render button_el,
    // detaching the old node and making listeners on it unreachable.
    const element = view.el;

    const initialDelay = orographerTagNumber(button, "orographer-click-hold-delay:", 260);
    const repeatDelay = orographerTagNumber(button, "orographer-click-hold-repeat:", 95);
    const state = {
        delayTimer: null,
        repeatTimer: null,
        holding: false,
    };

    function clearTimers() {
        if (state.delayTimer) {
            window.clearTimeout(state.delayTimer);
            state.delayTimer = null;
        }
        if (state.repeatTimer) {
            window.clearInterval(state.repeatTimer);
            state.repeatTimer = null;
        }
        state.holding = false;
    }

    function repeatClick() {
        if (button.disabled) { clearTimers(); return; }
        const BtnClick = root.Bokeh.Models("ButtonClick");
        button.trigger_event(new BtnClick({origin: button}));
    }

    function startHold(event) {
        if (event && event.button != null && event.button !== 0) return;
        clearTimers();
        state.holding = true;
        state.delayTimer = window.setTimeout(function () {
            if (!state.holding) return;
            repeatClick();
            state.repeatTimer = window.setInterval(repeatClick, repeatDelay);
        }, initialDelay);
    }

    element.style.touchAction = "manipulation";
    element.addEventListener("pointerdown", startHold);
    element.addEventListener("pointerup", clearTimers);
    element.addEventListener("pointerleave", clearTimers);
    element.addEventListener("pointercancel", clearTimers);
    element.addEventListener("lostpointercapture", clearTimers);
    window.addEventListener("blur", clearTimers);
    window.orographerButtonHoldRegistry[button.id] = true;
    return true;
}

function installOrographerButtonHolds() {
    if (!root || !root.Bokeh || !root.Bokeh.documents) return false;
    let pendingCount = 0;
    let taggedCount = 0;
    const docs = root.Bokeh.documents;
    if (!docs.length) return false;
    for (let docIndex = 0; docIndex < docs.length; docIndex += 1) {
        const allModels = docs[docIndex]._all_models;
        const models = allModels ? allModels.values() : [];
        for (const model of models) {
            if (!orographerHasTag(model, "orographer-click-hold")) continue;
            taggedCount += 1;
            if (!orographerInstallButtonHold(model)) pendingCount += 1;
        }
    }
    return taggedCount > 0 && pendingCount === 0;
}

(async function () {
    try {
        const response = await fetch(jsonFile);
        if (!response.ok) {
            throw new Error(
                "Failed to load JSON file: " + response.status + " " + response.statusText
            );
        }
        const stream = response.body.pipeThrough(new DecompressionStream("gzip"));
        const streamResponse = await new Response(stream);
        const data = expandOrographerCompact(await streamResponse.json());
        const docs_json = data.docs_json;
        const render_items = data.render_items;
        await root.Bokeh.embed.embed_items(docs_json, render_items);
        let holdTries = 0;
        const maxHoldTries = 80;
        function tryInstallButtonHolds() {
            holdTries += 1;
            if (installOrographerButtonHolds()) return;
            if (holdTries < maxHoldTries) setTimeout(tryInstallButtonHolds, 100);
        }
        let tries = 0;
        const maxTries = 80;
        function tryUpdateHeader() {
            tries += 1;
            if (updateResponsiveHeaderByModel()) return;
            if (tries < maxTries) setTimeout(tryUpdateHeader, 100);
        }
        tryInstallButtonHolds();
        setTimeout(tryUpdateHeader, 200);
        window.addEventListener("resize", function () {
            updateResponsiveHeaderByModel();
        });
    } catch (error) {
        const errorDiv = document.createElement("div");
        errorDiv.style.cssText =
            "position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%); " +
            "background: #f44336; color: white; padding: 20px; border-radius: 5px; " +
            "z-index: 10000; font-family: Arial, sans-serif; max-width: 500px;";
        errorDiv.innerHTML =
            "<h3 style=\"margin-top: 0;\">Error Loading Plot Data</h3>" +
            "<p style=\"margin: 10px 0;\">" + error.message + "</p>" +
            "<p style=\"font-size: 12px; margin-top: 10px;\">" +
            "Make sure you are serving files from a web server " +
            "(not opening file:// directly).</p>" +
            "<p style=\"font-size: 12px;\">Try: <code>orographer deploy --outdir ./output</code> " +
            "or any static file server.</p>" +
            "<p style=\"font-size: 12px; margin-top: 10px;\">" +
            "Then open: <code>http://localhost:8000/__OUTPUT_BASENAME__</code></p>";
        document.body.appendChild(errorDiv);
        console.error("Error loading JSON:", error);
    }
})();
