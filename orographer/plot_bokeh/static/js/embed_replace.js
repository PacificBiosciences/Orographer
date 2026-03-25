// Load JSON data from external .json.gz file and decompress in browser
const jsonFile = "__JSON_FILENAME__";
const headerBreakpoint = 825;
const centerOffsetWide = __CENTER_OFFSET_WIDE__;
const centerOffsetNarrow = 0;
const centerOffsetFullWidthAt = 1500;

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
        const data = await streamResponse.json();
        const docs_json = data.docs_json;
        const render_items = data.render_items;
        root.Bokeh.embed.embed_items(docs_json, render_items);
        let tries = 0;
        const maxTries = 80;
        function tryUpdateHeader() {
            tries += 1;
            if (updateResponsiveHeaderByModel()) return;
            if (tries < maxTries) setTimeout(tryUpdateHeader, 100);
        }
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
            "Make sure you are serving files from a web server (not opening file:// directly).</p>" +
            "<p style=\"font-size: 12px;\">Try: <code>orographer deploy --output-dir ./output</code> " +
            "or any static file server.</p>" +
            "<p style=\"font-size: 12px; margin-top: 10px;\">" +
            "Then open: <code>http://localhost:8000/__OUTPUT_BASENAME__</code></p>";
        document.body.appendChild(errorDiv);
        console.error("Error loading JSON:", error);
    }
})();
