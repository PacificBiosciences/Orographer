function parseTokens(raw) {
    const tokens = [];
    const separator = new RegExp("[\\s,]+");
    String(raw).split(separator).forEach(function (token) {
        if (token) tokens.push(token);
    });
    return tokens;
}

function setSelectedIndices(source, indices) {
    source.selected.indices = indices;
}

function clearTranscriptSelection() {
    transcript_sources.forEach(function (source) {
        setSelectedIndices(source, []);
        source.change.emit();
    });
    read_sources.forEach(function (source) {
        source.selected.indices = [];
        source.change.emit();
    });
}

function applyTranscriptSearch(tokens) {
    const wanted = {};
    tokens.forEach(function (token) {
        wanted[token] = true;
    });
    const matched = {};
    let transcriptRows = 0;

    transcript_sources.forEach(function (source) {
        const data = source.data;
        const matches = [];
        (data.transcript_id || []).forEach(function (transcriptId, rowIndex) {
            const transcriptName = data.transcript_name ? data.transcript_name[rowIndex] : "";
            if (wanted[transcriptId] || wanted[transcriptName]) {
                matches.push(rowIndex);
                matched[transcriptId] = true;
                if (transcriptName) matched[transcriptName] = true;
                transcriptRows += 1;
            }
        });
        setSelectedIndices(source, matches);
        source.change.emit();
    });

    read_sources.forEach(function (source) {
        const data = source.data;
        const matches = [];
        (data.transcript_id || []).forEach(function (transcriptId, rowIndex) {
            if (wanted[transcriptId]) {
                matches.push(rowIndex);
            }
        });
        source.selected.indices = matches;
        source.change.emit();
    });

    const missing = [];
    tokens.forEach(function (token) {
        if (!matched[token]) missing.push(token);
    });
    return { transcriptRows: transcriptRows, missing: missing };
}

const modal = document.getElementById("alignmentModal");
const content = document.getElementById("modalContent");
if (!modal || !content) {
    return;
}

content.innerHTML = "";

const title = document.createElement("div");
title.textContent = "Select Transcripts";
title.style.fontWeight = "bold";
title.style.fontSize = "16px";
title.style.marginBottom = "8px";
content.appendChild(title);

const help = document.createElement("div");
help.textContent = "Enter one or more exact transcript IDs or names, "
    + "separated by commas or spaces.";
help.style.fontSize = "12px";
help.style.color = "#555";
help.style.marginBottom = "8px";
content.appendChild(help);

const input = document.createElement("textarea");
input.placeholder = "ENST000001 transcript_name";
input.style.width = "100%";
input.style.height = "74px";
input.style.boxSizing = "border-box";
input.style.fontFamily = "Arial, sans-serif";
input.style.fontSize = "13px";
input.style.marginBottom = "8px";
content.appendChild(input);

const buttonRow = document.createElement("div");
buttonRow.style.display = "flex";
buttonRow.style.gap = "8px";
buttonRow.style.marginBottom = "8px";
content.appendChild(buttonRow);

const highlightButton = document.createElement("button");
highlightButton.textContent = "Highlight";
highlightButton.style.padding = "6px 10px";
highlightButton.style.cursor = "pointer";
buttonRow.appendChild(highlightButton);

const status = document.createElement("div");
status.style.fontSize = "12px";
status.style.color = "#555";
status.style.minHeight = "18px";
content.appendChild(status);

function highlightTranscripts() {
    const tokens = parseTokens(input.value);
    if (!tokens.length) {
        clearTranscriptSelection();
        status.textContent = "";
        return;
    }
    const result = applyTranscriptSearch(tokens);
    status.textContent = "Matched " + result.transcriptRows + " transcript glyph";
    if (result.transcriptRows != 1) status.textContent += "s";
    if (result.missing.length) {
        status.textContent += ". Missing: " + result.missing.join(", ");
    }
    modal.style.display = "none";
}

highlightButton.onclick = highlightTranscripts;
input.onkeydown = function (event) {
    if (event.key == "Enter") {
        if (!event.shiftKey) {
            event.preventDefault();
            highlightTranscripts();
        }
    }
};

if (window.fitAlignmentModalToPlot) window.fitAlignmentModalToPlot();
modal.style.display = "flex";
setTimeout(function () { input.focus(); }, 0);
