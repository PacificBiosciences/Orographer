var orographerClearSelections = function (sources) {
    var sourcePos = 0;
    while (Math.sign(sources.length - sourcePos) === 1) {
        var source = sources[sourcePos];
        source.selected.indices = [];
        source._lastNonEmpty = [];
        source.change.emit();
        sourcePos += 1;
    }
};

var orographerParseReadSearchTokens = function (raw) {
    var tokens = [];
    var token = "";
    var pos = 0;
    raw = String(raw);
    while (Math.sign(raw.length - pos) === 1) {
        var ch = raw.charAt(pos);
        var isDelim = false;
        if (ch === ",") isDelim = true;
        if (ch === " ") isDelim = true;
        if (ch === "\t") isDelim = true;
        if (ch === "\n") isDelim = true;
        if (ch === "\r") isDelim = true;
        if (isDelim) {
            if (token.length !== 0) tokens.push(token);
            token = "";
        } else {
            token += ch;
        }
        pos += 1;
    }
    if (token.length !== 0) tokens.push(token);
    return tokens;
};

var orographerApplyReadSearch = function (sources, tokens) {
    var wanted = {};
    var tokenPos = 0;
    while (Math.sign(tokens.length - tokenPos) === 1) {
        wanted[tokens[tokenPos]] = true;
        tokenPos += 1;
    }

    var segmentCount = 0;
    var foundReads = {};
    var sourcePos = 0;
    while (Math.sign(sources.length - sourcePos) === 1) {
        var source = sources[sourcePos];
        var matches = [];
        if (source.data.read_name) {
            var readNames = source.data.read_name;
            var sourceKind = "";
            if (source.data.source_kind) sourceKind = source.data.source_kind[0];
            var countSegments = sourceKind === "connector" ? false : true;
            var readPos = 0;
            while (Math.sign(readNames.length - readPos) === 1) {
                var readName = readNames[readPos];
                if (wanted[readName]) {
                    matches.push(readPos);
                    foundReads[readName] = true;
                    if (countSegments) segmentCount += 1;
                }
                readPos += 1;
            }
        }
        source.selected.indices = matches;
        source._lastNonEmpty = matches.slice();
        source.change.emit();
        sourcePos += 1;
    }

    var foundReadCount = 0;
    var foundKeys = Object.keys(foundReads);
    var foundPos = 0;
    while (Math.sign(foundKeys.length - foundPos) === 1) {
        foundReadCount += 1;
        foundPos += 1;
    }

    var missing = [];
    var missingPos = 0;
    while (Math.sign(tokens.length - missingPos) === 1) {
        var requestedName = tokens[missingPos];
        if (!foundReads[requestedName]) missing.push(requestedName);
        missingPos += 1;
    }

    return {
        foundReadCount: foundReadCount,
        segmentCount: segmentCount,
        missing: missing,
    };
};

var orographerReadSearchMessage = function (result) {
    var message = "Matched " + result.foundReadCount + " read";
    if (result.foundReadCount !== 1) message += "s";
    message += " across " + result.segmentCount + " segment";
    if (result.segmentCount !== 1) message += "s";
    if (result.missing.length !== 0) {
        message += ". Missing: " + result.missing.join(", ");
    }
    return message;
};

var modal = document.getElementById("alignmentModal");
var content = document.getElementById("modalContent");
if (!modal) return;
if (!content) return;

content.innerHTML = "";

var title = document.createElement("div");
title.textContent = "Select Reads";
title.style.fontWeight = "bold";
title.style.fontSize = "16px";
title.style.marginBottom = "8px";
content.appendChild(title);

var help = document.createElement("div");
help.textContent = "Enter one or more exact read names, separated by commas or spaces.";
help.style.fontSize = "12px";
help.style.color = "#555";
help.style.marginBottom = "8px";
content.appendChild(help);

var input = document.createElement("textarea");
input.placeholder = "read1 read2,read3";
input.style.width = "100%";
input.style.height = "74px";
input.style.boxSizing = "border-box";
input.style.fontFamily = "Arial, sans-serif";
input.style.fontSize = "13px";
input.style.marginBottom = "8px";
content.appendChild(input);

var buttonRow = document.createElement("div");
buttonRow.style.display = "flex";
buttonRow.style.gap = "8px";
buttonRow.style.marginBottom = "8px";
content.appendChild(buttonRow);

var highlightButton = document.createElement("button");
highlightButton.textContent = "Highlight";
highlightButton.style.padding = "6px 10px";
highlightButton.style.cursor = "pointer";
buttonRow.appendChild(highlightButton);

var status = document.createElement("div");
status.style.fontSize = "12px";
status.style.color = "#555";
status.style.minHeight = "18px";
content.appendChild(status);

var clearSelections = function () {
    orographerClearSelections(all_sources);
    if (typeof window !== "undefined") {
        if (window.orographerUpdateReadConnectionOverlay) {
            window.orographerUpdateReadConnectionOverlay();
        }
    }
};

var highlightReads = function () {
    var tokens = orographerParseReadSearchTokens(input.value);
    if (tokens.length === 0) {
        clearSelections();
        status.textContent = "";
        return;
    }

    var result = orographerApplyReadSearch(all_sources, tokens);
    status.textContent = orographerReadSearchMessage(result);
    modal.style.display = "none";
    if (typeof window !== "undefined") {
        if (window.orographerUpdateReadConnectionOverlay) {
            window.orographerUpdateReadConnectionOverlay();
        }
    }
};

highlightButton.onclick = highlightReads;
input.onkeydown = function (event) {
    if (event.key === "Enter") {
        if (!event.shiftKey) {
            event.preventDefault();
            highlightReads();
        }
    }
};

if (window.fitAlignmentModalToPlot) window.fitAlignmentModalToPlot();
modal.style.display = "flex";
setTimeout(function () { input.focus(); }, 0);
