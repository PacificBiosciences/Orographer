export default function (args) {
const source = args.source;
const indices = source.selected.indices || [];
if (!indices.length) {
    return;
}

const idx = indices[indices.length - 1];
const data = source.data;
const modal = document.getElementById("alignmentModal");
const modalContent = document.getElementById("modalContent");
if (!modal) {
    return;
}
if (!modalContent) {
    return;
}

const LT = String.fromCharCode(60);
const GT = String.fromCharCode(62);

function value(column, fallback) {
    if (!data[column]) return fallback;
    const raw = data[column][idx];
    if (raw == null) return fallback;
    if (String(raw) == "") return fallback;
    return raw;
}

function countValue(column, fallback) {
    const raw = value(column, fallback);
    const numeric = Number(raw);
    if (!Number.isFinite(numeric)) return raw;
    return numeric.toLocaleString("en-US");
}

let html = "";
html += LT + "p" + GT + LT + "strong" + GT + "Feature:" + LT + "/strong" + GT
    + LT + "br" + GT + value("feature_label", "Unknown") + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Feature Coordinates:" + LT + "/strong" + GT
    + LT + "br" + GT + value("feature_coordinates", "Unknown") + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Gene Name:" + LT + "/strong" + GT
    + LT + "br" + GT + value("gene_name", "Unknown") + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Transcript Name:" + LT + "/strong" + GT
    + LT + "br" + GT + value("transcript_name", value("transcript_id", "Unknown"))
    + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Transcript Coordinates:" + LT
    + "/strong" + GT + LT + "br" + GT + value("coordinates", "Unknown")
    + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Transcript Span:" + LT
    + "/strong" + GT + LT + "br" + GT + value("span", "Unknown") + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Exons:" + LT + "/strong" + GT
    + LT + "br" + GT + value("exon_count", "0") + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Introns:" + LT + "/strong" + GT
    + LT + "br" + GT + value("intron_count", "0") + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Strand Orientation:" + LT
    + "/strong" + GT + LT + "br" + GT + value("strand", "unknown") + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Reads Assigned to Transcript:" + LT
    + "/strong" + GT + LT + "br" + GT + countValue("assigned_read_count", "0")
    + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Transcripts in Gene:" + LT
    + "/strong" + GT + LT + "br" + GT + value("gene_transcript_count", "0")
    + LT + "/p" + GT;

modalContent.innerHTML = html;
if (window.fitAlignmentModalToPlot) window.fitAlignmentModalToPlot();
modal.style.display = "flex";
source.selected.indices = [];
}
