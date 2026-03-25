const currentIndices = source.selected.indices;
if (!currentIndices) return;
if (currentIndices.length === 0) return;

const clickedIdx = currentIndices[currentIndices.length - 1];
const data = source.data;

const cdsText = function (col, idx) {
    if (col == null) return "";
    const v = col[idx];
    if (v === undefined) return "";
    if (v === null) return "";
    const isNum = typeof v === "number";
    const nanNum = isNum ? Number.isNaN(v) : false;
    if (nanNum) return "N/A";
    return String(v);
};

const geneName = cdsText(data.gene_name, clickedIdx);
const strand = cdsText(data.gene_strand, clickedIdx);
const exonNumber = cdsText(data.exon_number, clickedIdx);

const modal = document.getElementById("alignmentModal");
const modalContent = document.getElementById("modalContent");
if (!modal) return;
if (!modalContent) return;

const LT = String.fromCharCode(60);
const GT = String.fromCharCode(62);
let html = "";
html += LT + "p" + GT + LT + "strong" + GT + "Gene Name:" + LT + "/strong" + GT
    + LT + "br" + GT + geneName + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Transcript Strand:" + LT + "/strong" + GT
    + LT + "br" + GT + strand + LT + "/p" + GT;
html += LT + "p" + GT + LT + "strong" + GT + "Exon Number:" + LT + "/strong" + GT
    + LT + "br" + GT + exonNumber + LT + "/p" + GT;
modalContent.innerHTML = html;
if (window.fitAlignmentModalToPlot) window.fitAlignmentModalToPlot();
modal.style.display = "flex";

source.selected.indices = [];
