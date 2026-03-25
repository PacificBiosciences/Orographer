const currentIndices = source.selected.indices;
if (!currentIndices) return;
if (currentIndices.length === 0) return;

const clickedIdx = currentIndices[currentIndices.length - 1];
const data = source.data;
const coords = data.coordinates[clickedIdx] || "";
const variantType = data.variant_type[clickedIdx] || "";
const altAllele = data.alt_allele[clickedIdx] || "";
const haplotypes = data.haplotypes[clickedIdx] || "None";
const sampleLabel = data.sample_label ? data.sample_label[clickedIdx] : null;

const modal = document.getElementById("alignmentModal");
const modalContent = document.getElementById("modalContent");

let canShowModal = false;
if (modal) {
    canShowModal = modalContent ? true : false;
}
if (canShowModal) {
    const LT = String.fromCharCode(60);
    const GT = String.fromCharCode(62);
    let html = "";
    if (sampleLabel) {
        const sampleH = String.fromCharCode(83, 97, 109, 112, 108, 101, 58);
        html += LT + String.fromCharCode(112) + GT + LT
            + String.fromCharCode(115, 116, 114, 111, 110, 103) + GT + sampleH
            + LT + String.fromCharCode(47, 115, 116, 114, 111, 110, 103) + GT
            + LT + String.fromCharCode(98, 114) + GT + sampleLabel
            + LT + String.fromCharCode(47, 112) + GT;
    }
    html += LT + "p" + GT + LT + "strong" + GT + "Coordinates:" + LT + "/strong" + GT
        + LT + "br" + GT + coords + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Variant Type:" + LT + "/strong" + GT
        + LT + "br" + GT + variantType + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Alt Allele:" + LT + "/strong" + GT
        + LT + "br" + GT + altAllele + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Haplotypes:" + LT + "/strong" + GT
        + LT + "br" + GT + haplotypes + LT + "/p" + GT;
    modalContent.innerHTML = html;
    if (window.fitAlignmentModalToPlot) window.fitAlignmentModalToPlot();
    modal.style.display = "flex";
}

source.selected.indices = [];
