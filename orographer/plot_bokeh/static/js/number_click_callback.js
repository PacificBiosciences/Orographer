const currentIndices = source.selected.indices;
if (!currentIndices) {
    return;
}
if (!currentIndices.length) {
    return;
}

const clickedIdx = currentIndices[currentIndices.length - 1];
const data = source.data;
const readName = data.read_name[clickedIdx] || "";
const alignNum = data.alignment_number[clickedIdx] || "";
const strand = data.strand[clickedIdx] || "";
const coords = data.coordinates[clickedIdx] || "";
const hap = data.haplotype[clickedIdx] || "";
let sampleLabel = null;
if (data.sample_label) {
    sampleLabel = data.sample_label[clickedIdx];
}
let inclusionReason = "";
if (data.inclusion_reason) {
    inclusionReason = data.inclusion_reason[clickedIdx];
}

const modal = document.getElementById("alignmentModal");
const modalContent = document.getElementById("modalContent");

let canShowModal = false;
if (modal) {
    if (modalContent) {
        canShowModal = true;
    }
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
    html += LT + "p" + GT + LT + "strong" + GT + "Read Name:" + LT + "/strong" + GT
        + LT + "br" + GT + readName + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Haplotype:" + LT + "/strong" + GT
        + LT + "br" + GT + hap + LT + "/p" + GT;
    if (inclusionReason) {
        html += LT + "p" + GT + LT + "strong" + GT
            + "Primary inclusion reason:" + LT + "/strong" + GT
            + LT + "br" + GT + inclusionReason + LT + "/p" + GT;
    }
    html += LT + "p" + GT + LT + "strong" + GT + "Alignment #:" + LT + "/strong" + GT
        + LT + "br" + GT + alignNum + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Strand:" + LT + "/strong" + GT
        + LT + "br" + GT + strand + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Coordinates:" + LT + "/strong" + GT
        + LT + "br" + GT + coords + LT + "/p" + GT;

    let alignmentNumbers = [];
    let alignmentCoordinates = [];
    if (data.all_alignment_numbers) {
        alignmentNumbers = data.all_alignment_numbers[clickedIdx] || [];
    }
    if (data.all_alignment_coordinates) {
        alignmentCoordinates = data.all_alignment_coordinates[clickedIdx] || [];
    }
    if (!alignmentNumbers.length) {
        all_sources.forEach(function (labelSource) {
            const labelData = labelSource.data;
            if (!labelData.read_name) {
                return;
            }
            labelData.read_name.forEach(function (candidateReadName, rowIndex) {
                if (String(candidateReadName).localeCompare(String(readName))) {
                    return;
                }
                alignmentNumbers.push(labelData.alignment_number[rowIndex] || "");
                alignmentCoordinates.push(labelData.coordinates[rowIndex] || "");
            });
        });
    }

    let explicitChimeric = true;
    if (String(inclusionReason).localeCompare("Chimeric")) {
        explicitChimeric = false;
    }
    const hasMultipleRows = Math.sign(alignmentNumbers.length - 2) + 1;
    let isChimericRead = hasMultipleRows;
    if (explicitChimeric) {
        isChimericRead = true;
    }
    let showAlignmentRows = false;
    if (isChimericRead) {
        if (alignmentNumbers.length) {
            showAlignmentRows = true;
        }
    }
    if (showAlignmentRows) {
        html += LT + "p" + GT + LT + "strong" + GT
            + "All alignment coordinates:" + LT + "/strong" + GT + LT + "/p" + GT;
        html += LT + "ol" + GT;
        alignmentNumbers.forEach(function (alignmentNumber, rowIndex) {
            html += LT + "li" + GT + "Alignment " + alignmentNumber + ": "
                + alignmentCoordinates[rowIndex] + LT + "/li" + GT;
        });
        html += LT + "/ol" + GT;
    }
    modalContent.innerHTML = html;
    if (window.fitAlignmentModalToPlot) window.fitAlignmentModalToPlot();
    modal.style.display = "flex";
}
source.selected.indices = [];
