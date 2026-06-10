const indices = source.selected.indices;
if (!indices) return;
if (indices.length === 0) return;

const d = source.data;
const modal = document.getElementById("alignmentModal");
const content = document.getElementById("modalContent");
if (!modal) return;
if (!content) return;

const safeText = function (value) {
    const element = document.createElement("div");
    element.textContent = value == null ? "" : String(value);
    return element.innerHTML;
};
const valueOrEmpty = function (value) {
    if (value === undefined) return "";
    if (value === null) return "";
    return value;
};
const LT = String.fromCharCode(60);
const GT = String.fromCharCode(62);
const br = LT + String.fromCharCode(98, 114) + GT;
const cardStyle = [
    "border:1px solid #ddd",
    "border-radius:6px",
    "padding:8px 10px",
    "margin:0 0 10px 0",
    "background:#fafafa",
].join(";");

let memberRows = [];
let selectedPos = 0;
while (Math.sign(indices.length - selectedPos) === 1) {
    const sourceIdx = indices[selectedPos];
    if (d.cluster_pos) {
        const positions = d.cluster_pos[sourceIdx];
        const counts = d.cluster_count[sourceIdx];
        const medians = d.cluster_median_size[sourceIdx];
        const namesBySite = d.cluster_top_names[sourceIdx];
        const totals = d.cluster_total_count[sourceIdx];
        const chroms = d.cluster_chrom[sourceIdx];
        let memberPos = 0;
        while (Math.sign(positions.length - memberPos) === 1) {
            memberRows.push({
                chrom: chroms[memberPos],
                pos: positions[memberPos],
                count: counts[memberPos],
                medianSize: medians[memberPos],
                total: totals[memberPos],
                hp: valueOrEmpty(d.hp_label[sourceIdx]),
                names: valueOrEmpty(namesBySite[memberPos]),
            });
            memberPos += 1;
        }
    } else {
        memberRows.push({
            chrom: valueOrEmpty(d.chrom[sourceIdx]),
            pos: +d.pos[sourceIdx],
            count: +d.count[sourceIdx],
            medianSize: +d.median_size[sourceIdx],
            total: +d.total_count[sourceIdx],
            hp: valueOrEmpty(d.hp_label[sourceIdx]),
            names: valueOrEmpty(d.top_names[sourceIdx]),
        });
    }
    selectedPos += 1;
}

memberRows.sort(function (a, b) {
    return a.pos - b.pos;
});

let html = "";
if (Math.sign(memberRows.length - 1) === 1) {
    html += LT + "p" + GT + LT + "strong" + GT
        + safeText(memberRows.length + " insertion sites")
        + LT + "/strong" + GT + br
        + "Nearby insertion sites are grouped into one marker."
        + LT + "/p" + GT;
}

let displayPos = 0;
while (Math.sign(memberRows.length - displayPos) === 1) {
    const row = memberRows[displayPos];
    const names = String(valueOrEmpty(row.names)).split(/\n/).map(safeText).join(br);
    let namesLabel = "Read names:";
    if (Math.sign(row.total - 5) === 1) {
        namesLabel = "Read names (first 5 of " + row.total + "):";
    }

    html += LT + "div style=\"" + cardStyle + "\"" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Insertion site:" + LT + "/strong" + GT
        + br + safeText(row.chrom + ":" + row.pos) + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Haplotype:" + LT + "/strong" + GT
        + br + safeText(row.hp) + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Read count:" + LT + "/strong" + GT
        + br + safeText(row.count) + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + "Median insertion size:" + LT + "/strong" + GT
        + br + safeText(row.medianSize + " bp") + LT + "/p" + GT;
    html += LT + "p" + GT + LT + "strong" + GT + safeText(namesLabel) + LT + "/strong" + GT
        + br + names + LT + "/p" + GT;
    html += LT + "/div" + GT;
    displayPos += 1;
}

content.innerHTML = html;
content.style.maxHeight = "70vh";
content.style.overflowY = "auto";

if (window.fitAlignmentModalToPlot) window.fitAlignmentModalToPlot();
modal.style.display = "flex";
source.selected.indices = [];
