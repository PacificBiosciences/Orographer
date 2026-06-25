const modal = document.getElementById("alignmentModal");
const content = document.getElementById("modalContent");
const title = document.getElementById("alignmentModalTitle");
const dialog = document.getElementById("alignmentModalDialog");
if (!modal) return;
if (!content) return;

content.innerHTML = "";
if (title) title.textContent = title_text;
if (dialog) dialog.style.width = "780px";

content.style.overflowY = "auto";
content.style.maxHeight = "85vh";

function buildNotesList(items) {
    const ul = document.createElement("ul");
    ul.style.fontSize = "13px";
    ul.style.color = "#555555";
    ul.style.margin = "0 0 10px 0";
    ul.style.paddingLeft = "18px";
    for (const item of items) {
        const li = document.createElement("li");
        li.appendChild(item);
        ul.appendChild(li);
    }
    return ul;
}

function addSection(imgUrl, headingText, notesEl) {
    if (headingText) {
        const heading = document.createElement("h3");
        heading.style.fontSize = "14px";
        heading.style.fontWeight = "600";
        heading.style.color = "#1f1f1f";
        heading.style.margin = "16px 0 6px 0";
        heading.appendChild(document.createTextNode(headingText));
        content.appendChild(heading);
    }
    if (notesEl) {
        content.appendChild(notesEl);
    }
    const image = document.createElement("img");
    image.src = imgUrl;
    image.alt = "Reference self-identity dotplot";
    image.style.display = "block";
    image.style.width = "720px";
    image.style.maxWidth = "82vw";
    image.style.height = "auto";
    image.style.background = "#ffffff";
    content.appendChild(image);
}

const identEl = document.createElement("span");
identEl.appendChild(document.createTextNode("Identity is shown in the main figure as "));
const codeEl = document.createElement("code");
codeEl.appendChild(document.createTextNode("Ident"));
identEl.appendChild(codeEl);

const mainNotes = buildNotesList([
    document.createTextNode(region_label),
    identEl,
    document.createTextNode("Full region is binned into 512 bins of size " + Number(bin_size_bp).toLocaleString() + " bp, marked if bins match at 85% identity"),
]);
addSection(image_url, null, mainNotes);

for (const ind of individual_images) {
    const divider = document.createElement("hr");
    divider.style.border = "none";
    divider.style.borderTop = "1px solid #e5e7eb";
    divider.style.margin = "20px 0";
    content.appendChild(divider);
    const indNotes = buildNotesList([
        document.createTextNode(ind.label),
        document.createTextNode("Full region is binned into 512 bins of size " + Number(ind.bin_size_bp).toLocaleString() + " bp, marked if bins match at 85% identity"),
    ]);
    addSection(ind.url, ind.title, indNotes);
}

modal.style.display = "flex";
