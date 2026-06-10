const modal = document.getElementById("alignmentModal");
const content = document.getElementById("modalContent");
const title = document.getElementById("alignmentModalTitle");
const dialog = document.getElementById("alignmentModalDialog");
if (!modal) return;
if (!content) return;

content.innerHTML = "";
if (title) title.textContent = title_text;
if (dialog) dialog.style.width = "760px";

const label = document.createElement("div");
label.style.fontSize = "13px";
label.style.color = "#555555";
label.style.marginBottom = "10px";
label.appendChild(document.createTextNode(region_label));

const image = document.createElement("img");
image.src = image_url;
image.alt = "Reference self-identity dotplot";
image.style.display = "block";
image.style.width = "720px";
image.style.maxWidth = "82vw";
image.style.height = "auto";
image.style.background = "#ffffff";

content.appendChild(label);
content.appendChild(image);
modal.style.display = "flex";
