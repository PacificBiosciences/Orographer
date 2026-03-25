// Modal close handlers and fitAlignmentModalToPlot
const modal = document.getElementById("alignmentModal");
const closeBtn = document.getElementById("closeModal");
window.fitAlignmentModalToPlot = function () {
    const wrapper = document.getElementById("alignmentModalWrapper");
    const dialog = document.getElementById("alignmentModalDialog");
    const root = document.querySelector(".bk-root");
    if (!wrapper || !dialog) return;
    const maxH = root
        ? Math.max(120, 0.9 * root.getBoundingClientRect().height)
        : 0.9 * window.innerHeight;
    dialog.style.maxHeight = "";
    dialog.style.transform = "";
    const naturalH = dialog.scrollHeight;
    const naturalW = dialog.offsetWidth || 400;
    if (naturalH > maxH) {
        if (maxH > 0) {
            const scale = maxH / naturalH;
            dialog.style.transform = "scale(" + scale + ")";
            dialog.style.transformOrigin = "top left";
            wrapper.style.width = naturalW * scale + "px";
            wrapper.style.height = naturalH * scale + "px";
            wrapper.style.overflow = "hidden";
        }
    } else {
        dialog.style.transformOrigin = "";
        wrapper.style.width = "";
        wrapper.style.height = "";
        wrapper.style.overflow = "";
    }
};
if (closeBtn) closeBtn.onclick = function () { modal.style.display = "none"; };
window.onclick = function (e) {
    if (e.target === modal) modal.style.display = "none";
};
document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") modal.style.display = "none";
});
