// Modal close handlers and fitAlignmentModalToPlot
const modal = document.getElementById("alignmentModal");
const closeBtn = document.getElementById("closeModal");
const modalTitle = document.getElementById("alignmentModalTitle");
const modalDialog = document.getElementById("alignmentModalDialog");
window.resetAlignmentModal = function () {
    if (modalTitle) modalTitle.textContent = "Details";
    if (modalDialog) modalDialog.style.width = "400px";
};
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
function consumeModalEvent(event) {
    if (!event) return;
    if (event.preventDefault) event.preventDefault();
    if (event.stopPropagation) event.stopPropagation();
    if (event.stopImmediatePropagation) event.stopImmediatePropagation();
}
function stopModalEvent(event) {
    if (!event) return;
    if (event.stopPropagation) event.stopPropagation();
    if (event.stopImmediatePropagation) event.stopImmediatePropagation();
}
window.closeAlignmentModal = function (event) {
    consumeModalEvent(event);
    modal.style.display = "none";
    window.resetAlignmentModal();
};
if (modalDialog) {
    ["pointerdown", "pointerup", "mousedown", "mouseup", "click"].forEach(function (eventName) {
        modalDialog.addEventListener(eventName, stopModalEvent);
    });
}
if (closeBtn) {
    closeBtn.addEventListener("click", window.closeAlignmentModal);
}
window.addEventListener("click", function (e) {
    if (e.target === modal) window.closeAlignmentModal(e);
});
document.addEventListener("keydown", function (e) {
    if (e.key === "Escape") window.closeAlignmentModal();
});
