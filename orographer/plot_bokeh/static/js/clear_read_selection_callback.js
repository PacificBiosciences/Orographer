all_sources.forEach(function (source) {
    source.selected.indices = [];
    source._lastNonEmpty = [];
    source.change.emit();
});
if (typeof clear_button !== "undefined") {
    clear_button.disabled = true;
}
if (typeof window !== "undefined") {
    if (window.orographerUpdateReadConnectionOverlay) {
        window.orographerUpdateReadConnectionOverlay();
    }
}
