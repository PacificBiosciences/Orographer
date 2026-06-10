var selectableSources = [];
try {
    selectableSources = all_sources || [];
} catch (error) {
    selectableSources = [];
}
if (window) {
    window.orographerSelectableSources = selectableSources;
    if (window.orographerUpdateReadConnectionOverlay) {
        window.orographerUpdateReadConnectionOverlay();
    }
}
