function geneTrackMode() {
    const span = Math.abs(x_range.end - x_range.start) + 1;
    const detailActive = Math.sign(detail_span_bp - span) + 1;
    if (detailActive) {
        return "detail";
    }
    const mediumActive = Math.sign(medium_span_bp - span) + 1;
    return mediumActive ? "medium" : "overview";
}

function applyGeneTrackMode() {
    const mode = geneTrackMode();
    sources.forEach(function (source, index) {
        const data = source.data;
        const targetColumn = target_columns[index];
        const modeColumn = targetColumn + "_" + mode;
        if (!data[targetColumn]) {
            return;
        }
        if (!data[modeColumn]) {
            return;
        }
        data[targetColumn] = data[modeColumn].slice();
        source.change.emit();
    });
}

applyGeneTrackMode();
