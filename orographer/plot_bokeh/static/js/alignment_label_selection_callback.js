const selectedReads = {};
let selectedCount = 0;

selection_sources.forEach(function (source) {
    const indices = source.selected.indices || [];
    const readNames = [];
    if (source.data.read_name) {
        source.data.read_name.forEach(function (readName) {
            readNames.push(readName);
        });
    }
    indices.forEach(function (idx) {
        const readName = readNames[idx];
        if (readName) {
            if (!selectedReads[readName]) {
                selectedReads[readName] = true;
                selectedCount += 1;
            }
        }
    });
});

alignment_label_sources.forEach(function (source) {
    const readNames = [];
    if (source.data.read_name) {
        source.data.read_name.forEach(function (readName) {
            readNames.push(readName);
        });
    }
    const labelAlpha = [];
    readNames.forEach(function (readName, rowIndex) {
        const filterAlpha = source.data.read_filter_alpha
            ? source.data.read_filter_alpha[rowIndex]
            : 1;
        let alpha = filterAlpha ? 0.8 : 0.0;
        if (selectedCount) {
            alpha = 0.0;
            if (selectedReads[readName]) {
                alpha = filterAlpha ? 0.8 : 0.0;
            }
        }
        labelAlpha.push(alpha);
    });
    source.data.label_alpha = labelAlpha;
    source.change.emit();
});
