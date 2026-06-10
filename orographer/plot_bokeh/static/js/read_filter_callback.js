function readFilterSources() {
    const sources = read_filter_sources || [];
    const uniqueSources = [];
    const seen = new Set();
    sources.forEach(function (source) {
        if (!source) {
            return;
        }
        if (!source.data) {
            return;
        }
        if (seen.has(source)) {
            return;
        }
        seen.add(source);
        uniqueSources.push(source);
    });
    return uniqueSources;
}

function readFilterLayoutMode() {
    if (hide_non_split_checkbox.active) {
        return hide_non_multiregion_checkbox.active ? "split_multiregion" : "split";
    }
    return hide_non_multiregion_checkbox.active ? "multiregion" : "all";
}

function readFilterVisibility(data, rowIndex, layoutMode) {
    const columnName = "read_filter_visible_" + layoutMode;
    if (!data[columnName]) {
        return 1;
    }
    return data[columnName][rowIndex] ? 1 : 0;
}

function updateAlphaColumn(data, columnName, visibleValue, layoutMode) {
    if (!data[columnName]) {
        return;
    }
    const values = [];
    data[columnName].forEach(function (_oldValue, rowIndex) {
        const visible = readFilterVisibility(data, rowIndex, layoutMode);
        values.push(visible * visibleValue);
    });
    data[columnName] = values;
}

function shiftedLineRows(rows, oldValues, newValues) {
    const shifted = [];
    rows.forEach(function (row, rowIndex) {
        const delta = newValues[rowIndex] - oldValues[rowIndex];
        if (Array.isArray(row)) {
            const shiftedRow = [];
            row.forEach(function (value) {
                shiftedRow.push(value + delta);
            });
            shifted.push(shiftedRow);
        } else {
            shifted.push(row);
        }
    });
    return shifted;
}

function copyLineRows(rows) {
    const copied = [];
    rows.forEach(function (row) {
        copied.push(Array.isArray(row) ? row.slice() : row);
    });
    return copied;
}

function updateScalarLayoutColumn(data, columnName, layoutMode) {
    const modeColumn = columnName + "_" + layoutMode;
    if (!data[columnName]) {
        return;
    }
    if (!data[modeColumn]) {
        return;
    }
    const oldValues = data[columnName];
    const newValues = data[modeColumn].slice();
    if (columnName == "y") {
        if (data.chevron_ys) {
            data.chevron_ys = shiftedLineRows(data.chevron_ys, oldValues, newValues);
        }
    }
    data[columnName] = newValues;
}

function updateLineLayoutColumn(data, columnName, layoutMode) {
    const modeColumn = columnName + "_" + layoutMode;
    if (!data[columnName]) {
        return;
    }
    if (!data[modeColumn]) {
        return;
    }
    data[columnName] = copyLineRows(data[modeColumn]);
}

function updateLayoutColumns(data, layoutMode) {
    updateScalarLayoutColumn(data, "y", layoutMode);
    updateScalarLayoutColumn(data, "y0", layoutMode);
    updateScalarLayoutColumn(data, "y1", layoutMode);
    updateLineLayoutColumn(data, "ys", layoutMode);
}

function selectedReadNames(sources) {
    const reads = {};
    let count = 0;
    sources.forEach(function (source) {
        const data = source.data;
        if (!data.read_name) {
            return;
        }
        if (!data.source_kind) {
            return;
        }
        const indices = source.selected ? source.selected.indices || [] : [];
        indices.forEach(function (rowIndex) {
            const readName = data.read_name[rowIndex];
            if (readName) {
                if (!reads[readName]) {
                    reads[readName] = true;
                    count += 1;
                }
            }
        });
    });
    return { reads: reads, count: count };
}

function updateReadFilterYRanges(layoutMode) {
    if (!read_filter_y_ranges) {
        return;
    }
    read_filter_y_ranges.forEach(function (yRange, rangeIndex) {
        const boundsByMode = read_filter_y_bounds ? read_filter_y_bounds[rangeIndex] : null;
        if (!boundsByMode) {
            return;
        }
        const bounds = boundsByMode[layoutMode];
        if (!bounds) {
            return;
        }
        yRange.start = bounds[0];
        yRange.end = bounds[1];
    });
}

function applyReadFilter() {
    const layoutMode = readFilterLayoutMode();
    const sources = readFilterSources();
    const selected = selectedReadNames(sources);
    const changedSources = [];

    sources.forEach(function (source) {
        const data = source.data;
        updateLayoutColumns(data, layoutMode);
        updateAlphaColumn(data, "read_filter_alpha", 1, layoutMode);
        updateAlphaColumn(data, "arrow_line_alpha", 0.5, layoutMode);
        updateAlphaColumn(data, "arrowhead_alpha", 0.65, layoutMode);
        updateAlphaColumn(data, "arrow_selected_alpha", 1, layoutMode);
        updateAlphaColumn(data, "arrow_nonselected_alpha", 0.12, layoutMode);
        updateAlphaColumn(data, "connector_line_alpha", data.xs ? 0.68 : 0.9, layoutMode);
        updateAlphaColumn(data, "connector_selected_alpha", data.xs ? 0.82 : 1, layoutMode);
        updateAlphaColumn(data, "connector_nonselected_alpha", 0.18, layoutMode);
        updateAlphaColumn(data, "marker_fill_alpha", 0.8, layoutMode);
        updateAlphaColumn(data, "marker_line_alpha", 0.9, layoutMode);
        updateAlphaColumn(data, "line_alpha", 1, layoutMode);
        updateAlphaColumn(data, "separator_line_alpha", 0.3, layoutMode);
        updateAlphaColumn(data, "text_alpha", 1, layoutMode);

        if (data.label_alpha) {
            const labelAlpha = [];
            data.label_alpha.forEach(function (_oldValue, rowIndex) {
                const readName = data.read_name ? data.read_name[rowIndex] : "";
                let alpha = readFilterVisibility(data, rowIndex, layoutMode) ? 0.8 : 0;
                if (selected.count) {
                    alpha = 0;
                    if (selected.reads[readName]) {
                        alpha = readFilterVisibility(data, rowIndex, layoutMode) ? 0.8 : 0;
                    }
                }
                labelAlpha.push(alpha);
            });
            data.label_alpha = labelAlpha;
        }
        changedSources.push(source);
    });

    updateReadFilterYRanges(layoutMode);
    changedSources.forEach(function (source) {
        source.change.emit();
    });

    if (window.orographerUpdateReadConnectionOverlay) {
        window.orographerUpdateReadConnectionOverlay();
    }
}

window.orographerApplyReadFilter = applyReadFilter;
applyReadFilter();
