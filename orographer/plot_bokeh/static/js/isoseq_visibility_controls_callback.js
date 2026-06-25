export default function (args) {
const transcript_sources = args.transcript_sources;
const isoseq_label_sources = args.isoseq_label_sources;
const isoseq_intron_sources = args.isoseq_intron_sources;
const isoseq_intron_arrow_sources = args.isoseq_intron_arrow_sources;
const isoseq_feature_label_sources = args.isoseq_feature_label_sources;
const read_arrow_sources = args.read_arrow_sources;
const read_label_sources = args.read_label_sources;
const plot_figures = args.plot_figures;
const hide_unselected_isoforms_checkbox = args.hide_unselected_isoforms_checkbox;
const hide_unselected_reads_checkbox = args.hide_unselected_reads_checkbox;
const ISOFORM_SELECTED_Y = 0.55;
const ISOFORM_CONTEXT_START_Y = 1.25;
const ISOFORM_CONTEXT_STEP = 0.12;
const ISOFORM_CONTEXT_HEIGHT = 0.08;

function selectedReadNames(source) {
    const names = {};
    let selectedCount = 0;
    if (!source) return { names: names, selectedCount: selectedCount };
    const data = source.data;
    const selected = source.selected.indices || [];
    selected.forEach(function (rowIndex) {
        const readName = data.read_name ? data.read_name[rowIndex] : "";
        if (readName) {
            names[readName] = true;
            selectedCount += 1;
        }
    });
    return { names: names, selectedCount: selectedCount };
}

function selectedReadNamesFromSources(sources) {
    const names = {};
    let selectedCount = 0;
    sources.forEach(function (source) {
        const selectedReads = selectedReadNames(source);
        Object.keys(selectedReads.names).forEach(function (readName) {
            if (!names[readName]) {
                names[readName] = true;
                selectedCount += 1;
            }
        });
    });
    return { names: names, selectedCount: selectedCount };
}

function ensureOriginalY(source) {
    if (!source) return;
    const data = source.data;
    if (!data.y) return;
    const originalY = data.original_y || [];
    if (originalY.length - data.y.length) {
        data.original_y = data.y.map(function (value) {
            return value;
        });
    }
}

function ensureOriginalTranscriptGeometry(source) {
    if (!source) return;
    const data = source.data;
    const originalTop = data.original_top || [];
    if (originalTop.length - data.top.length) {
        data.original_top = data.top.map(function (value) {
            return value;
        });
    }
    const originalBottom = data.original_bottom || [];
    if (originalBottom.length - data.bottom.length) {
        data.original_bottom = data.bottom.map(function (value) {
            return value;
        });
    }
}

function ensureOriginalYs(source) {
    if (!source) return;
    const data = source.data;
    if (!data.ys) return;
    const originalYs = data.original_ys || [];
    if (originalYs.length - data.ys.length) {
        data.original_ys = data.ys.map(function (values) {
            return values.slice();
        });
    }
}

function estimateReadStep(source) {
    ensureOriginalY(source);
    if (!source) return 0.24;
    const data = source.data;
    if (!data.original_y) return 0.24;
    const seen = {};
    const yValues = [];
    data.original_y.forEach(function (yValue) {
        const key = String(yValue);
        if (!seen[key]) {
            seen[key] = true;
            yValues.push(yValue);
        }
    });
    yValues.sort(function (left, right) {
        return left - right;
    });
    let previousY = null;
    let bestStep = null;
    yValues.forEach(function (yValue) {
        if (previousY != null) {
            const diff = yValue - previousY;
            if (diff) {
                if (bestStep == null) {
                    bestStep = diff;
                } else if (Math.sign(bestStep - diff) == 1) {
                    bestStep = diff;
                }
            }
        }
        previousY = yValue;
    });
    return bestStep == null ? 0.24 : bestStep;
}

function packedYByRead(source, selectedReads) {
    ensureOriginalY(source);
    const yByRead = {};
    if (!source) return yByRead;
    const data = source.data;
    if (!data.original_y) return yByRead;
    const readRows = {};
    data.original_y.forEach(function (yValue, rowIndex) {
        const readName = data.read_name ? data.read_name[rowIndex] : "";
        if (!selectedReads.names[readName]) return;
        if (!readRows[readName]) {
            readRows[readName] = { readName: readName, y: yValue };
            return;
        }
        if (Math.sign(readRows[readName].y - yValue) == 1) {
            readRows[readName].y = yValue;
        }
    });
    const rows = Object.keys(readRows).map(function (readName) {
        return readRows[readName];
    });
    rows.sort(function (left, right) {
        const yDiff = left.y - right.y;
        if (yDiff) return yDiff;
        return left.readName.localeCompare(right.readName);
    });
    if (!rows.length) return yByRead;
    let topY = null;
    data.original_y.forEach(function (yValue) {
        if (topY == null) {
            topY = yValue;
        } else if (Math.sign(topY - yValue) == 1) {
            topY = yValue;
        }
    });
    if (topY == null) return yByRead;
    const readStep = estimateReadStep(source);
    rows.forEach(function (row, rank) {
        yByRead[row.readName] = topY + rank * readStep;
    });
    return yByRead;
}

function applyPackedY(source, yByRead, hideReads, hasReadSelection) {
    ensureOriginalY(source);
    if (!source) return;
    const data = source.data;
    if (!data.y) return;
    if (!data.original_y) return;
    data.y = data.original_y.map(function (originalY, rowIndex) {
        if (!hideReads) return originalY;
        if (!hasReadSelection) return originalY;
        const readName = data.read_name ? data.read_name[rowIndex] : "";
        return typeof yByRead[readName] == "undefined" ? originalY : yByRead[readName];
    });
}

function sortedIsoformRows(source, selectedTranscriptId) {
    ensureOriginalTranscriptGeometry(source);
    const data = source.data;
    const rows = {};
    data.transcript_id.forEach(function (transcriptId, rowIndex) {
        const idText = String(transcriptId);
        if (idText == String(selectedTranscriptId)) return;
        const center = (data.original_top[rowIndex] + data.original_bottom[rowIndex]) / 2;
        if (!rows[idText]) {
            rows[idText] = { transcriptId: idText, y: center };
            return;
        }
        if (Math.sign(rows[idText].y - center) == 1) {
            rows[idText].y = center;
        }
    });
    const sortedRows = Object.keys(rows).map(function (transcriptId) {
        return rows[transcriptId];
    });
    sortedRows.sort(function (left, right) {
        const yDiff = left.y - right.y;
        if (yDiff) return yDiff;
        return left.transcriptId.localeCompare(right.transcriptId);
    });
    return sortedRows;
}

function denseIsoformYByTranscript(source, selectedTranscriptId) {
    const data = source.data;
    const yByTranscript = {};
    const rows = sortedIsoformRows(source, selectedTranscriptId);
    const readStart = data.selected_read_y_start ? data.selected_read_y_start[0] : 0;
    const contextEnd = Math.max(ISOFORM_CONTEXT_START_Y, readStart - 0.3);
    const contextSpan = Math.max(0, contextEnd - ISOFORM_CONTEXT_START_Y);
    yByTranscript[String(selectedTranscriptId)] = ISOFORM_SELECTED_Y;
    rows.forEach(function (row, rank) {
        let center = ISOFORM_CONTEXT_START_Y + rank * ISOFORM_CONTEXT_STEP;
        if (rows.length == 1) {
            center = ISOFORM_CONTEXT_START_Y + contextSpan / 2;
        } else if (Math.sign(rows.length - 1) == 1) {
            center = ISOFORM_CONTEXT_START_Y + rank * contextSpan / (rows.length - 1);
        }
        yByTranscript[row.transcriptId] = center;
    });
    return yByTranscript;
}

function applyTranscriptGeometry(source, selectedState, compactIsoforms) {
    ensureOriginalTranscriptGeometry(source);
    const data = source.data;
    if (!compactIsoforms) {
        data.top = data.original_top.map(function (value) {
            return value;
        });
        data.bottom = data.original_bottom.map(function (value) {
            return value;
        });
        return;
    }
    const selectedTranscriptId = selectedState.selectedTranscriptId;
    const yByTranscript = denseIsoformYByTranscript(source, selectedTranscriptId);
    data.top = data.original_top.map(function (originalTop, rowIndex) {
        const transcriptId = String(data.transcript_id[rowIndex]);
        const center = yByTranscript[transcriptId];
        if (transcriptId == String(selectedTranscriptId)) {
            const height = data.original_bottom[rowIndex] - originalTop;
            return center - height / 2;
        }
        return center - ISOFORM_CONTEXT_HEIGHT / 2;
    });
    data.bottom = data.original_bottom.map(function (originalBottom, rowIndex) {
        const transcriptId = String(data.transcript_id[rowIndex]);
        const center = yByTranscript[transcriptId];
        if (transcriptId == String(selectedTranscriptId)) {
            const height = originalBottom - data.original_top[rowIndex];
            return center + height / 2;
        }
        return center + ISOFORM_CONTEXT_HEIGHT / 2;
    });
}

function applySourceYGeometry(source, mainSource, selectedState, compactIsoforms) {
    ensureOriginalY(source);
    if (!source) return;
    const data = source.data;
    if (!data.y) return;
    if (!data.original_y) return;
    if (!compactIsoforms) {
        data.y = data.original_y.map(function (value) {
            return value;
        });
        return;
    }
    const selectedTranscriptId = selectedState.selectedTranscriptId;
    const yByTranscript = denseIsoformYByTranscript(mainSource, selectedTranscriptId);
    data.y = data.original_y.map(function (originalY, rowIndex) {
        const transcriptId = String(data.transcript_id[rowIndex]);
        if (transcriptId == String(selectedTranscriptId)) return ISOFORM_SELECTED_Y;
        return yByTranscript[transcriptId] || originalY;
    });
}

function applyIntronGeometry(source, mainSource, selectedState, compactIsoforms) {
    ensureOriginalYs(source);
    if (!source) return;
    const data = source.data;
    if (!data.ys) return;
    if (!data.original_ys) return;
    if (!compactIsoforms) {
        data.ys = data.original_ys.map(function (values) {
            return values.slice();
        });
        return;
    }
    const selectedTranscriptId = selectedState.selectedTranscriptId;
    const yByTranscript = denseIsoformYByTranscript(mainSource, selectedTranscriptId);
    data.ys = data.original_ys.map(function (values, rowIndex) {
        const transcriptId = String(data.transcript_id[rowIndex]);
        const yValue = yByTranscript[transcriptId] || values[0];
        return [yValue, yValue];
    });
}

function uniqueReadCount(source) {
    if (!source) return 0;
    const names = {};
    const data = source.data;
    (data.read_name || []).forEach(function (readName) {
        if (readName) {
            names[readName] = true;
        }
    });
    return Object.keys(names).length;
}

function selectedReadCount(selectedReads) {
    if (!selectedReads) return 0;
    return Object.keys(selectedReads.names || {}).length;
}

function resizeSelectedView(source, plotFigure, readCount) {
    if (!source) return;
    if (!plotFigure) return;
    const count = Math.max(1, Number(readCount) || 0);
    const readStep = estimateReadStep(source);
    const data = source.data;
    let topY = null;
    const yValues = data.original_y || data.y || [];
    yValues.forEach(function (yValue) {
        if (topY == null) {
            topY = yValue;
        } else if (Math.sign(topY - yValue) == 1) {
            topY = yValue;
        }
    });
    if (topY == null) {
        topY = 1;
    }
    plotFigure.y_range.start = topY + count * readStep + readStep * 3;
    plotFigure.height = Math.min(900, Math.max(260, 160 + count * 9));
}

window.orographerHideUnselectedIsoforms = hide_unselected_isoforms_checkbox.active;
window.orographerHideUnselectedReads = hide_unselected_reads_checkbox.active;

function applyTranscriptVisibility(source) {
    const data = source.data;
    const selected = source.selected.indices || [];
    const hasSelection = Boolean(selected.length);
    const selectedIndex = hasSelection ? selected[selected.length - 1] : -1;
    const selectedTranscriptId = hasSelection ? data.transcript_id[selectedIndex] : null;
    const hideIsoforms = hide_unselected_isoforms_checkbox.active;
    const compactIsoforms = hasSelection ? !hideIsoforms : false;

    applyTranscriptGeometry(source, {
        hasSelection: hasSelection,
        selectedTranscriptId: selectedTranscriptId,
    }, compactIsoforms);
    if (data.alpha) {
        data.alpha = data.alpha.map(function (_oldValue, rowIndex) {
            const baseAlpha = data.base_alpha ? data.base_alpha[rowIndex] : 0.86;
            const isUnassigned = String(data.transcript_id[rowIndex]) == "UNASSIGNED";
            if (!hasSelection) {
                return hideIsoforms ? (isUnassigned ? baseAlpha : 0) : baseAlpha;
            }
            if (!hideIsoforms) {
                return String(data.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? baseAlpha
                    : 0.24;
            }
            const isSelected = String(data.transcript_id[rowIndex]) == String(selectedTranscriptId);
            return (isSelected ? true : isUnassigned) ? baseAlpha : 0;
        });
    }
    if (data.line_alpha) {
        data.line_alpha = data.line_alpha.map(function (_oldValue, rowIndex) {
            const baseAlpha = data.base_line_alpha ? data.base_line_alpha[rowIndex] : 1;
            const isUnassigned = String(data.transcript_id[rowIndex]) == "UNASSIGNED";
            if (!hasSelection) {
                return hideIsoforms ? (isUnassigned ? baseAlpha : 0) : baseAlpha;
            }
            if (!hideIsoforms) {
                return String(data.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? baseAlpha
                    : 0.28;
            }
            const isSelected = String(data.transcript_id[rowIndex]) == String(selectedTranscriptId);
            return (isSelected ? true : isUnassigned) ? baseAlpha : 0;
        });
    }
    source.change.emit();
    return {
        compactIsoforms: compactIsoforms,
        hasSelection: hasSelection,
        selectedTranscriptId: selectedTranscriptId,
    };
}

function applyFeatureVisibility(source, selectedState, alphaColumn, visibleAlpha, contextAlpha) {
    if (!source) return;
    const data = source.data;
    if (!data[alphaColumn]) return;
    const hideIsoforms = hide_unselected_isoforms_checkbox.active;
    data[alphaColumn] = data[alphaColumn].map(function (_oldValue, rowIndex) {
        const isUnassigned = String(data.transcript_id[rowIndex]) == "UNASSIGNED";
        if (!selectedState.hasSelection) {
            return hideIsoforms ? (isUnassigned ? visibleAlpha : 0) : visibleAlpha;
        }
        if (!hideIsoforms) {
            return String(data.transcript_id[rowIndex])
                == String(selectedState.selectedTranscriptId)
                ? visibleAlpha
                : contextAlpha;
        }
        const isSelected = String(data.transcript_id[rowIndex])
            == String(selectedState.selectedTranscriptId);
        return (isSelected ? true : isUnassigned) ? visibleAlpha : 0;
    });
    source.change.emit();
}

function applyReadVisibility(source, visibleAlphaByColumn, selectedReads, yByRead) {
    if (!source) return;
    const data = source.data;
    const hideReads = hide_unselected_reads_checkbox.active;
    const hasReadSelection = Boolean(selectedReads.selectedCount);
    applyPackedY(source, yByRead, hideReads, hasReadSelection);
    Object.keys(visibleAlphaByColumn).forEach(function (column) {
        if (!data[column]) return;
        const visibleAlpha = visibleAlphaByColumn[column];
        data[column] = data[column].map(function (_oldValue, rowIndex) {
            if (!hideReads) return visibleAlpha;
            if (!hasReadSelection) return 0;
            const readName = data.read_name ? data.read_name[rowIndex] : "";
            return selectedReads.names[readName] ? visibleAlpha : 0;
        });
    });
    source.change.emit();
}

transcript_sources.forEach(function (source, sourceIndex) {
    const selectedState = applyTranscriptVisibility(source);
    applySourceYGeometry(
        isoseq_label_sources[sourceIndex],
        source,
        selectedState,
        selectedState.compactIsoforms
    );
    applyIntronGeometry(
        isoseq_intron_sources[sourceIndex],
        source,
        selectedState,
        selectedState.compactIsoforms
    );
    applySourceYGeometry(
        isoseq_intron_arrow_sources[sourceIndex],
        source,
        selectedState,
        selectedState.compactIsoforms
    );
    applySourceYGeometry(
        isoseq_feature_label_sources[sourceIndex],
        source,
        selectedState,
        selectedState.compactIsoforms
    );
    applyFeatureVisibility(isoseq_label_sources[sourceIndex], selectedState, "alpha", 1, 0.18);
    applyFeatureVisibility(
        isoseq_intron_sources[sourceIndex],
        selectedState,
        "line_alpha",
        0.75,
        0.22
    );
    applyFeatureVisibility(
        isoseq_intron_arrow_sources[sourceIndex],
        selectedState,
        "alpha",
        0.85,
        0.22
    );
    applyFeatureVisibility(
        isoseq_feature_label_sources[sourceIndex],
        selectedState,
        "label_alpha",
        0.8,
        0
    );
});

const selectedReads = selectedReadNamesFromSources(
    [].concat(read_arrow_sources, read_label_sources)
);
read_arrow_sources.forEach(function (source, sourceIndex) {
    const yByRead = packedYByRead(source, selectedReads);
    applyReadVisibility(source, {
        arrow_line_alpha: 0.5,
        arrowhead_alpha: 0.65,
        arrow_nonselected_alpha: 0.12,
    }, selectedReads, yByRead);
    applyReadVisibility(
        read_label_sources[sourceIndex],
        { label_alpha: 0.8 },
        selectedReads,
        yByRead
    );
    if (hide_unselected_reads_checkbox.active) {
        if (selectedReads.selectedCount) {
            resizeSelectedView(
                source,
                plot_figures[sourceIndex],
                selectedReadCount(selectedReads)
            );
        } else {
            resizeSelectedView(source, plot_figures[sourceIndex], 0);
        }
    } else {
        resizeSelectedView(source, plot_figures[sourceIndex], uniqueReadCount(source));
    }
});
}
