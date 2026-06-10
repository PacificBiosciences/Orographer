function orographerArrowHitIndex(data, x, y) {
    if (x === undefined || y === undefined || x === null || y === null) {
        return null;
    }

    const x0s = data.x0;
    const x1s = data.x1;
    const ys = data.y;
    let minX = Infinity;
    let maxX = -Infinity;
    x0s.forEach(function (xVal) {
        if (Math.sign(minX - xVal) === 1) minX = xVal;
        if (Math.sign(xVal - maxX) === 1) maxX = xVal;
    });
    x1s.forEach(function (xVal) {
        if (Math.sign(minX - xVal) === 1) minX = xVal;
        if (Math.sign(xVal - maxX) === 1) maxX = xVal;
    });
    const spanX = Math.max(1e-9, maxX - minX);
    const xTol = spanX / 400.0;
    const yTol = 0.25;

    let bestIdx = null;
    let bestDist = Infinity;

    x0s.forEach(function (segX0, i) {
        const segX1 = x1s[i];
        const segY = ys[i];
        let minSegX, maxSegX;
        if (Math.sign(segX0 - segX1) === 1) {
            minSegX = segX1 - xTol;
            maxSegX = segX0 + xTol;
        } else {
            minSegX = segX0 - xTol;
            maxSegX = segX1 + xTol;
        }

        if (Math.sign(minSegX - x) === 1 || Math.sign(x - maxSegX) === 1) return;

        const dy = Math.abs(y - segY);
        if (Math.sign(dy - yTol) === 1) return;

        const dist = dy;
        if (Math.sign(bestDist - dist) === 1) {
            bestDist = dist;
            bestIdx = i;
        }
    });

    return bestIdx;
}

function orographerSourceSpanX(data) {
    const values = [];
    if (data.x0) values.push.apply(values, data.x0);
    if (data.x1) values.push.apply(values, data.x1);
    if (data.xs) {
        data.xs.forEach(function (row) {
            values.push.apply(values, row);
        });
    }
    if (values.length === 0) {
        return 1;
    }
    let minX = Infinity;
    let maxX = -Infinity;
    values.forEach(function (xVal) {
        if (Math.sign(minX - xVal) === 1) minX = xVal;
        if (Math.sign(xVal - maxX) === 1) maxX = xVal;
    });
    return Math.max(1e-9, maxX - minX);
}

function orographerDistanceToScaledSegment(x, y, x0, y0, x1, y1, xTol, yTol) {
    const scaledX = x / xTol;
    const scaledY = y / yTol;
    const scaledX0 = x0 / xTol;
    const scaledY0 = y0 / yTol;
    const scaledX1 = x1 / xTol;
    const scaledY1 = y1 / yTol;
    const dx = scaledX1 - scaledX0;
    const dy = scaledY1 - scaledY0;
    const len2 = dx * dx + dy * dy;
    let t = 0;
    if (len2 !== 0) {
        t = ((scaledX - scaledX0) * dx + (scaledY - scaledY0) * dy) / len2;
        if (Math.sign(0 - t) === 1) t = 0;
        if (Math.sign(t - 1) === 1) t = 1;
    }
    const projX = scaledX0 + t * dx;
    const projY = scaledY0 + t * dy;
    const distX = scaledX - projX;
    const distY = scaledY - projY;
    return Math.sqrt(distX * distX + distY * distY);
}

function orographerConnectorLineHitIndex(data, x, y) {
    if (!data.xs || !data.ys) {
        return null;
    }
    if (x === undefined || y === undefined || x === null || y === null) {
        return null;
    }
    const spanX = orographerSourceSpanX(data);
    const xTol = spanX / 180.0;
    const yTol = 0.45;
    let bestIdx = null;
    let bestDist = Infinity;
    data.xs.forEach(function (rowXs, rowIndex) {
        const rowYs = data.ys[rowIndex];
        if (!rowYs) return;
        rowXs.forEach(function (segX1, pointIndex) {
            if (pointIndex === 0) return;
            const segY1 = rowYs[pointIndex];
            const segX0 = rowXs[pointIndex - 1];
            const segY0 = rowYs[pointIndex - 1];
            const dist = orographerDistanceToScaledSegment(
                x,
                y,
                segX0,
                segY0,
                segX1,
                segY1,
                xTol,
                yTol,
            );
            if (Math.sign(bestDist - dist) === 1) {
                bestDist = dist;
                bestIdx = rowIndex;
            }
        });
    });
    return Math.sign(1 - bestDist) + 1 ? bestIdx : null;
}

function orographerTappedReadHit(source, allSources, x, y) {
    const arrowIdx = orographerArrowHitIndex(source.data, x, y);
    if (arrowIdx !== null) {
        return { hitSource: source, hitIndex: arrowIdx };
    }
    let bestHit = null;
    allSources.forEach(function (candidateSource) {
        if (bestHit) return;
        const lineIdx = orographerConnectorLineHitIndex(candidateSource.data, x, y);
        if (lineIdx !== null) {
            bestHit = { hitSource: candidateSource, hitIndex: lineIdx };
        }
    });
    return bestHit;
}

function orographerMatchingReadIndices(data, readName) {
    const matchingIndices = [];
    if (data.read_name) {
        data.read_name.forEach(function (candidateReadName, idx) {
            if (candidateReadName === readName) {
                matchingIndices.push(idx);
            }
        });
    }
    return matchingIndices;
}

function orographerSetReadSelectionAcrossSources(allSources, clickedReadName, shouldSelect) {
    allSources.forEach(function (otherSource) {
        const matchingIndices = orographerMatchingReadIndices(
            otherSource.data,
            clickedReadName,
        );
        const currentIndices = otherSource.selected.indices
            ? otherSource.selected.indices.slice()
            : [];
        let newIndices = [];
        if (shouldSelect) {
            newIndices = currentIndices.slice();
            matchingIndices.forEach(function (idx) {
                if (newIndices.indexOf(idx) === -1) {
                    newIndices.push(idx);
                }
            });
        } else {
            newIndices = currentIndices.filter(function (idx) {
                return matchingIndices.indexOf(idx) === -1;
            });
        }
        otherSource.selected.indices = newIndices;
        otherSource._lastNonEmpty = newIndices.slice();
        otherSource.change.emit();
    });
}

if (!source._lastNonEmpty) source._lastNonEmpty = [];

const hit = orographerTappedReadHit(source, all_sources, cb_obj.x, cb_obj.y);

if (hit === null) {
    return;
}

const clickedReadName = hit.hitSource.data.read_name[hit.hitIndex];
const currentSelected = hit.hitSource.selected.indices
    ? hit.hitSource.selected.indices.slice()
    : [];
const isSelected = currentSelected.indexOf(hit.hitIndex) !== -1;

orographerSetReadSelectionAcrossSources(all_sources, clickedReadName, !isSelected);
if (typeof window !== "undefined") {
    if (window.orographerUpdateReadConnectionOverlay) {
        window.orographerUpdateReadConnectionOverlay();
    }
}
