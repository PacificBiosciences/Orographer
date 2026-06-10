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

function orographerMatchingReadIndices(data, readName, fallbackIdx) {
    const matchingIndices = [];
    if (data.read_name) {
        data.read_name.forEach(function (candidateReadName, idx) {
            if (candidateReadName === readName) {
                matchingIndices.push(idx);
            }
        });
    }
    if (matchingIndices.length === 0) {
        matchingIndices.push(fallbackIdx);
    }
    return matchingIndices;
}

function orographerToggleIndices(previousIndices, matchingIndices) {
    const anySelected = matchingIndices.some(function (idx) {
        return previousIndices.indexOf(idx) !== -1;
    });
    if (anySelected) {
        return previousIndices.filter(function (idx) {
            return matchingIndices.indexOf(idx) === -1;
        });
    }

    const next = previousIndices.slice();
    matchingIndices.forEach(function (idx) {
        if (next.indexOf(idx) === -1) {
            next.push(idx);
        }
    });
    return next;
}

if (!source._lastNonEmpty) source._lastNonEmpty = [];

const bestIdx = orographerArrowHitIndex(source.data, cb_obj.x, cb_obj.y);

if (bestIdx === null) {
    return;
}

const clickedReadName = source.data.read_name ? source.data.read_name[bestIdx] : null;
const prev = source._lastNonEmpty.slice();
const matchingIndices = orographerMatchingReadIndices(source.data, clickedReadName, bestIdx);
const next = orographerToggleIndices(prev, matchingIndices);

source.selected.indices = next;
source._lastNonEmpty = next.slice();
source.change.emit();
