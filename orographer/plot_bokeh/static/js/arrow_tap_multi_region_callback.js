// Initialize storage
if (!source._lastNonEmpty) source._lastNonEmpty = [];

// Data arrays
const x0s = source.data.x0;
const x1s = source.data.x1;
const ys  = source.data.y;

// Tap event coordinates in data space
const x = cb_obj.x;
const y = cb_obj.y;

// Validate coordinates
if (x === undefined || y === undefined || x === null || y === null) {
    return;
}

// Compute overall span to set a loose tolerance
let minX = Infinity;
let maxX = -Infinity;
x0s.forEach(function(xVal) {
    if (Math.sign(minX - xVal) === 1) minX = xVal;
    if (Math.sign(xVal - maxX) === 1) maxX = xVal;
});
x1s.forEach(function(xVal) {
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

// If no arrow hit, do nothing (preserve current selection)
if (bestIdx === null) {
    return;
}

// Get the read name of the clicked arrow
const clickedReadName = source.data.read_name[bestIdx];

// Toggle logic: check if this arrow is already selected
const currentSelected = source.selected.indices ? source.selected.indices.slice() : [];
const isSelected = currentSelected.indexOf(bestIdx) !== -1;

// If clicking on a selected arrow, deselect all matching reads across all regions
// Otherwise, select all matching reads across all regions
if (isSelected) {
    // Deselect: clear all matching reads
    all_sources.forEach(function (otherSource) {
        const matchingIndices = [];
        if (otherSource.data.read_name) {
            otherSource.data.read_name.forEach(function (readName, j) {
                if (readName === clickedReadName) {
                    matchingIndices.push(j);
                }
            });
        }
        const currentIndices = otherSource.selected.indices
            ? otherSource.selected.indices.slice()
            : [];
        const newIndices = currentIndices.filter(function (idx) {
            return matchingIndices.indexOf(idx) === -1;
        });
        otherSource.selected.indices = newIndices;
        otherSource._lastNonEmpty = newIndices.slice();
        otherSource.change.emit();
    });
} else {
    // Select: highlight all matching reads across all regions
    all_sources.forEach(function (otherSource) {
        const matchingIndices = [];
        if (otherSource.data.read_name) {
            otherSource.data.read_name.forEach(function (readName, j) {
                if (readName === clickedReadName) {
                    matchingIndices.push(j);
                }
            });
        }
        const currentIndices = otherSource.selected.indices
            ? otherSource.selected.indices.slice()
            : [];
        const combined = currentIndices.slice();
        matchingIndices.forEach(function (idx) {
            if (combined.indexOf(idx) === -1) {
                combined.push(idx);
            }
        });
        otherSource.selected.indices = combined;
        otherSource._lastNonEmpty = combined.slice();
        otherSource.change.emit();
    });
}
