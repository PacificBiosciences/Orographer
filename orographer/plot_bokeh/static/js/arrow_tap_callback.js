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

// Select all segments for this read (same read_name), not just the clicked one
const clickedReadName = source.data.read_name ? source.data.read_name[bestIdx] : null;
const matchingIndices = [];
if (source.data.read_name) {
    source.data.read_name.forEach(function(readName, j) {
        if (readName === clickedReadName) {
            matchingIndices.push(j);
        }
    });
}
// If no read_name data, fall back to single-segment selection
if (matchingIndices.length === 0) {
    matchingIndices.push(bestIdx);
}

const prev = source._lastNonEmpty.slice();
let next = prev.slice();

// Check if any segment of this read is currently selected
const anySelected = matchingIndices.some(function(idx) { return prev.indexOf(idx) !== -1; });
if (anySelected) {
    next = prev.filter(function(idx) { return matchingIndices.indexOf(idx) === -1; });
} else {
    matchingIndices.forEach(function(idx) {
        if (next.indexOf(idx) === -1) {
            next.push(idx);
        }
    });
}

source.selected.indices = next;
source._lastNonEmpty = next.slice();
source.change.emit();
