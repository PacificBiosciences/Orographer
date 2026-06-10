export default function (args, obj) {
const source = args.source;
const intron_source = args.intron_source;
const x = obj.x;
const y = obj.y;
if (typeof x == "undefined") return;
if (typeof y == "undefined") return;
if (x == null) return;
if (y == null) return;

function isLess(a, b) {
    return Math.sign(a - b) == -1;
}

function isGreater(a, b) {
    return Math.sign(a - b) == 1;
}

function isBetween(value, left, right) {
    if (isLess(value, left)) return false;
    if (isGreater(value, right)) return false;
    return true;
}

function setSelectedIndices(indices) {
    source.selected.indices = indices;
    source.selected.change.emit();
}

const data = source.data;
let hitIndex = -1;

let rowIndex = 0;
while (Math.sign((data.left || []).length - rowIndex) == 1) {
    const left = data.left[rowIndex];
    const right = data.right[rowIndex];
    const top = data.top[rowIndex];
    const bottom = data.bottom[rowIndex];
    if (isBetween(x, left, right)) {
        if (isBetween(y, top, bottom)) {
            hitIndex = rowIndex;
            break;
        }
    }
    rowIndex += 1;
}

if (hitIndex == -1) {
    if (intron_source) {
        const intronData = intron_source.data;
        rowIndex = 0;
        while (Math.sign((intronData.xs || []).length - rowIndex) == 1) {
            const xs = intronData.xs[rowIndex];
            const ys = intronData.ys[rowIndex];
            if (xs) {
                if (ys) {
                    const left = Math.min(xs[0], xs[1]);
                    const right = Math.max(xs[0], xs[1]);
                    const intronY = ys[0];
                    if (isBetween(x, left, right)) {
                        if (!isGreater(Math.abs(y - intronY), 0.28)) {
                            const transcriptId = intronData.transcript_id[rowIndex];
                            hitIndex = (data.transcript_id || []).indexOf(transcriptId);
                            break;
                        }
                    }
                }
            }
            rowIndex += 1;
        }
    }
}

if (hitIndex == -1) {
    let transcriptHeight = 0;
    if (data.transcript_area_height) {
        if (data.transcript_area_height.length) {
            transcriptHeight = data.transcript_area_height[0];
        }
    }
    if (!isGreater(y, transcriptHeight)) {
        setSelectedIndices([]);
        source.change.emit();
    }
    return;
}

const selected = source.selected.indices || [];
const currentIndex = selected.length ? selected[selected.length - 1] : -1;
if (currentIndex == hitIndex) {
    setSelectedIndices([]);
} else {
    setSelectedIndices([hitIndex]);
}
source.change.emit();
}
