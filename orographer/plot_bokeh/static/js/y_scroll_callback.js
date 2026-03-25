// Pan the y-range on mouse wheel (scroll up/down to view alignments above or below).
// cb_obj.delta: positive = wheel down, negative = wheel up (show content above).
const delta = cb_obj.delta;
const span = y_range.end - y_range.start;
const step = -0.2 * span * (delta > 0 ? 1 : -1);
let newStart = y_range.start + step;
let newEnd = y_range.end + step;
// Clamp to full extent so we don't scroll past the data
if (newStart < y_min) {
    newStart = y_min;
    newEnd = Math.min(y_max, y_min + span);
}
if (newEnd > y_max) {
    newEnd = y_max;
    newStart = Math.max(y_min, y_max - span);
}
y_range.start = newStart;
y_range.end = newEnd;
y_range.change.emit();
