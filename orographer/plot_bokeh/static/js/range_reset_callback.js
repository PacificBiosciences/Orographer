x_range.start = x_start;
x_range.end = x_end;
const n = y_ranges.length;
for (let i = 0; i !== n; i++) {
    y_ranges[i].start = y_bounds[i][0];
    y_ranges[i].end = y_bounds[i][1];
}
all_sources.forEach(function (source) {
    source.selected.indices = [];
    if (source._lastNonEmpty) {
        source._lastNonEmpty = [];
    }
    source.change.emit();
});
