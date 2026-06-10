let selectedCount = 0;
all_sources.forEach(function (source) {
    const indices = source.selected.indices || [];
    selectedCount += indices.length;
});
eval(
    "clear_button.disabled "
    + String.fromCharCode(61)
    + selectedCount
    + " ? false : true;"
);
