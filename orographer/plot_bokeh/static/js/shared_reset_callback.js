all_sources.forEach(function(source) {
    source.selected.indices = [];
    if (source._lastNonEmpty) {
        source._lastNonEmpty = [];
    }
    source.change.emit();
});
