const x = cb_obj.x;
const showCursorGuide = cursor_guide_checkbox.active ? Number.isFinite(x) : false;
cursor_guide_spans.forEach(function (span) {
    span.location = x;
    span.visible = showCursorGuide;
});
