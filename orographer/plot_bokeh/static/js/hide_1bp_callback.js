// Toggle visibility of 1bp INDEL renderers when checkbox is toggled.
// When showing again, apply same LOD as variant_lod_callback (threshold 1000).
const hide = cb_obj.active;
one_bp_renderers.forEach(function (r) {
    r.visible = !hide;
});
if (!hide) {
    const visible_range = x_range.end - x_range.start;
    const threshold = 1000;
    const diff = threshold - visible_range;
    const show_text = Math.sign(diff) === 1;
    one_bp_markers.forEach(function (r) {
        r.visible = !show_text;
    });
    one_bp_texts.forEach(function (r) {
        r.visible = show_text;
    });
    one_bp_segments.forEach(function (r) {
        r.visible = true;
    });
}
