const hide = cb_obj.active;
phase_set_marker_renderers.forEach(function (renderer) {
    renderer.visible = !hide;
});
