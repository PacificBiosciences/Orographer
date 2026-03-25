const hide = cb_obj.active;
alignment_label_renderers.forEach(function (r) {
    r.visible = !hide;
});
