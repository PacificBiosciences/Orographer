export default function (args) {
    var x_range = args.x_range;
    var marker_renderers = args.marker_renderers;
    var text_renderers = args.text_renderers;
    var one_bp_renderers = args.one_bp_renderers;
    var hide_1bp_checkbox = args.hide_1bp_checkbox;

    var visible_range = x_range.end - x_range.start;
    var threshold = 1000;
    var diff = threshold - visible_range;
    var show_text = Math.sign(diff) === 1;

    marker_renderers.forEach(function (renderer) {
        renderer.visible = !show_text;
    });

    text_renderers.forEach(function (renderer) {
        renderer.visible = show_text;
    });

    var one_bp_defined = typeof one_bp_renderers !== "undefined";
    var hide_active = hide_1bp_checkbox ? hide_1bp_checkbox.active : false;
    var do_hide = one_bp_defined ? one_bp_renderers ? hide_active : false : false;
    if (do_hide) {
        one_bp_renderers.forEach(function (renderer) {
            renderer.visible = false;
        });
    }
}
