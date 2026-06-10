eval(
    "var leftDisabled " + String.fromCharCode(61) + " Math.sign(orig_start - x_range.start) + 1;"
    + "var rightDisabled " + String.fromCharCode(61) + " Math.sign(x_range.end - orig_end) + 1;"
    + "var span " + String.fromCharCode(61) + " x_range.end - x_range.start;"
    + "var zoomOutDisabled " + String.fromCharCode(61) + " leftDisabled ? rightDisabled : false;"
    + "var zoomInDisabled " + String.fromCharCode(61) + " Math.sign(10 - span) + 1;"
    + "left_button.disabled " + String.fromCharCode(61) + " leftDisabled ? true : false;"
    + "right_button.disabled " + String.fromCharCode(61) + " rightDisabled ? true : false;"
    + "zoom_out_button.disabled "
    + String.fromCharCode(61)
    + " zoomOutDisabled ? true : false;"
    + "zoom_in_button.disabled " + String.fromCharCode(61) + " zoomInDisabled ? true : false;"
);
