const text = coord_input.value;
const cleaned = text.replace(new RegExp("[\\s,]", "g"), "");
const re = new RegExp(":?(\\d+)(?:-(\\d+))?$");
const match = cleaned.match(re);
let isValid = false;

if (match) {
    let start = parseInt(match[1], 10);
    let end = match[2] ? parseInt(match[2], 10) : start + 1;

    if (!isNaN(start)) {
        if (!isNaN(end)) {
            if (end === start) {
                end = start + 1;
            }
            if (Math.sign(end - start) === 1) {
                let overlapsAtAll = false;
                if (Math.sign(orig_end - start) === 1) {
                    if (Math.sign(end - orig_start) === 1) {
                        overlapsAtAll = true;
                    }
                }

                if (overlapsAtAll) {
                    const regionSize = end - start;
                    if (Math.sign(10 - regionSize) === 1) {
                        const center = (start + end) / 2;
                        start = Math.floor(center - 5);
                        end = Math.ceil(center + 5);

                        // Shift expanded region to fit within bounds when it overflows.
                        if (Math.sign(start - orig_start) === -1) {
                            start = orig_start;
                            end = orig_start + 10;
                            if (Math.sign(end - orig_end) === 1) {
                                end = orig_end;
                                start = orig_end - 10;
                                if (Math.sign(start - orig_start) === -1) {
                                    start = orig_start;
                                }
                            }
                        } else if (Math.sign(end - orig_end) === 1) {
                            end = orig_end;
                            start = orig_end - 10;
                            if (Math.sign(start - orig_start) === -1) {
                                start = orig_start;
                            }
                        }
                    }
                }

                const startInBounds = Math.sign(start - orig_start) !== -1;
                const endInBounds = Math.sign(orig_end - end) !== -1;
                const bothInBounds = startInBounds ? endInBounds : false;

                if (bothInBounds) {
                    x_range.start = start;
                    x_range.end = end;
                    isValid = true;
                    error_div.text = "";
                    error_div.styles = {
                        "font-size": "12px",
                        "font-family": "Arial, sans-serif",
                        "color": "#ff0000",
                        "text-align": "left",
                        "padding-left": "8px",
                        "min-height": "0",
                        "height": "0",
                        "overflow": "hidden",
                        "margin": "0",
                        "padding-top": "0",
                        "padding-bottom": "0",
                    };
                }
            }
        }
    }
}

if (!isValid) {
    error_div.text = "Coordinate not found in target region";
    error_div.styles = {
        "font-size": "12px",
        "font-family": "Arial, sans-serif",
        "color": "#ff0000",
        "text-align": "left",
        "padding-left": "8px",
        "margin": "0",
        "padding-top": "2px",
        "padding-bottom": "0",
        "overflow": "visible",
        "min-height": "0",
        "height": "auto",
    };
}
