// Guard against programmatic updates during slot-swap; the swap callback sets this flag.
if (window.orographerSwapInProgress) { return; }

// Read navigation bounds from updatable Div models so slot-swap can update them.
var _s = parseInt(orig_start_div.text, 10);
var _e = parseInt(orig_end_div.text, 10);

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
                if (Math.sign(_e - start) === 1) {
                    if (Math.sign(end - _s) === 1) {
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
                        if (Math.sign(start - _s) === -1) {
                            start = _s;
                            end = _s + 10;
                            if (Math.sign(end - _e) === 1) {
                                end = _e;
                                start = _e - 10;
                                if (Math.sign(start - _s) === -1) {
                                    start = _s;
                                }
                            }
                        } else if (Math.sign(end - _e) === 1) {
                            end = _e;
                            start = _e - 10;
                            if (Math.sign(start - _s) === -1) {
                                start = _s;
                            }
                        }
                    }
                }

                const startInBounds = Math.sign(start - _s) !== -1;
                const endInBounds = Math.sign(_e - end) !== -1;
                const bothInBounds = startInBounds ? endInBounds : false;

                if (bothInBounds) {
                    x_range.start = start;
                    x_range.end = end;
                    x_range.change.emit();
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
