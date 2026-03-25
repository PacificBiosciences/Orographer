export default function (args, obj) {
    const start = Math.round(obj.start);
    const end = Math.round(obj.end);
    const chrom = args.chrom;
    const coord_input = args.coord_input;
    const error_div = args.error_div;
    const view_size_div = args.view_size_div;

    function formatNumber(n) {
        const re = new RegExp("\\B(?=(\\d{3})+(?!\\d))", "g");
        return n.toString().replace(re, ",");
    }

    const newValue = chrom + ":" + formatNumber(start) + "-" + formatNumber(end);
    coord_input.value = newValue;
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
    const span = end - start;
    const over1k = Math.sign(span - 1000) + 1;
    const over1m = Math.sign(span - 1000000) + 1;
    const sizeStr = over1k
        ? (over1m ? (span / 1000000).toFixed(1) + " Mb" : (span / 1000).toFixed(1) + " kb")
        : span + " bp";
    if (view_size_div) {
        view_size_div.text = sizeStr;
        view_size_div.change.emit();
    }
}
