export default function (args) {
const gene_select = args.gene_select;
const transcript_sources = args.transcript_sources;
const read_sources = args.read_sources;
const selectedGene = gene_select.value || "ALL";

function rowVisible(data, rowIndex) {
    if (selectedGene == "ALL") {
        return true;
    }
    if (!data.gene_id) {
        return false;
    }
    return String(data.gene_id[rowIndex]) == String(selectedGene);
}

transcript_sources.forEach(function (source) {
    const data = source.data;
    if (data.alpha) {
        data.alpha = data.alpha.map(function (_oldValue, rowIndex) {
            const baseAlpha = data.base_alpha ? data.base_alpha[rowIndex] : 0.86;
            return rowVisible(data, rowIndex) ? baseAlpha : 0;
        });
    }
    if (data.line_alpha) {
        data.line_alpha = data.line_alpha.map(function (_oldValue, rowIndex) {
            const baseAlpha = data.base_line_alpha ? data.base_line_alpha[rowIndex] : 0.75;
            return rowVisible(data, rowIndex) ? baseAlpha : 0;
        });
    }
    source.change.emit();
});

read_sources.forEach(function (source) {
    const data = source.data;
    if (data.arrow_line_alpha) {
        data.arrow_line_alpha = data.arrow_line_alpha.map(function (_oldValue, rowIndex) {
            return rowVisible(data, rowIndex) ? 0.5 : 0;
        });
    }
    if (data.arrowhead_alpha) {
        data.arrowhead_alpha = data.arrowhead_alpha.map(function (_oldValue, rowIndex) {
            return rowVisible(data, rowIndex) ? 0.65 : 0;
        });
    }
    if (data.arrow_selected_alpha) {
        data.arrow_selected_alpha = data.arrow_selected_alpha.map(function (_oldValue, rowIndex) {
            return rowVisible(data, rowIndex) ? 1 : 0;
        });
    }
    if (data.arrow_nonselected_alpha) {
        data.arrow_nonselected_alpha = data.arrow_nonselected_alpha.map(
            function (_oldValue, rowIndex) {
                return rowVisible(data, rowIndex)
                    ? 0.12
                    : 0;
            }
        );
    }
    if (data.label_alpha) {
        data.label_alpha = data.label_alpha.map(function (_oldValue, rowIndex) {
            return rowVisible(data, rowIndex) ? 0.8 : 0;
        });
    }
    if (data.alpha) {
        data.alpha = data.alpha.map(function (_oldValue, rowIndex) {
            return rowVisible(data, rowIndex) ? 0.85 : 0;
        });
    }
    if (data.line_alpha) {
        data.line_alpha = data.line_alpha.map(function (_oldValue, rowIndex) {
            return rowVisible(data, rowIndex) ? 0.75 : 0;
        });
    }
    source.change.emit();
});
}
