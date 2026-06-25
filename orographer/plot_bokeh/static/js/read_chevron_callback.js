const data = source.data;
if (!data.x0 || !data.x1 || !data.y) return;

function chevronLength() {
    const width = plot_figure.inner_width || plot_figure.width || 1;
    const rangeSpan = Math.abs(x_range.end - x_range.start);
    return Math.max(1, rangeSpan / Math.max(1, width) * target_px);
}

function chevronYOffset() {
    const height = plot_figure.inner_height || plot_figure.height || 1;
    const rangeSpan = Math.abs(y_range.end - y_range.start);
    return Math.max(0.001, rangeSpan / Math.max(1, height) * target_y_px);
}

function chevronVertices(x0Value, x1Value, yValue, targetLength, targetYOffset, rowIndex) {
    const span = Math.abs(x1Value - x0Value);
    if (span == 0) return { xs: [], ys: [] };
    const direction = Math.sign(x1Value - x0Value) == 1 ? 1 : -1;
    const length = Math.min(Math.max(1, targetLength), span / 2);
    const rowTipFraction = data.chevron_tip_fraction
        ? data.chevron_tip_fraction[rowIndex]
        : tip_fraction;
    const endTip = x0Value + direction * span * rowTipFraction;
    const endBase = endTip - direction * length;
    return {
        xs: [endBase, endTip, endBase],
        ys: [
            yValue + targetYOffset,
            yValue,
            yValue - targetYOffset,
        ],
    };
}

const targetLength = chevronLength();
const targetYOffset = chevronYOffset();
const chevronXs = [];
const chevronYs = [];
data.x0.forEach(function (x0Value, rowIndex) {
    const vertices = chevronVertices(
        x0Value,
        data.x1[rowIndex],
        data.y[rowIndex],
        targetLength,
        targetYOffset,
        rowIndex
    );
    chevronXs.push(vertices.xs);
    chevronYs.push(vertices.ys);
});
data.chevron_xs = chevronXs;
data.chevron_ys = chevronYs;
source.change.emit();
