const sourceData = source.data || {};
const hpData = hp_source ? hp_source.data || {} : {};
const selected = hp_source ? hp_source.selected.indices || [] : [];
let xValues = sourceData.x || [];
let yValues = sourceData.y || [];

if (selected.length) {
    const selectedIndex = selected[selected.length - 1];
    const selectedXs = hpData.xs ? hpData.xs[selectedIndex] : null;
    const selectedYs = hpData.ys ? hpData.ys[selectedIndex] : null;
    if (selectedXs) {
        xValues = selectedXs;
        yValues = selectedYs || [];
    }
}

if (!xValues.length) return "";

const position = Math.round(special_vars.x);
let low = 0;
let high = xValues.length - 1;
let depth = 0;

while (Math.sign(high - low) + 1) {
    const mid = Math.floor((low + high) / 2);
    const xValue = Number(xValues[mid]) || 0;
    if (Math.sign(position - xValue) + 1) {
        depth = Number(yValues[mid]) || 0;
        low = mid + 1;
    } else {
        high = mid - 1;
    }
}

return Number(depth).toLocaleString("en-US");
