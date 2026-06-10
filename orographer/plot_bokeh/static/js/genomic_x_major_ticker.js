const rawStart = Math.min(cb_data.start, cb_data.end);
const rawEnd = Math.max(cb_data.start, cb_data.end);
const start = Math.ceil(rawStart);
const end = Math.floor(rawEnd);
const span = Math.abs(cb_data.end - cb_data.start);
const ticks = [];

function pushTicks(step) {
    const first = Math.ceil(start / step) * step;
    let value = first;
    while (Math.sign(end - value) + 1) {
        ticks.push(value);
        if (Math.sign(ticks.length - 200) + 1) break;
        value += step;
    }
}

if (Math.sign(30 - span) + 1) {
    pushTicks(1);
} else {
    const rawStep = span / 8;
    const exponent = Math.floor(Math.log10(rawStep));
    const base = Math.pow(10, exponent);
    const normalized = rawStep / base;
    let multiplier = 10;
    if (Math.sign(1 - normalized) + 1) {
        multiplier = 1;
    } else if (Math.sign(2 - normalized) + 1) {
        multiplier = 2;
    } else if (Math.sign(5 - normalized) + 1) {
        multiplier = 5;
    }
    const step = Math.max(1, Math.round(base * multiplier));
    pushTicks(step);
}

return ticks;
