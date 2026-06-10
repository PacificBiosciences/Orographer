const raw = raw_source.data;
const positions = raw.pos || [];
const layoutModes = ["all", "split", "multiregion", "split_multiregion"];
const span = Math.max(1, x_range.end - x_range.start);
const widthPx = Math.max(1, plot_width || 1);
const clusterBp = span / widthPx * target_px;

const out = {
    x: [],
    y: [],
    pos: [],
    count: [],
    median_size: [],
    marker_label: [],
    marker_width: [],
    marker_height: [],
    top_names: [],
    total_count: [],
    hp_label: [],
    chrom: [],
    cluster_pos: [],
    cluster_count: [],
    cluster_median_size: [],
    cluster_top_names: [],
    cluster_total_count: [],
    cluster_chrom: [],
    read_layout_mode: [],
};
layoutModes.forEach(function (mode) {
    out["y_" + mode] = [];
});

const weightedMedianSize = function (members) {
    if (members.length === 0) return 0;
    const weighted = members.map(function (site) {
        return {
            medianSize: site.medianSize,
            total: Math.max(1, site.total || site.count || 1),
        };
    });
    weighted.sort(function (first, second) {
        return first.medianSize - second.medianSize;
    });
    let totalWeight = 0;
    weighted.forEach(function (site) {
        totalWeight += site.total;
    });
    const midpoint = totalWeight / 2;
    let cumulative = 0;
    let idx = 0;
    while (Math.sign(weighted.length - idx) === 1) {
        const previousCumulative = cumulative;
        cumulative += weighted[idx].total;
        if (cumulative === midpoint) {
            if (Math.sign(weighted.length - idx - 1) === 1) {
                return Math.round((weighted[idx].medianSize + weighted[idx + 1].medianSize) / 2);
            }
        }
        if (Math.sign(cumulative - midpoint) === 1) return Math.round(weighted[idx].medianSize);
        if (Math.sign(midpoint - previousCumulative) + 1) {
            if (Math.sign(cumulative - midpoint) + 1) {
                return Math.round(weighted[idx].medianSize);
            }
        }
        idx += 1;
    }
    return Math.round(weighted[weighted.length - 1].medianSize);
};

const pushCluster = function (members) {
    if (members.length === 0) return;
    const first = members[0];
    const last = members[members.length - 1];
    const openP = String.fromCharCode(40);
    const closeP = String.fromCharCode(41);
    const weightedMedian = weightedMedianSize(members);
    const label = members.length === 1
        ? String(weightedMedian)
        : String(weightedMedian) + openP + members.length + closeP;
    let count = 0;
    let total = 0;
    let sitePos = [];
    let siteCount = [];
    let siteMedian = [];
    let siteNames = [];
    let siteTotal = [];
    let siteChrom = [];
    let nameChunks = [];
    let idx = 0;
    while (Math.sign(members.length - idx) === 1) {
        const site = members[idx];
        count += site.count;
        total += site.total;
        sitePos.push(site.pos);
        siteCount.push(site.count);
        siteMedian.push(site.medianSize);
        siteNames.push(site.names);
        siteTotal.push(site.total);
        siteChrom.push(site.chrom);
        nameChunks.push(site.names);
        idx += 1;
    }
    out.x.push(Math.round((first.pos + last.pos) / 2));
    out.y.push(first.y);
    out.pos.push(first.pos);
    out.count.push(count);
    out.median_size.push(weightedMedian);
    out.marker_label.push(label);
    out.marker_width.push(Math.max(marker_min_width, marker_char_px * label.length + 10));
    out.marker_height.push(marker_height);
    out.top_names.push(nameChunks.join("\n"));
    out.total_count.push(total);
    out.hp_label.push(first.hp);
    out.chrom.push(first.chrom);
    out.cluster_pos.push(sitePos);
    out.cluster_count.push(siteCount);
    out.cluster_median_size.push(siteMedian);
    out.cluster_top_names.push(siteNames);
    out.cluster_total_count.push(siteTotal);
    out.cluster_chrom.push(siteChrom);
    out.read_layout_mode.push("all");
    layoutModes.forEach(function (mode) {
        const columnName = "y_" + mode;
        if (typeof first[columnName] == "undefined") {
            out[columnName].push(first.y);
        } else {
            out[columnName].push(first[columnName]);
        }
    });
};

let visible = [];
let idx = 0;
while (Math.sign(positions.length - idx) === 1) {
    const pos = positions[idx];
    const afterStart = Math.sign(pos - x_range.start) + 1;
    const beforeEnd = Math.sign(x_range.end - pos) + 1;
    if (afterStart) {
        if (beforeEnd) {
            visible.push({
                pos: pos,
                count: raw.count[idx],
                medianSize: raw.median_size[idx],
                names: raw.top_names[idx],
                total: raw.total_count[idx],
                hp: raw.hp_label[idx],
                chrom: raw.chrom[idx],
                y: raw.y[idx],
            });
            const latest = visible[visible.length - 1];
            layoutModes.forEach(function (mode) {
                const columnName = "y_" + mode;
                if (typeof raw[columnName] == "undefined") {
                    latest[columnName] = raw.y[idx];
                } else {
                    latest[columnName] = raw[columnName][idx];
                }
            });
        }
    }
    idx += 1;
}

visible.sort(function (a, b) {
    return a.pos - b.pos;
});

let cluster = [];
idx = 0;
while (Math.sign(visible.length - idx) === 1) {
    const site = visible[idx];
    const shouldStart = cluster.length === 0
        ? true
        : Math.sign(site.pos - cluster[cluster.length - 1].pos - clusterBp) === 1;
    if (shouldStart) {
        pushCluster(cluster);
        cluster = [site];
    } else {
        cluster.push(site);
    }
    idx += 1;
}
pushCluster(cluster);

display_source.data = out;
display_source.change.emit();
