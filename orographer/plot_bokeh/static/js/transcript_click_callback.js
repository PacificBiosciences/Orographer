export default function (args) {
const source = args.source;
const arrow_source = args.arrow_source;
const label_source = args.label_source;
const transcript_label_source = args.transcript_label_source;
const intron_source = args.intron_source;
const intron_arrow_source = args.intron_arrow_source;
const feature_label_source = args.feature_label_source;
const coverage_source = args.coverage_source;
const total_coverage_source = args.total_coverage_source;
const coverage_figure = args.coverage_figure;
const plot_figure = args.plot_figure;
const page_status_div = args.page_status_div;
const page_first_button = args.page_first_button;
const page_prev_button = args.page_prev_button;
const page_next_button = args.page_next_button;
const page_last_button = args.page_last_button;
const page_target = args.page_target;
const page_delta = args.page_delta;
const indices = source.selected.indices || [];
const data = source.data;
const DEFAULT_PAGE_SIZE = 100;
const ISOFORM_SELECTED_Y = 0.55;
const ISOFORM_CONTEXT_START_Y = 1.25;
const ISOFORM_CONTEXT_STEP = 0.12;
const ISOFORM_CONTEXT_HEIGHT = 0.08;
const READ_CHEVRON_TARGET_PX = 14;
const READ_CHEVRON_TARGET_Y_PX = 5;
const READ_CHEVRON_TIP_FRACTION = 0.8;
const READ_CHEVRON_Y_OFFSET = 0.035;
const pageStateKey = plot_figure ? plot_figure.id : "isoseq";
if (!window.orographerIsoseqPageState) window.orographerIsoseqPageState = {};
if (!window.orographerIsoseqPageState[pageStateKey]) {
    window.orographerIsoseqPageState[pageStateKey] = {
        rowIndex: null,
        page: 0,
        requestedPage: 0,
        requestSerial: 0,
        pageSize: DEFAULT_PAGE_SIZE,
        assignedReadCount: 0,
    };
}
const pageState = window.orographerIsoseqPageState[pageStateKey];

function emptyLike(targetSource) {
    const emptyData = {};
    if (!targetSource) return emptyData;
    Object.keys(targetSource.data).forEach(function (key) {
        emptyData[key] = [];
    });
    return emptyData;
}

function replaceColumns(targetSource, incoming) {
    if (!targetSource) return;
    if (!incoming) {
        targetSource.data = emptyLike(targetSource);
        targetSource.selected.indices = [];
        targetSource.change.emit();
        return;
    }
    const nextData = {};
    Object.keys(incoming).forEach(function (key) {
        nextData[key] = incoming[key] || [];
    });
    targetSource.data = nextData;
    targetSource.selected.indices = [];
    targetSource.change.emit();
}

function currentChevronLength() {
    const xRange = plot_figure ? plot_figure.x_range : null;
    const width = plot_figure ? plot_figure.inner_width || plot_figure.width || 1 : 1;
    if (!xRange) return 1;
    const rangeSpan = Math.abs(xRange.end - xRange.start);
    return Math.max(1, rangeSpan / Math.max(1, width) * READ_CHEVRON_TARGET_PX);
}

function currentChevronYOffset() {
    const yRange = plot_figure ? plot_figure.y_range : null;
    const height = plot_figure ? plot_figure.inner_height || plot_figure.height || 1 : 1;
    if (!yRange) return READ_CHEVRON_Y_OFFSET;
    const rangeSpan = Math.abs(yRange.end - yRange.start);
    return Math.max(0.001, rangeSpan / Math.max(1, height) * READ_CHEVRON_TARGET_Y_PX);
}

function readChevronVertices(x0Value, x1Value, yValue) {
    const span = Math.abs(x1Value - x0Value);
    if (span == 0) return { xs: [], ys: [] };
    const yOffset = currentChevronYOffset();
    const direction = Math.sign(x1Value - x0Value) == 1 ? 1 : -1;
    const tipOffset = span * READ_CHEVRON_TIP_FRACTION;
    const length = Math.min(Math.max(1, currentChevronLength()), tipOffset);
    const endTip = x0Value + direction * tipOffset;
    const endBase = endTip - direction * length;
    return {
        xs: [endBase, endTip, endBase],
        ys: [
            yValue + yOffset,
            yValue,
            yValue - yOffset,
        ],
    };
}

function pushArrowChevronData(arrowData, x0Value, x1Value, yValue) {
    if (!arrowData.chevron_xs) arrowData.chevron_xs = [];
    if (!arrowData.chevron_ys) arrowData.chevron_ys = [];
    const vertices = readChevronVertices(x0Value, x1Value, yValue);
    arrowData.chevron_xs.push(vertices.xs);
    arrowData.chevron_ys.push(vertices.ys);
}

function updateCoverageRange(payload) {
    if (!coverage_figure) return;
    if (!coverage_figure.y_range) return;
    let maxDepth = 0;
    const yValues = payload ? payload.y || [] : [];
    yValues.forEach(function (value) {
        maxDepth = Math.max(maxDepth, Number(value) || 0);
    });
    coverage_figure.y_range.start = 0;
    coverage_figure.y_range.end = Math.max(1, maxDepth * 1.05);
}

function formatBlockCoordinates(chrom, x0Value, x1Value) {
    const start = Math.min(x0Value, x1Value);
    const end = Math.max(x0Value, x1Value);
    const span = Math.max(0, end - start);
    const openP = String.fromCharCode(40);
    const closeP = String.fromCharCode(41);
    const prefix = chrom ? chrom + ":" : "";
    return prefix + String(start) + "-" + String(end) + " "
        + openP + String(span) + " bp" + closeP;
}

function expandCompactReadPayload(payload) {
    if (!payload) return payload;
    if (payload.schema != "isoseq_compact_v1") return payload;
    const reads = payload.reads || {};
    const blocks = payload.blocks || {};
    const chrom = payload.chrom || "";
    const rowCount = (blocks.x0 || []).length;
    const arrowData = {
        x0: [],
        x1: [],
        y: [],
        y0: [],
        y1: [],
        color: [],
        read_name: [],
        segment_id: [],
        region_idx: [],
        row_index: [],
        alignment_order: [],
        fwd_read_start: [],
        fwd_read_end: [],
        source_kind: [],
        has_split_alignment: [],
        has_multiregion_connection: [],
        read_filter_alpha: [],
        arrow_line_alpha: [],
        arrowhead_alpha: [],
        arrow_selected_alpha: [],
        arrow_nonselected_alpha: [],
        gene_id: [],
        gene_name: [],
        transcript_id: [],
        group_id: [],
        angle: [],
        chevron_tip_fraction: [],
        chevron_xs: [],
        chevron_ys: [],
    };
    const labelData = {
        x: [],
        y: [],
        text: [],
        read_name: [],
        alignment_number: [],
        strand: [],
        coordinates: [],
        haplotype: [],
        sample_label: [],
        inclusion_reason: [],
        all_alignment_numbers: [],
        all_alignment_coordinates: [],
        has_split_alignment: [],
        has_multiregion_connection: [],
        gene_id: [],
        gene_name: [],
        transcript_id: [],
        annotation_label: [],
        read_filter_alpha: [],
        label_alpha: [],
    };
    (blocks.x0 || []).forEach(function (x0Value, rowIndex) {
        const readIndex = blocks.read_index[rowIndex] || 0;
        const x1Value = blocks.x1[rowIndex];
        const yValue = blocks.y[rowIndex];
        const alignmentNumber = blocks.alignment_order[rowIndex];
        arrowData.x0.push(x0Value);
        arrowData.x1.push(x1Value);
        arrowData.y.push(yValue);
        arrowData.y0.push(blocks.y0[rowIndex]);
        arrowData.y1.push(blocks.y1[rowIndex]);
        arrowData.color.push(blocks.color[rowIndex]);
        arrowData.read_name.push(reads.name[readIndex]);
        arrowData.segment_id.push(blocks.segment_id[rowIndex]);
        arrowData.region_idx.push(0);
        arrowData.row_index.push(0);
        arrowData.alignment_order.push(alignmentNumber);
        arrowData.fwd_read_start.push(blocks.fwd_read_start[rowIndex]);
        arrowData.fwd_read_end.push(blocks.fwd_read_end[rowIndex]);
        arrowData.source_kind.push("arrow");
        arrowData.has_split_alignment.push(false);
        arrowData.has_multiregion_connection.push(false);
        arrowData.read_filter_alpha.push(1);
        arrowData.arrow_line_alpha.push(0.5);
        arrowData.arrowhead_alpha.push(0.65);
        arrowData.arrow_selected_alpha.push(1);
        arrowData.arrow_nonselected_alpha.push(0.12);
        arrowData.chevron_tip_fraction.push(0.8);
        arrowData.gene_id.push(reads.gene_id[readIndex]);
        arrowData.gene_name.push(reads.gene_name[readIndex]);
        arrowData.transcript_id.push(reads.transcript_id[readIndex]);
        arrowData.group_id.push(reads.group_id[readIndex]);
        arrowData.angle.push(Math.sign(x1Value - x0Value) == 1 ? -Math.PI / 2 : Math.PI / 2);
        pushArrowChevronData(arrowData, x0Value, x1Value, yValue);
        labelData.x.push((x0Value + x1Value) / 2);
        labelData.y.push(yValue);
        labelData.text.push(String(alignmentNumber));
        labelData.read_name.push(reads.name[readIndex]);
        labelData.alignment_number.push(alignmentNumber);
        labelData.strand.push(reads.strand[readIndex]);
        labelData.coordinates.push(formatBlockCoordinates(chrom, x0Value, x1Value));
        labelData.haplotype.push(reads.haplotype[readIndex]);
        labelData.sample_label.push(reads.sample_label[readIndex]);
        labelData.inclusion_reason.push("");
        labelData.all_alignment_numbers.push(reads.all_alignment_numbers[readIndex]);
        labelData.all_alignment_coordinates.push(reads.all_alignment_coordinates[readIndex]);
        labelData.has_split_alignment.push(false);
        labelData.has_multiregion_connection.push(false);
        labelData.gene_id.push(reads.gene_id[readIndex]);
        labelData.gene_name.push(reads.gene_name[readIndex]);
        labelData.transcript_id.push(reads.transcript_id[readIndex]);
        labelData.annotation_label.push(reads.annotation_label ? reads.annotation_label[readIndex] : "");
        labelData.read_filter_alpha.push(1);
        labelData.label_alpha.push(0.8);
    });
    payload.arrow_data = arrowData;
    payload.label_data = labelData;
    payload.loaded_read_count = payload.loaded_read_count || rowCount;
    return payload;
}

function haplotypeLabel(haplotypeValue) {
    return haplotypeValue ? "HP:" + String(haplotypeValue) : "Unassigned";
}

function haplotypeColor(haplotypeValue) {
    const palette = [
        "#1F77B4",
        "#FF7F0E",
        "#2CA02C",
        "#D62728",
        "#9467BD",
        "#8C564B",
        "#E377C2",
        "#7F7F7F",
        "#BCBD22",
        "#17BECF",
    ];
    if (!haplotypeValue) return "#888888";
    const index = Number(haplotypeValue) - 1;
    return palette[index % palette.length];
}

function expandShardedReadPayload(manifest, group, page, records) {
    const readHeight = 0.24;
    const arrowData = emptyLike(arrow_source);
    const labelData = emptyLike(label_source);
    records.forEach(function (record, readRank) {
        const readName = record[0];
        const geneId = record[1];
        const geneName = record[2];
        const transcriptId = record[3];
        const groupId = record[4];
        const haplotypeValue = record[5];
        const isFwd = record[6];
        const fwdReadStart = record[7];
        const fwdReadEnd = record[8];
        const blocks = record[9] || [];
        const yValue = manifest.selected_read_y_start + readRank * readHeight + readHeight / 2;
        const fullBlockRows = [];
        blocks.forEach(function (block) {
            const start = block[0];
            const end = block[1];
            fullBlockRows.push({
                number: block[2],
                coordinate: formatBlockCoordinates(manifest.chrom || "", start, end),
            });
            if (Math.sign(end - manifest.coordinate_start) != 1) return;
            if (Math.sign(manifest.coordinate_end - start) != 1) return;
            const plotStart = Math.max(start, manifest.coordinate_start);
            const plotEnd = Math.min(end, manifest.coordinate_end);
            const x0Value = isFwd ? plotStart : plotEnd;
            const x1Value = isFwd ? plotEnd : plotStart;
            const alignmentNumber = block[2];
            const segmentId = "r" + String(manifest.region_idx)
                + ":row" + String(manifest.row_index)
                + ":" + readName
                + ":1:" + String(manifest.chrom)
                + ":" + String(start)
                + "-" + String(end)
                + ":block" + String(alignmentNumber);
            arrowData.x0.push(x0Value);
            arrowData.x1.push(x1Value);
            arrowData.y.push(yValue);
            if (arrowData.original_y) arrowData.original_y.push(yValue);
            arrowData.y0.push(yValue);
            arrowData.y1.push(yValue);
            arrowData.color.push(haplotypeColor(haplotypeValue));
            arrowData.read_name.push(readName);
            arrowData.layout_read_name.push(readName);
            arrowData.segment_id.push(segmentId);
            arrowData.region_idx.push(manifest.region_idx);
            arrowData.row_index.push(manifest.row_index);
            arrowData.alignment_order.push(alignmentNumber);
            arrowData.fwd_read_start.push(fwdReadStart);
            arrowData.fwd_read_end.push(fwdReadEnd);
            arrowData.source_kind.push("arrow");
            arrowData.has_split_alignment.push(false);
            arrowData.has_multiregion_connection.push(false);
            arrowData.read_filter_alpha.push(1);
            arrowData.arrow_line_alpha.push(0.5);
            arrowData.arrowhead_alpha.push(0.65);
            arrowData.arrow_selected_alpha.push(1);
            arrowData.arrow_nonselected_alpha.push(0.12);
            arrowData.chevron_tip_fraction.push(0.8);
            arrowData.gene_id.push(geneId);
            arrowData.gene_name.push(geneName);
            arrowData.transcript_id.push(transcriptId);
            arrowData.group_id.push(groupId);
            if (arrowData.read_filter_visible_all) arrowData.read_filter_visible_all.push(true);
            if (arrowData.read_filter_visible_split) arrowData.read_filter_visible_split.push(true);
            if (arrowData.read_filter_visible_multiregion) {
                arrowData.read_filter_visible_multiregion.push(true);
            }
            if (arrowData.read_filter_visible_split_multiregion) {
                arrowData.read_filter_visible_split_multiregion.push(true);
            }
            arrowData.angle.push(Math.sign(x1Value - x0Value) == 1 ? -Math.PI / 2 : Math.PI / 2);
            pushArrowChevronData(arrowData, x0Value, x1Value, yValue);
            labelData.x.push((plotStart + plotEnd) / 2);
            labelData.y.push(yValue);
            if (labelData.original_y) labelData.original_y.push(yValue);
            labelData.text.push(String(alignmentNumber));
            labelData.read_name.push(readName);
            labelData.layout_read_name.push(readName);
            labelData.alignment_number.push(alignmentNumber);
            labelData.strand.push(isFwd ? "Forward (+)" : "Reverse (-)");
            labelData.coordinates.push(
                formatBlockCoordinates(manifest.chrom || "", plotStart, plotEnd)
            );
            labelData.haplotype.push(haplotypeLabel(haplotypeValue));
            labelData.sample_label.push(manifest.sample_label || "");
            labelData.inclusion_reason.push("");
            labelData.all_alignment_numbers.push([]);
            labelData.all_alignment_coordinates.push([]);
            labelData.has_split_alignment.push(false);
            labelData.has_multiregion_connection.push(false);
            labelData.gene_id.push(geneId);
            labelData.gene_name.push(geneName);
            labelData.transcript_id.push(transcriptId);
            labelData.annotation_label.push(manifest.annotation_label || "");
            labelData.read_filter_alpha.push(1);
            labelData.label_alpha.push(0.8);
        });
        fullBlockRows.sort(function (left, right) {
            return left.number - right.number;
        });
        const numbers = fullBlockRows.map(function (row) {
            return row.number;
        });
        const coordinates = fullBlockRows.map(function (row) {
            return row.coordinate;
        });
        let rowIndex = 0;
        while (Math.sign(labelData.read_name.length - rowIndex) == 1) {
            if (labelData.read_name[rowIndex] == readName) {
                labelData.all_alignment_numbers[rowIndex] = numbers;
                labelData.all_alignment_coordinates[rowIndex] = coordinates;
            }
            rowIndex += 1;
        }
    });
    return {
        arrow_data: arrowData,
        label_data: labelData,
        page: page,
        page_size: manifest.page_size,
        loaded_read_count: records.length,
        assigned_read_count: group.assigned_read_count,
    };
}

function parseJsonResponse(response, url) {
    if (url.slice(-3) != ".gz") return response.json();
    if (typeof DecompressionStream == "undefined") {
        throw new Error("This browser cannot decompress gzip IsoSeq chunks");
    }
    return response.blob()
        .then(function (blob) {
            const stream = blob.stream().pipeThrough(new DecompressionStream("gzip"));
            return new Response(stream).text();
        })
        .then(function (text) {
            return JSON.parse(text);
        });
}

function inClosedRange(value, low, high) {
    const aboveLow = Math.sign(value - low) + 1;
    const belowHigh = Math.sign(high - value) + 1;
    return aboveLow ? Boolean(belowHigh) : false;
}

function unpackMessagePack(buffer) {
    const bytes = new Uint8Array(buffer);
    const decoder = new TextDecoder();
    let offset = 0;

    function readByte() {
        const value = bytes[offset];
        offset += 1;
        return value;
    }

    function readUnsigned(byteCount) {
        let value = 0;
        let index = 0;
        while (Math.sign(byteCount - index) == 1) {
            value = value * 256 + readByte();
            index += 1;
        }
        return value;
    }

    function readSigned(byteCount) {
        const value = readUnsigned(byteCount);
        const maxPositive = Math.pow(2, byteCount * 8 - 1);
        const fullRange = Math.pow(2, byteCount * 8);
        return Math.sign(value - maxPositive) + 1 ? value - fullRange : value;
    }

    function readString(length) {
        const start = offset;
        offset += length;
        return decoder.decode(bytes.slice(start, offset));
    }

    function readArray(length) {
        const values = [];
        let index = 0;
        while (Math.sign(length - index) == 1) {
            values.push(readValue());
            index += 1;
        }
        return values;
    }

    function readMap(length) {
        const values = {};
        let index = 0;
        while (Math.sign(length - index) == 1) {
            const key = readValue();
            values[key] = readValue();
            index += 1;
        }
        return values;
    }

    function readValue() {
        const prefix = readByte();
        if (inClosedRange(prefix, 0, 127)) return prefix;
        if (inClosedRange(prefix, 128, 143)) return readMap(prefix - 128);
        if (inClosedRange(prefix, 144, 159)) return readArray(prefix - 144);
        if (inClosedRange(prefix, 160, 191)) return readString(prefix - 160);
        if (inClosedRange(prefix, 224, 255)) return prefix - 256;
        if (prefix == 192) return null;
        if (prefix == 194) return false;
        if (prefix == 195) return true;
        if (prefix == 204) return readUnsigned(1);
        if (prefix == 205) return readUnsigned(2);
        if (prefix == 206) return readUnsigned(4);
        if (prefix == 207) return readUnsigned(8);
        if (prefix == 208) return readSigned(1);
        if (prefix == 209) return readSigned(2);
        if (prefix == 210) return readSigned(4);
        if (prefix == 211) return readSigned(8);
        if (prefix == 217) return readString(readUnsigned(1));
        if (prefix == 218) return readString(readUnsigned(2));
        if (prefix == 219) return readString(readUnsigned(4));
        if (prefix == 220) return readArray(readUnsigned(2));
        if (prefix == 221) return readArray(readUnsigned(4));
        if (prefix == 222) return readMap(readUnsigned(2));
        if (prefix == 223) return readMap(readUnsigned(4));
        throw new Error("Unsupported MessagePack prefix " + String(prefix));
    }

    return readValue();
}

function parseBinaryResponse(response, url) {
    if (url.slice(-3) != ".gz") return response.arrayBuffer();
    if (typeof DecompressionStream == "undefined") {
        throw new Error("This browser cannot decompress gzip IsoSeq shards");
    }
    return response.blob()
        .then(function (blob) {
            const stream = blob.stream().pipeThrough(new DecompressionStream("gzip"));
            return new Response(stream).arrayBuffer();
        });
}

function loadIsoseqManifest(url) {
    if (!window.orographerIsoseqChunkCache) window.orographerIsoseqChunkCache = {};
    if (window.orographerIsoseqChunkCache[url]) {
        return Promise.resolve(window.orographerIsoseqChunkCache[url]);
    }
    return fetch(url)
        .then(function (response) {
            if (!response.ok) {
                throw new Error("Could not load IsoSeq read manifest: " + url);
            }
            return parseJsonResponse(response, url);
        })
        .then(function (manifest) {
            window.orographerIsoseqChunkCache[url] = manifest;
            prefetchAllReadAssignments(manifest.all_assignments_url);
            return manifest;
        });
}

function prefetchAllReadAssignments(assignmentsUrl) {
    if (!assignmentsUrl) return;
    if (!window.orographerAllReadAssignments) window.orographerAllReadAssignments = {};
    if (window.orographerAllReadAssignments[assignmentsUrl]) return;
    window.orographerAllReadAssignments[assignmentsUrl] = "loading";
    fetch(assignmentsUrl)
        .then(function (response) {
            if (!response.ok) return null;
            return parseJsonResponse(response, assignmentsUrl);
        })
        .then(function (payload) {
            if (payload) {
                window.orographerAllReadAssignments[assignmentsUrl] = payload.assignments || {};
            } else {
                window.orographerAllReadAssignments[assignmentsUrl] = {};
            }
        })
        .catch(function () {
            window.orographerAllReadAssignments[assignmentsUrl] = {};
        });
}

function shardUrl(manifestUrl, shardFile) {
    const slashIndex = manifestUrl.lastIndexOf("/");
    if (slashIndex == -1) return shardFile;
    return manifestUrl.slice(0, slashIndex + 1) + shardFile;
}

function loadIsoseqShard(manifestUrl, manifest, shardId) {
    const url = shardUrl(manifestUrl, manifest.shards[shardId]);
    if (!window.orographerIsoseqChunkCache) window.orographerIsoseqChunkCache = {};
    if (window.orographerIsoseqChunkCache[url]) {
        return Promise.resolve(window.orographerIsoseqChunkCache[url]);
    }
    return fetch(url)
        .then(function (response) {
            if (!response.ok) {
                throw new Error("Could not load IsoSeq read shard: " + url);
            }
            return parseBinaryResponse(response, url);
        })
        .then(function (buffer) {
            const shard = unpackMessagePack(buffer);
            window.orographerIsoseqChunkCache[url] = shard;
            return shard;
        });
}

function uniqueReadCount(targetSource) {
    if (!targetSource) return 0;
    const targetData = targetSource.data;
    const names = {};
    (targetData.read_name || []).forEach(function (readName) {
        if (readName) {
            names[readName] = true;
        }
    });
    return Object.keys(names).length;
}

function selectedReadCount(selectedReads) {
    if (!selectedReads) return 0;
    return Object.keys(selectedReads.names || {}).length;
}

function resizeSelectedView(readCount) {
    if (!plot_figure) return;
    const count = Math.max(1, Number(readCount) || 0);
    const readStep = estimateReadStep(arrow_source);
    let topY = null;
    if (arrow_source) {
        const yValues = arrow_source.data.original_y || arrow_source.data.y || [];
        yValues.forEach(function (yValue) {
            if (topY == null) {
                topY = yValue;
            } else if (Math.sign(topY - yValue) == 1) {
                topY = yValue;
            }
        });
    }
    if (topY == null) {
        topY = 1;
    }
    plot_figure.y_range.start = topY + count * readStep + readStep * 3;
    plot_figure.height = Math.min(900, Math.max(260, 160 + count * 9));
}

function selectedReadNames() {
    const names = {};
    let selectedCount = 0;
    if (arrow_source) {
        const arrowData = arrow_source.data;
        const selected = arrow_source.selected.indices || [];
        selected.forEach(function (rowIndex) {
            const readName = arrowData.read_name ? arrowData.read_name[rowIndex] : "";
            if (readName) {
                names[readName] = true;
                selectedCount += 1;
            }
        });
    }
    if (label_source) {
        const labelData = label_source.data;
        const selected = label_source.selected.indices || [];
        selected.forEach(function (rowIndex) {
            const readName = labelData.read_name ? labelData.read_name[rowIndex] : "";
            if (readName) {
                names[readName] = true;
                selectedCount += 1;
            }
        });
    }
    return { names: names, selectedCount: selectedCount };
}

function ensureOriginalY(targetSource) {
    if (!targetSource) return;
    const targetData = targetSource.data;
    if (!targetData.y) return;
    const originalY = targetData.original_y || [];
    if (originalY.length - targetData.y.length) {
        targetData.original_y = targetData.y.map(function (value) {
            return value;
        });
    }
}

function ensureOriginalTranscriptGeometry() {
    const originalTop = data.original_top || [];
    if (originalTop.length - data.top.length) {
        data.original_top = data.top.map(function (value) {
            return value;
        });
    }
    const originalBottom = data.original_bottom || [];
    if (originalBottom.length - data.bottom.length) {
        data.original_bottom = data.bottom.map(function (value) {
            return value;
        });
    }
}

function ensureOriginalYs(targetSource) {
    if (!targetSource) return;
    const targetData = targetSource.data;
    if (!targetData.ys) return;
    const originalYs = targetData.original_ys || [];
    if (originalYs.length - targetData.ys.length) {
        targetData.original_ys = targetData.ys.map(function (values) {
            return values.slice();
        });
    }
}

function estimateReadStep(targetSource) {
    ensureOriginalY(targetSource);
    if (!targetSource) return 0.24;
    const targetData = targetSource.data;
    if (!targetData.original_y) return 0.24;
    const seen = {};
    const yValues = [];
    targetData.original_y.forEach(function (yValue) {
        const key = String(yValue);
        if (!seen[key]) {
            seen[key] = true;
            yValues.push(yValue);
        }
    });
    yValues.sort(function (left, right) {
        return left - right;
    });
    let previousY = null;
    let bestStep = null;
    yValues.forEach(function (yValue) {
        if (previousY != null) {
            const diff = yValue - previousY;
            if (diff) {
                if (bestStep == null) {
                    bestStep = diff;
                } else if (Math.sign(bestStep - diff) == 1) {
                    bestStep = diff;
                }
            }
        }
        previousY = yValue;
    });
    return bestStep == null ? 0.24 : bestStep;
}

function packedYByRead(targetSource, selectedReads) {
    ensureOriginalY(targetSource);
    const yByRead = {};
    if (!targetSource) return yByRead;
    const targetData = targetSource.data;
    if (!targetData.original_y) return yByRead;
    const readRows = {};
    targetData.original_y.forEach(function (yValue, rowIndex) {
        const readName = targetData.read_name ? targetData.read_name[rowIndex] : "";
        if (!selectedReads.names[readName]) return;
        if (!readRows[readName]) {
            readRows[readName] = { readName: readName, y: yValue };
            return;
        }
        if (Math.sign(readRows[readName].y - yValue) == 1) {
            readRows[readName].y = yValue;
        }
    });
    const rows = Object.keys(readRows).map(function (readName) {
        return readRows[readName];
    });
    rows.sort(function (left, right) {
        const yDiff = left.y - right.y;
        if (yDiff) return yDiff;
        return left.readName.localeCompare(right.readName);
    });
    if (!rows.length) return yByRead;
    let topY = null;
    targetData.original_y.forEach(function (yValue) {
        if (topY == null) {
            topY = yValue;
        } else if (Math.sign(topY - yValue) == 1) {
            topY = yValue;
        }
    });
    if (topY == null) return yByRead;
    const readStep = estimateReadStep(targetSource);
    rows.forEach(function (row, rank) {
        yByRead[row.readName] = topY + rank * readStep;
    });
    return yByRead;
}

function sortedIsoformRows(selectedTranscriptId) {
    ensureOriginalTranscriptGeometry();
    const rows = {};
    data.transcript_id.forEach(function (transcriptId, rowIndex) {
        const idText = String(transcriptId);
        if (idText == String(selectedTranscriptId)) return;
        const center = (data.original_top[rowIndex] + data.original_bottom[rowIndex]) / 2;
        if (!rows[idText]) {
            rows[idText] = { transcriptId: idText, y: center };
            return;
        }
        if (Math.sign(rows[idText].y - center) == 1) {
            rows[idText].y = center;
        }
    });
    const sortedRows = Object.keys(rows).map(function (transcriptId) {
        return rows[transcriptId];
    });
    sortedRows.sort(function (left, right) {
        const yDiff = left.y - right.y;
        if (yDiff) return yDiff;
        return left.transcriptId.localeCompare(right.transcriptId);
    });
    return sortedRows;
}

function denseIsoformYByTranscript(selectedTranscriptId) {
    const yByTranscript = {};
    const rows = sortedIsoformRows(selectedTranscriptId);
    const readStart = data.selected_read_y_start ? data.selected_read_y_start[0] : 0;
    const contextEnd = Math.max(ISOFORM_CONTEXT_START_Y, readStart - 0.3);
    const contextSpan = Math.max(0, contextEnd - ISOFORM_CONTEXT_START_Y);
    yByTranscript[String(selectedTranscriptId)] = ISOFORM_SELECTED_Y;
    rows.forEach(function (row, rank) {
        let center = ISOFORM_CONTEXT_START_Y + rank * ISOFORM_CONTEXT_STEP;
        if (rows.length == 1) {
            center = ISOFORM_CONTEXT_START_Y + contextSpan / 2;
        } else if (Math.sign(rows.length - 1) == 1) {
            center = ISOFORM_CONTEXT_START_Y + rank * contextSpan / (rows.length - 1);
        }
        yByTranscript[row.transcriptId] = center;
    });
    return yByTranscript;
}

function applyTranscriptGeometry(selectedTranscriptId, compactIsoforms) {
    ensureOriginalTranscriptGeometry();
    if (!compactIsoforms) {
        data.top = data.original_top.map(function (value) {
            return value;
        });
        data.bottom = data.original_bottom.map(function (value) {
            return value;
        });
        return;
    }
    const yByTranscript = denseIsoformYByTranscript(selectedTranscriptId);
    data.top = data.original_top.map(function (originalTop, rowIndex) {
        const transcriptId = String(data.transcript_id[rowIndex]);
        const center = yByTranscript[transcriptId];
        if (transcriptId == String(selectedTranscriptId)) {
            const height = data.original_bottom[rowIndex] - originalTop;
            return center - height / 2;
        }
        return center - ISOFORM_CONTEXT_HEIGHT / 2;
    });
    data.bottom = data.original_bottom.map(function (originalBottom, rowIndex) {
        const transcriptId = String(data.transcript_id[rowIndex]);
        const center = yByTranscript[transcriptId];
        if (transcriptId == String(selectedTranscriptId)) {
            const height = originalBottom - data.original_top[rowIndex];
            return center + height / 2;
        }
        return center + ISOFORM_CONTEXT_HEIGHT / 2;
    });
}

function applySourceYGeometry(targetSource, selectedTranscriptId, compactIsoforms) {
    ensureOriginalY(targetSource);
    if (!targetSource) return;
    const targetData = targetSource.data;
    if (!targetData.y) return;
    if (!targetData.original_y) return;
    if (!compactIsoforms) {
        targetData.y = targetData.original_y.map(function (value) {
            return value;
        });
        return;
    }
    const yByTranscript = denseIsoformYByTranscript(selectedTranscriptId);
    targetData.y = targetData.original_y.map(function (originalY, rowIndex) {
        const transcriptId = String(targetData.transcript_id[rowIndex]);
        if (transcriptId == String(selectedTranscriptId)) return ISOFORM_SELECTED_Y;
        return yByTranscript[transcriptId] || originalY;
    });
}

function applyIntronGeometry(selectedTranscriptId, compactIsoforms) {
    ensureOriginalYs(intron_source);
    if (!intron_source) return;
    const intronData = intron_source.data;
    if (!intronData.ys) return;
    if (!intronData.original_ys) return;
    if (!compactIsoforms) {
        intronData.ys = intronData.original_ys.map(function (values) {
            return values.slice();
        });
        return;
    }
    const yByTranscript = denseIsoformYByTranscript(selectedTranscriptId);
    intronData.ys = intronData.original_ys.map(function (values, rowIndex) {
        const transcriptId = String(intronData.transcript_id[rowIndex]);
        const yValue = yByTranscript[transcriptId] || values[0];
        return [yValue, yValue];
    });
}

function applyPackedY(targetSource, yByRead, hideReads, hasReadSelection) {
    ensureOriginalY(targetSource);
    if (!targetSource) return;
    const targetData = targetSource.data;
    if (!targetData.y) return;
    if (!targetData.original_y) return;
    targetData.y = targetData.original_y.map(function (originalY, rowIndex) {
        if (!hideReads) return originalY;
        if (!hasReadSelection) return originalY;
        const readName = targetData.read_name ? targetData.read_name[rowIndex] : "";
        return typeof yByRead[readName] == "undefined" ? originalY : yByRead[readName];
    });
}

function hideUnselectedReadsActive() {
    return window.orographerHideUnselectedReads == true;
}

function hideUnselectedIsoformsActive() {
    return window.orographerHideUnselectedIsoforms == true;
}

function applyReadVisibility() {
    const hideReads = hideUnselectedReadsActive();
    const selectedReads = selectedReadNames();
    const hasReadSelection = Boolean(selectedReads.selectedCount);
    const yByRead = packedYByRead(arrow_source, selectedReads);
    if (arrow_source) {
        const arrowData = arrow_source.data;
        applyPackedY(arrow_source, yByRead, hideReads, hasReadSelection);
        if (arrowData.arrow_line_alpha) {
            arrowData.arrow_line_alpha = arrowData.arrow_line_alpha.map(
                function (_oldValue, rowIndex) {
                    if (!hideReads) return 0.5;
                    if (!hasReadSelection) return 0.5;
                    const readName = arrowData.read_name ? arrowData.read_name[rowIndex] : "";
                    return selectedReads.names[readName] ? 0.5 : 0;
                }
            );
        }
        if (arrowData.arrowhead_alpha) {
            arrowData.arrowhead_alpha = arrowData.arrowhead_alpha.map(
                function (_oldValue, rowIndex) {
                    if (!hideReads) return 0.65;
                    if (!hasReadSelection) return 0.65;
                    const readName = arrowData.read_name ? arrowData.read_name[rowIndex] : "";
                    return selectedReads.names[readName] ? 0.65 : 0;
                }
            );
        }
        if (arrowData.arrow_nonselected_alpha) {
            arrowData.arrow_nonselected_alpha = arrowData.arrow_nonselected_alpha.map(
                function (_oldValue, rowIndex) {
                    if (!hideReads) return 0.12;
                    if (!hasReadSelection) return 0.12;
                    const readName = arrowData.read_name ? arrowData.read_name[rowIndex] : "";
                    return selectedReads.names[readName] ? 0.12 : 0;
                }
            );
        }
        arrow_source.change.emit();
    }
    if (label_source) {
        const labelData = label_source.data;
        applyPackedY(label_source, yByRead, hideReads, hasReadSelection);
        if (labelData.label_alpha) {
            labelData.label_alpha = labelData.label_alpha.map(function (_oldValue, rowIndex) {
                if (!hideReads) return 0.8;
                if (!hasReadSelection) return 0.8;
                const readName = labelData.read_name ? labelData.read_name[rowIndex] : "";
                return selectedReads.names[readName] ? 0.8 : 0;
            });
            label_source.change.emit();
        }
    }
    if (hideReads) {
        if (hasReadSelection) {
            resizeSelectedView(selectedReadCount(selectedReads));
        } else {
            resizeSelectedView(uniqueReadCount(arrow_source));
        }
    } else {
        resizeSelectedView(uniqueReadCount(arrow_source));
    }
}

function setTranscriptVisibility(selectedTranscriptId) {
    const selected = selectedTranscriptId != null;
    const hideIsoforms = hideUnselectedIsoformsActive();
    const compactIsoforms = selected ? !hideIsoforms : false;
    applyTranscriptGeometry(selectedTranscriptId, compactIsoforms);
    if (data.alpha) {
        data.alpha = data.alpha.map(function (_oldValue, rowIndex) {
            const baseAlpha = data.base_alpha ? data.base_alpha[rowIndex] : 0.86;
            if (!selected) return baseAlpha;
            if (!hideIsoforms) {
                return String(data.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? baseAlpha
                    : 0.24;
            }
            return String(data.transcript_id[rowIndex]) == String(selectedTranscriptId)
                ? baseAlpha
                : 0;
        });
    }
    if (data.line_alpha) {
        data.line_alpha = data.line_alpha.map(function (_oldValue, rowIndex) {
            const baseAlpha = data.base_line_alpha ? data.base_line_alpha[rowIndex] : 1;
            if (!selected) return baseAlpha;
            if (!hideIsoforms) {
                return String(data.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? baseAlpha
                    : 0.28;
            }
            return String(data.transcript_id[rowIndex]) == String(selectedTranscriptId)
                ? baseAlpha
                : 0;
        });
    }
    source.change.emit();

    if (transcript_label_source) {
        applySourceYGeometry(transcript_label_source, selectedTranscriptId, compactIsoforms);
        const labelData = transcript_label_source.data;
        if (labelData.alpha) {
            labelData.alpha = labelData.alpha.map(function (_oldValue, rowIndex) {
                const baseAlpha = labelData.base_alpha ? labelData.base_alpha[rowIndex] : 1;
                if (!selected) return baseAlpha;
                if (!hideIsoforms) {
                    return String(labelData.transcript_id[rowIndex]) == String(selectedTranscriptId)
                        ? baseAlpha
                        : 0.18;
                }
                return String(labelData.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? baseAlpha
                    : 0;
            });
            transcript_label_source.change.emit();
        }
    }

    if (intron_source) {
        applyIntronGeometry(selectedTranscriptId, compactIsoforms);
        const intronData = intron_source.data;
        if (intronData.line_alpha) {
            intronData.line_alpha = intronData.line_alpha.map(function (_oldValue, rowIndex) {
                if (!selected) return 0.75;
                if (!hideIsoforms) {
                    return String(intronData.transcript_id[rowIndex])
                        == String(selectedTranscriptId)
                        ? 0.75
                        : 0.22;
                }
                return String(intronData.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? 0.75
                    : 0;
            });
            intron_source.change.emit();
        }
    }

    if (intron_arrow_source) {
        applySourceYGeometry(intron_arrow_source, selectedTranscriptId, compactIsoforms);
        const arrowData = intron_arrow_source.data;
        if (arrowData.alpha) {
            arrowData.alpha = arrowData.alpha.map(function (_oldValue, rowIndex) {
                if (!selected) return 0.85;
                if (!hideIsoforms) {
                    return String(arrowData.transcript_id[rowIndex]) == String(selectedTranscriptId)
                        ? 0.85
                        : 0.22;
                }
                return String(arrowData.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? 0.85
                    : 0;
            });
            intron_arrow_source.change.emit();
        }
    }

    if (feature_label_source) {
        applySourceYGeometry(feature_label_source, selectedTranscriptId, compactIsoforms);
        const featureData = feature_label_source.data;
        if (featureData.label_alpha) {
            featureData.label_alpha = featureData.label_alpha.map(function (_oldValue, rowIndex) {
                if (!selected) return 0.8;
                if (!hideIsoforms) {
                    return String(featureData.transcript_id[rowIndex])
                        == String(selectedTranscriptId)
                        ? 0.8
                        : 0;
                }
                return String(featureData.transcript_id[rowIndex]) == String(selectedTranscriptId)
                    ? 0.8
                    : 0;
            });
            feature_label_source.change.emit();
        }
    }
}

function restoreTranscriptOnlyView() {
    setTranscriptVisibility(null);
    replaceColumns(arrow_source, null);
    replaceColumns(label_source, null);
    if (total_coverage_source) {
        replaceColumns(coverage_source, total_coverage_source.data);
        updateCoverageRange(total_coverage_source.data);
    }
    if (coverage_figure) coverage_figure.visible = true;
    if (plot_figure) {
        if (data.transcript_area_height) {
            if (data.transcript_area_height.length) {
        plot_figure.y_range.start = data.transcript_area_height[0];
        plot_figure.height = data.transcript_plot_height[0];
            }
        }
    }
    pageState.rowIndex = null;
    pageState.page = 0;
    pageState.requestedPage = 0;
    pageState.requestSerial += 1;
    pageState.assignedReadCount = 0;
    updatePageControls(null);
}

function updatePageControls(payload) {
    const noSelection = pageState.rowIndex == null;
    const currentFirstPage = pageState.page == 0;
    const assigned = payload ? payload.assigned_read_count || 0 : pageState.assignedReadCount;
    const pageSize = pageState.pageSize || DEFAULT_PAGE_SIZE;
    const nextStart = (pageState.page + 1) * pageSize;
    const lastPageReached = Boolean(Math.sign(nextStart - assigned) + 1);
    if (typeof page_status_div != "undefined") {
        if (page_status_div) {
            let text = "";
            if (payload) {
                const loaded = payload.loaded_read_count || 0;
                const assignedReads = payload.assigned_read_count || loaded;
                const startRead = loaded ? pageState.page * pageSize + 1 : 0;
                const endRead = Math.min(assignedReads, pageState.page * pageSize + loaded);
                const startText = Number(startRead).toLocaleString("en-US");
                const endText = Number(endRead).toLocaleString("en-US");
                const assignedText = Number(assignedReads).toLocaleString("en-US");
                text = loaded
                    ? startText + "-" + endText + " of " + assignedText + " reads"
                    : "";
            }
            page_status_div.text = text;
        }
    }
    if (typeof page_first_button != "undefined") {
        page_first_button.disabled = Boolean(noSelection || currentFirstPage);
    }
    if (typeof page_prev_button != "undefined") {
        page_prev_button.disabled = Boolean(noSelection || currentFirstPage);
    }
    if (typeof page_next_button != "undefined") {
        page_next_button.disabled = Boolean(noSelection || lastPageReached);
    }
    if (typeof page_last_button != "undefined") {
        page_last_button.disabled = Boolean(noSelection || lastPageReached);
    }
}

function loadShardedReadPage(manifestUrl, chunkKey, page) {
    return loadIsoseqManifest(manifestUrl)
        .then(function (manifest) {
            const group = manifest.groups ? manifest.groups[chunkKey] : null;
            if (!group) return null;
            const pageSize = manifest.page_size || DEFAULT_PAGE_SIZE;
            const startIndex = page * pageSize;
            const endIndex = Math.min(group.read_ids.length, startIndex + pageSize);
            const readIds = group.read_ids.slice(startIndex, endIndex);
            const shardIds = {};
            readIds.forEach(function (readId) {
                const shardId = Math.floor(readId / manifest.shard_size);
                shardIds[shardId] = true;
            });
            const shardPromises = Object.keys(shardIds).map(function (shardIdText) {
                return loadIsoseqShard(manifestUrl, manifest, Number(shardIdText));
            });
            return Promise.all(shardPromises).then(function (shards) {
                const recordsById = {};
                shards.forEach(function (shard) {
                    const shardStart = shard.shard_id * manifest.shard_size;
                    (shard.reads || []).forEach(function (record, offset) {
                        recordsById[shardStart + offset] = record;
                    });
                });
                const records = readIds.map(function (readId) {
                    return recordsById[readId];
                }).filter(function (record) {
                    return Boolean(record);
                });
                return expandShardedReadPayload(manifest, group, page, records);
            });
        });
}

function loadReadChunk(rowIndex, page) {
    if (!data.chunk_url) return Promise.resolve(null);
    const url = data.chunk_url[rowIndex];
    if (!url) return Promise.resolve(null);
    const chunkKey = data.chunk_key ? data.chunk_key[rowIndex] : "";
    if (url.indexOf("_manifest.json.gz") + 1) {
        return loadShardedReadPage(url, chunkKey, page).catch(function (error) {
            console.error(error);
            return null;
        });
    }
    if (!window.orographerIsoseqChunkCache) window.orographerIsoseqChunkCache = {};
    const pageToken = "__PAGE__";
    const amp = String.fromCharCode(38);
    const hasPageToken = url.indexOf(pageToken) + 1;
    const separator = url.indexOf("?") == -1 ? "?" : amp;
    const pageUrl = hasPageToken
        ? url.replace(pageToken, String(page))
        : url + separator + "page=" + page + amp + "page_size=" + pageState.pageSize;
    if (window.orographerIsoseqChunkCache[pageUrl]) {
        const cachedPayload = window.orographerIsoseqChunkCache[pageUrl];
        if (chunkKey) {
            return Promise.resolve(
                cachedPayload.chunks ? cachedPayload.chunks[chunkKey] : cachedPayload
            );
        }
        return Promise.resolve(cachedPayload);
    }
    return fetch(pageUrl)
        .then(function (response) {
            if (!response.ok) {
                throw new Error("Could not load IsoSeq read chunk: " + pageUrl);
            }
            return parseJsonResponse(response, pageUrl);
        })
        .then(function (payload) {
            payload = expandCompactReadPayload(payload);
            window.orographerIsoseqChunkCache[pageUrl] = payload;
            if (chunkKey) {
                return payload.chunks ? payload.chunks[chunkKey] : payload;
            }
            return payload;
        })
        .catch(function (error) {
            console.error(error);
            return null;
        });
}

function loadCoverageChunk(rowIndex) {
    if (!data.coverage_url) return Promise.resolve(null);
    const coverageUrl = data.coverage_url[rowIndex];
    if (!coverageUrl) return Promise.resolve(null);
    if (!window.orographerIsoseqChunkCache) window.orographerIsoseqChunkCache = {};
    if (window.orographerIsoseqChunkCache[coverageUrl]) {
        return Promise.resolve(window.orographerIsoseqChunkCache[coverageUrl]);
    }
    return fetch(coverageUrl)
        .then(function (response) {
            if (!response.ok) {
                throw new Error("Could not load IsoSeq coverage chunk: " + coverageUrl);
            }
            return parseJsonResponse(response, coverageUrl);
        })
        .then(function (payload) {
            window.orographerIsoseqChunkCache[coverageUrl] = payload;
            return payload;
        })
        .catch(function (error) {
            console.error(error);
            return null;
        });
}

function renderReadPage(rowIndex, page) {
    pageState.requestSerial += 1;
    const requestSerial = pageState.requestSerial;
    pageState.requestedPage = page;
    const pagePromise = loadReadChunk(rowIndex, page);
    const coveragePromise = loadCoverageChunk(rowIndex);
    Promise.all([pagePromise, coveragePromise]).then(function (loaded) {
        if (requestSerial != pageState.requestSerial) return;
        const payload = loaded[0];
        const coveragePayload = loaded[1];
        if (payload) {
            const coverageData =
                payload.coverage_data ||
                (coveragePayload ? coveragePayload.coverage_data : null);
            pageState.rowIndex = rowIndex;
            pageState.page = page;
            pageState.pageSize = payload.page_size || DEFAULT_PAGE_SIZE;
            pageState.assignedReadCount = payload.assigned_read_count || 0;
            replaceColumns(arrow_source, payload.arrow_data);
            replaceColumns(label_source, payload.label_data);
            replaceColumns(coverage_source, coverageData);
            updateCoverageRange(coverageData);
            applyReadVisibility();
            updatePageControls(payload);
            if (coverage_figure) coverage_figure.visible = true;
        }
    });
}

if (typeof page_target != "undefined") {
    if (pageState.rowIndex == null) return;
    const pageSize = pageState.pageSize || DEFAULT_PAGE_SIZE;
    const pageCount = Math.ceil(pageState.assignedReadCount / pageSize);
    const targetPage = page_target == "last" ? Math.max(0, pageCount - 1) : 0;
    renderReadPage(pageState.rowIndex, targetPage);
} else if (typeof page_delta != "undefined") {
    if (pageState.rowIndex == null) return;
    const pageSize = pageState.pageSize || DEFAULT_PAGE_SIZE;
    const pageCount = Math.ceil(pageState.assignedReadCount / pageSize);
    const lastPage = Math.max(0, pageCount - 1);
    const basePage =
        typeof pageState.requestedPage == "number" ? pageState.requestedPage : pageState.page;
    const nextPage = Math.min(lastPage, Math.max(0, basePage + page_delta));
    if (nextPage == pageState.requestedPage) return;
    renderReadPage(pageState.rowIndex, nextPage);
} else {
    if (!indices.length) {
        restoreTranscriptOnlyView();
        return;
    }

    const idx = indices[indices.length - 1];
    const selectedTranscriptId = data.transcript_id[idx] || "";
    setTranscriptVisibility(selectedTranscriptId);
    replaceColumns(arrow_source, null);
    replaceColumns(label_source, null);
    replaceColumns(coverage_source, null);
    updateCoverageRange(null);
    updatePageControls(null);
    if (coverage_figure) coverage_figure.visible = false;

    if (plot_figure) {
        if (data.selected_view_height) {
            plot_figure.y_range.start = data.selected_view_height[idx];
            plot_figure.height = data.selected_plot_height[idx];
        }
    }

    renderReadPage(idx, 0);
}
}
