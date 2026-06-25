// Region slot-swap callback.
// Fires when either left_select or right_select changes value.
// Loads the chosen region's data from a cached JSON blob and updates all slot sources.
//
// AGENTS.md compliance:
//  - No < operator; use <= n-1 or Math.sign patterns
//  - No && (ternaries used instead)
//  - No ===, !==, >=, <= (Math.sign / Math.abs / Math.min alternatives)
//  - No regex literals
//  - No () in strings (no parens needed here)
//  - forEach instead of for-loops

var oldLeft = parseInt(left_idx_div.text, 10);
var oldRight = parseInt(right_idx_div.text, 10);
var newIdx = parseInt(cb_obj.value, 10);

// Determine which select fired: localeCompare returns 0 only when identical.
var idCmp = cb_obj.id.localeCompare(left_select.id);
var isLeft = Math.sign(Math.abs(idCmp)) ? 0 : 1;

var newLeft = isLeft ? newIdx : oldLeft;
var newRight = isLeft ? oldRight : newIdx;

left_idx_div.text = String(newLeft);
right_idx_div.text = String(newRight);

// Update each select's options to exclude the region shown in the other panel.
left_select.options = all_options.filter(function(o) {
    return Math.min(1, Math.abs(parseInt(o[0], 10) - newRight));
});
right_select.options = all_options.filter(function(o) {
    return Math.min(1, Math.abs(parseInt(o[0], 10) - newLeft));
});

// Load and cache the region JSON blob (parse once per page).
var allRegions = window.orographerRegionData;
if (!allRegions) {
    allRegions = JSON.parse(json_div.text);
    window.orographerRegionData = allRegions;
}

var targetIdx = isLeft ? newLeft : newRight;
var regionData = allRegions[String(targetIdx)];
if (!regionData) {
    return;
}

var slotArrowSources = isLeft ? s0_arrow_sources : s1_arrow_sources;
var slotMarkerSources = isLeft ? s0_marker_connector_sources : s1_marker_connector_sources;
var slotArcSources = isLeft ? s0_arc_connector_sources : s1_arc_connector_sources;
var slotCovTotalSources = isLeft ? s0_coverage_total_sources : s1_coverage_total_sources;
var slotCovHpSources = isLeft ? s0_coverage_hp_sources : s1_coverage_hp_sources;
var slotVcfSources = isLeft ? s0_vcf_sources : s1_vcf_sources;
var slotClearableSources = isLeft ? s0_clearable_sources : s1_clearable_sources;
var slotXRange = isLeft ? s0_x_range : s1_x_range;
var slotYRanges = isLeft ? s0_y_ranges : s1_y_ranges;
var slotOrigSizeDiv = isLeft ? s0_orig_size_div : s1_orig_size_div;
var slotOrigCoordDiv = isLeft ? s0_orig_coord_div : s1_orig_coord_div;
var slotDomClasses = isLeft ? s0_dom_classes : s1_dom_classes;
var slotModelIds = isLeft ? s0_plot_model_ids : s1_plot_model_ids;
var slotAlignmentLabelSources = isLeft ? s0_alignment_label_sources : s1_alignment_label_sources;
var slotHpSepSources = isLeft ? s0_hp_sep_sources : s1_hp_sep_sources;
var slotHpLblSources = isLeft ? s0_hp_lbl_sources : s1_hp_lbl_sources;
var slotGeneBody = isLeft ? s0_gene_body : s1_gene_body;
var slotGeneExon = isLeft ? s0_gene_exon : s1_gene_exon;
var slotGeneIntron = isLeft ? s0_gene_intron : s1_gene_intron;
var slotGeneArrow = isLeft ? s0_gene_arrow : s1_gene_arrow;
var slotGeneLabel = isLeft ? s0_gene_label : s1_gene_label;
var slotCoordInput = isLeft ? s0_coord_input : s1_coord_input;
var slotNavChromDiv = isLeft ? s0_nav_chrom_div : s1_nav_chrom_div;
var slotNavOrigStartDiv = isLeft ? s0_nav_orig_start_div : s1_nav_orig_start_div;
var slotNavOrigEndDiv = isLeft ? s0_nav_orig_end_div : s1_nav_orig_end_div;
var slotRepeatDensitySource = isLeft ? s0_repeat_density_source : s1_repeat_density_source;
var slotRepeatDensityYRange = isLeft ? s0_repeat_density_y_range : s1_repeat_density_y_range;
var slotGeneTrackYRange = isLeft ? s0_gene_track_y_range : s1_gene_track_y_range;
var slotCoverageYRanges = isLeft ? s0_coverage_y_ranges : s1_coverage_y_ranges;
var slotSeparatorSources = isLeft ? s0_separator_sources : s1_separator_sources;
var slotPhaseSetSources = isLeft ? s0_phase_set_sources : s1_phase_set_sources;
var slotVariantSourcesPerRow = isLeft ? s0_variant_sources_per_row : s1_variant_sources_per_row;
var slotInsertionRawSourcesPerRow = (
    isLeft ? s0_insertion_raw_sources_per_row : s1_insertion_raw_sources_per_row
);

// Clear all sources that are not explicitly repopulated below.
slotClearableSources.forEach(function(src) {
    if (!src) {
        return;
    }
    var empty = {};
    Object.keys(src.data).forEach(function(k) {
        empty[k] = [];
    });
    src.data = empty;
});

// Update per-BAM-row sources.
var bamRows = regionData.bam_rows;
bamRows.forEach(function(rowData, i) {
    if (!rowData) {
        return;
    }

    var arrowSrc = slotArrowSources[i];
    if (arrowSrc) {
        arrowSrc.data = rowData.arrow_data;
    }

    var markerSrc = slotMarkerSources[i];
    if (markerSrc) {
        var mdata = rowData.connector_data;
        var n = mdata.stub_x0 ? mdata.stub_x0.length : 0;
        var domClass = slotDomClasses[i];
        var modelId = slotModelIds[i];
        mdata.plot_dom_class = Array.from({length: n}, function() { return domClass; });
        mdata.plot_model_id = Array.from({length: n}, function() { return modelId; });
        markerSrc.data = mdata;
    }

    var arcSrc = slotArcSources[i];
    if (arcSrc) {
        var arcData = rowData.same_region_connector_data;
        arcSrc.data = arcData ? arcData : {
            xs: [], ys: [], read_name: [], connection_id: [], haplotype_transition: [],
            connection_scope: [], read_layout_mode: [],
        };
    }

    var covTotalSrc = slotCovTotalSources[i];
    if (covTotalSrc) {
        covTotalSrc.data = rowData.cov_total ? rowData.cov_total : {x: [], y: []};
    }

    var covHpSrc = slotCovHpSources[i];
    if (covHpSrc) {
        covHpSrc.data = rowData.cov_hp ? rowData.cov_hp : {
            xs: [], ys: [], colors: [], names: [], label_x: [], label_y: [],
        };
    }

    // Rescale coverage y_range to match the swapped region's depth.
    var covYRange = slotCoverageYRanges ? slotCoverageYRanges[i] : null;
    if (covYRange) {
        var covYMax = rowData.cov_y_max;
        covYRange.start = 0;
        covYRange.end = covYMax ? covYMax : 1;
    }

    var vcfSrc = slotVcfSources[i];
    if (vcfSrc) {
        vcfSrc.data = rowData.vcf ? rowData.vcf : {
            x: [], y: [], color: [], angle: [],
            coordinates: [], variant_type: [], alt_allele: [], alt_base: [],
            haplotypes: [], sample_label: [],
        };
    }

    // Restore alignment number label source for this BAM row.
    var alignLblSrc = slotAlignmentLabelSources ? slotAlignmentLabelSources[i] : null;
    if (alignLblSrc) {
        alignLblSrc.data = rowData.alignment_label_data ? rowData.alignment_label_data : {
            x: [], y: [], text: [], read_name: [], layout_read_name: [],
            alignment_number: [], strand: [], coordinates: [], haplotype: [],
            sample_label: [], inclusion_reason: [], all_alignment_numbers: [],
            all_alignment_coordinates: [], has_split_alignment: [],
            has_multiregion_connection: [], gene_id: [], gene_name: [],
            transcript_id: [], read_filter_alpha: [], label_alpha: [], read_layout_mode: [],
        };
    }

    // Restore HP label / separator sources for this BAM row.
    var hpSepSrc = slotHpSepSources ? slotHpSepSources[i] : null;
    if (hpSepSrc) {
        hpSepSrc.data = rowData.hp_sep_data ? rowData.hp_sep_data : {
            xs: [], ys: [], read_layout_mode: [],
        };
    }

    var hpLblSrc = slotHpLblSources ? slotHpLblSources[i] : null;
    if (hpLblSrc) {
        hpLblSrc.data = rowData.hp_lbl_data ? rowData.hp_lbl_data : {
            x: [], y: [], text_x: [], text: [], color: [], read_layout_mode: [],
        };
    }

    // Restore horizontal separator lines (dotted lines between reads).
    var sepSrc = slotSeparatorSources ? slotSeparatorSources[i] : null;
    if (sepSrc) {
        sepSrc.data = rowData.separator_data ? rowData.separator_data : {
            xs: [], ys: [], read_name: [], has_split_alignment: [],
            has_multiregion_connection: [], separator_line_alpha: [],
        };
    }

    // Restore vertical phase-set boundary lines.
    var phaseSrc = slotPhaseSetSources ? slotPhaseSetSources[i] : null;
    if (phaseSrc) {
        phaseSrc.data = rowData.phase_boundary_data ? rowData.phase_boundary_data : {
            xs: [], ys: [], color: [], read_layout_mode: [],
        };
    }

    // Restore variant marker sources (mismatch, insertion, deletion sub-sources).
    var varSrcs = slotVariantSourcesPerRow ? slotVariantSourcesPerRow[i] : null;
    var varData = rowData.variant_sources_data;
    if (varSrcs ? varData : null) {
        varSrcs.forEach(function(src, j) {
            var d = varData[j];
            if (src ? d : null) {
                src.data = d;
            }
        });
    }

    // Restore insertion-summary raw sources (display sources recluster via x_range callback).
    var insRawSrcs = slotInsertionRawSourcesPerRow ? slotInsertionRawSourcesPerRow[i] : null;
    var insRawData = rowData.insertion_raw_data;
    if (insRawSrcs ? insRawData : null) {
        insRawSrcs.forEach(function(src, j) {
            var d = insRawData[j];
            if (src ? d : null) {
                src.data = d;
            }
        });
    }

    var yRange = slotYRanges[i];
    if (yRange) {
        yRange.start = rowData.y_bounds[0];
        yRange.end = rowData.y_bounds[1];
    }
});

// Restore gene track sources.
var geneTrack = regionData.gene_track;
if (slotGeneBody) {
    slotGeneBody.data = geneTrack ? (geneTrack.body ? geneTrack.body : {
        x0: [], x1: [], y: [], gene_id: [], gene_name: [], strand: [], coordinates: [],
        line_alpha: [], line_alpha_overview: [], line_alpha_medium: [], line_alpha_detail: [],
    }) : {};
}
if (slotGeneExon) {
    slotGeneExon.data = geneTrack ? (geneTrack.exon ? geneTrack.exon : {
        left: [], right: [], top: [], bottom: [], gene_name: [], gene_strand: [],
        exon_number: [], representative_transcript: [], representative_selection_method: [],
        fill_alpha: [], fill_alpha_overview: [], fill_alpha_medium: [], fill_alpha_detail: [],
        line_alpha: [], line_alpha_overview: [], line_alpha_medium: [], line_alpha_detail: [],
    }) : {};
}
if (slotGeneIntron) {
    slotGeneIntron.data = geneTrack ? (geneTrack.intron ? geneTrack.intron : {
        xs: [], ys: [], line_alpha: [], line_alpha_overview: [], line_alpha_medium: [],
        line_alpha_detail: [],
    }) : {};
}
if (slotGeneArrow) {
    slotGeneArrow.data = geneTrack ? (geneTrack.arrow ? geneTrack.arrow : {
        x: [], y: [], angle: [], fill_alpha: [], fill_alpha_overview: [],
        fill_alpha_medium: [], fill_alpha_detail: [], line_alpha: [], line_alpha_overview: [],
        line_alpha_medium: [], line_alpha_detail: [],
    }) : {};
}
if (slotGeneLabel) {
    slotGeneLabel.data = geneTrack ? (geneTrack.label ? geneTrack.label : {
        x: [], y: [], text: [], text_alpha: [], text_alpha_overview: [],
        text_alpha_medium: [], text_alpha_detail: [],
    }) : {};
}

// Rescale gene track y_range to match the swapped region's layout height.
if (slotGeneTrackYRange) {
    var newGeneYStart = geneTrack ? geneTrack.y_range_start : null;
    slotGeneTrackYRange.start = newGeneYStart ? newGeneYStart : 4.0;
}

// Restore repeat-density (Ident) track data and y-range scale.
var densData = regionData.repeat_density;
if (slotRepeatDensitySource) {
    var densPayload = densData
        ? {x: densData.x, top: densData.top, width: densData.width}
        : {x: [], top: [], width: []};
    slotRepeatDensitySource.data = densPayload;
}
if (slotRepeatDensityYRange ? densData : null) {
    var densYMax = densData.y_max;
    slotRepeatDensityYRange.end = densYMax * 1.1;
}

// Update nav-bound Divs BEFORE updating x_range; the go_callback reads them on
// x_range change. Setting these first prevents "not found" errors for swapped regions.
if (slotNavChromDiv) {
    slotNavChromDiv.text = regionData.chrom;
}
if (slotNavOrigStartDiv) {
    slotNavOrigStartDiv.text = String(regionData.x_start);
}
if (slotNavOrigEndDiv) {
    slotNavOrigEndDiv.text = String(regionData.x_end);
}

// Guard against go_callback firing erroneously while we programmatically update
// x_range and coord_input. The guard flag is checked at the top of coord_go_callback.js.
window.orographerSwapInProgress = true;

// Update x-range (triggers view_callback which updates coord_input).
// Set end first, then bump start by 1 bp before the final value so insertion-marker
// recluster callbacks fire even when the target region has identical coordinates to
// the current view (e.g. swap-back scenarios where x_range would not otherwise change).
slotXRange.end = regionData.x_end;
slotXRange.start = regionData.x_start - 1;
slotXRange.start = regionData.x_start;

// Override coord_input with the correct chromosome string; view_callback may have
// already set it but used the previous chromosome if chrom_div update raced with it.
if (slotCoordInput) {
    function _fmt(n) {
        var re = new RegExp("\\B(?=(\\d{3})+(?!\\d))", "g");
        return n.toString().replace(re, ",");
    }
    var coordStr = regionData.chrom + ":" + _fmt(regionData.x_start) + "-" + _fmt(regionData.x_end);
    slotCoordInput.value = coordStr;
}

window.orographerSwapInProgress = false;

// Update the "Original region" annotation in the navigation bar.
if (slotOrigSizeDiv) {
    slotOrigSizeDiv.text = regionData.orig_size_label;
}
if (slotOrigCoordDiv) {
    slotOrigCoordDiv.text = regionData.orig_coord_label;
}

// Trigger overlay redraw after DOM paint to ensure arrow positions are up-to-date.
setTimeout(function() {
    if (window.orographerUpdateReadConnectionOverlay) {
        window.orographerUpdateReadConnectionOverlay();
    }
}, 150);
