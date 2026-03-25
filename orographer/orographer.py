import logging

from .bam_parser import (
    collect_all_alignments_for_reads,
    fetch_all_alignments,
    validate_bam_file,
)
from .gtf_parser import parse_annotation_file
from .plot_bokeh import plot_reads_bokeh
from .utils import (
    COMPLEX_SV_REGION_TYPE,
    PARAPHASE_REGION_TYPE,
    ProcessingPaths,
    Region,
    parse_coordinate,
    validate_output_dir,
)
from .vcf_parser import parse_vcf_file

logger = logging.getLogger(__name__)


def log_final_information(bam_file, region: Region):
    logger.debug("=" * 50)
    logger.debug("BAM COORDINATE PARSER RESULTS")
    logger.debug("=" * 50)
    logger.debug(f"BAM file: {bam_file}")
    logger.debug(f"Chromosome: {region.chromosome}")
    logger.debug(f"Start position: {region.start}")
    logger.debug(f"End position: {region.end}")
    logger.debug(f"Coordinate range: {region.end - region.start} bp")
    logger.debug("=" * 50)


def log_segments(segments_by_read, message_prefix=""):
    """Log segments per read (dict is haplotype-keyed and ordered from BAM parser)."""
    if segments_by_read:
        logger.debug(f"{message_prefix}Reads found:")
        for i, read_name in enumerate(sorted(segments_by_read.keys()), 1):
            logger.debug(f"  {i}. Read: {read_name}")
        logger.debug("=" * 50)
        logger.debug(f"{message_prefix}COLLECTING ALL SEGMENTS")
        logger.debug("=" * 50)
        for read_name, segments in segments_by_read.items():
            logger.debug(f"Read: {read_name}")
            for idx, segment in enumerate(segments, 1):
                strand = "+" if segment.is_fwd_strand else "-"
                chrom = segment.chrom
                ref_start = segment.pos
                ref_end = segment.end
                q_start = segment.fwd_read_start
                q_end = segment.fwd_read_end
                logger.debug(
                    f"  {idx}. q_start={q_start}, q_end={q_end}, strand={strand}, "
                    f"ref={chrom}:{ref_start}-{ref_end}"
                )


def process_complex_sv(paths: ProcessingPaths, region: Region):
    """Analyze split alignments in the region; return segment and annotation data."""
    log_final_information(paths.bam_file, region)

    gene_annotations = parse_annotation_file(paths.gtf_file, region)
    vcf_variants = parse_vcf_file(paths.vcf_file, region)

    split_segments_by_read = fetch_all_alignments(
        paths.bam_file,
        region,
        only_split=True,
        reference_path=paths.reference_path,
    )
    logger.debug(
        f"Found {len(split_segments_by_read)} reads with SA tags in region "
        f"{region.chromosome}:{region.start}-{region.end}"
    )

    if split_segments_by_read:
        all_segments_by_read = collect_all_alignments_for_reads(split_segments_by_read)
        log_segments(all_segments_by_read, message_prefix="Split ")
        return {
            "segments_by_read": all_segments_by_read,
            "chromosome": region.chromosome,
            "start": region.start,
            "end": region.end,
            "coordinate_str": region.coordinate_str,
            "gene_annotations": gene_annotations,
            "vcf_variants": vcf_variants,
        }
    logger.debug("No split reads found in the specified region.")
    return {
        "segments_by_read": {},
        "chromosome": region.chromosome,
        "start": region.start,
        "end": region.end,
        "coordinate_str": region.coordinate_str,
        "gene_annotations": gene_annotations,
        "vcf_variants": vcf_variants,
    }


def process_paraphase(paths: ProcessingPaths, region: Region):
    """Fetch all alignments in the region; return segment and annotation data."""
    log_final_information(paths.bam_file, region)

    gene_annotations = parse_annotation_file(paths.gtf_file, region)
    vcf_variants = parse_vcf_file(paths.vcf_file, region)

    all_segments_by_read = fetch_all_alignments(
        paths.bam_file,
        region,
        reference_path=paths.reference_path,
    )
    logger.debug(
        f"Found {len(all_segments_by_read)} reads in region "
        f"{region.chromosome}:{region.start}-{region.end}"
    )

    if all_segments_by_read:
        log_segments(all_segments_by_read, message_prefix="")
        return {
            "segments_by_read": all_segments_by_read,
            "chromosome": region.chromosome,
            "start": region.start,
            "end": region.end,
            "coordinate_str": region.coordinate_str,
            "gene_annotations": gene_annotations,
            "vcf_variants": vcf_variants,
        }
    logger.debug("No reads found in the specified region.")
    return {
        "segments_by_read": {},
        "chromosome": region.chromosome,
        "start": region.start,
        "end": region.end,
        "coordinate_str": region.coordinate_str,
        "gene_annotations": gene_annotations,
        "vcf_variants": vcf_variants,
    }


def orographer(
    region_type,
    bam_file,
    coordinate_strs,
    reference_path,
    output_config,
    gtf_file,
    vcf_file,
    other_bam_files,
    other_vcf_files,
    sample_label,
    other_sample_labels,
):
    """Run alignment processing and plot.

    output_config is an OutputConfig(output_dir, prefix). other_bam_files and
    other_vcf_files are lists (0-2 items). sample_label is the display label
    for the primary BAM (optional). other_sample_labels is a list of display
    labels for other BAMs, length 0 or len(other_bam_files).
    """
    bam_files = [bam_file] + (other_bam_files if other_bam_files else [])
    vcf_files = [vcf_file] + (other_vcf_files if other_vcf_files else [])
    while len(vcf_files) < len(bam_files):
        vcf_files.append(None)

    if len(bam_files) > 3:
        raise ValueError("Total number of BAMs cannot exceed 3.")
    if len(vcf_files) > len(bam_files):
        raise ValueError("Number of VCFs cannot exceed number of BAMs.")
    other_labels = other_sample_labels if other_sample_labels else []
    if other_labels and len(other_labels) != len(other_bam_files or []):
        raise ValueError(
            f"other_sample_labels must have length 0 or {len(other_bam_files)} "
            f"(number of other BAMs), got {len(other_labels)}."
        )

    for bam_path in bam_files:
        logger.debug(f"Validating BAM file: {bam_path}")
        validate_bam_file(bam_path)
    logger.debug("BAM file validation successful")

    logger.debug(f"Validating output directory: {output_config.output_dir}")
    validate_output_dir(output_config.output_dir)
    logger.debug("Output directory validation successful")

    if isinstance(coordinate_strs, str):
        coordinate_strs = [coordinate_strs]

    region_data_list = []
    for coordinate_str in coordinate_strs:
        logger.debug(f"Processing coordinate: {coordinate_str}")
        chromosome, start, end = parse_coordinate(coordinate_str)

        region = Region(chromosome, start, end, coordinate_str)
        rows_in_bam_order = []
        for bam_index in range(len(bam_files)):
            current_bam = bam_files[bam_index]
            current_vcf = vcf_files[bam_index] if bam_index < len(vcf_files) else None
            paths = ProcessingPaths(current_bam, reference_path, gtf_file, current_vcf)

            if region_type == COMPLEX_SV_REGION_TYPE:
                row_data = process_complex_sv(paths, region)
            elif region_type == PARAPHASE_REGION_TYPE:
                row_data = process_paraphase(paths, region)
            else:
                raise ValueError(f"Invalid region type: {region_type}")

            row_data["region_type"] = region_type
            if bam_index == 0:
                row_data["sample_label"] = sample_label
            elif bam_index <= len(other_labels):
                row_data["sample_label"] = other_labels[bam_index - 1]
            else:
                row_data["sample_label"] = None
            rows_in_bam_order.append(row_data)

        # Display order: first other (top), second other (middle), primary (bottom)
        if len(rows_in_bam_order) == 3:
            bam_rows_display_order = [
                rows_in_bam_order[1],
                rows_in_bam_order[2],
                rows_in_bam_order[0],
            ]
        elif len(rows_in_bam_order) == 2:
            bam_rows_display_order = [rows_in_bam_order[1], rows_in_bam_order[0]]
        else:
            bam_rows_display_order = [rows_in_bam_order[0]]

        # Use first row for shared region metadata
        first_row = rows_in_bam_order[0]
        region_data_list.append(
            {
                "region": region,
                "gene_annotations": first_row["gene_annotations"],
                "bam_rows": bam_rows_display_order,
            }
        )

    if region_data_list:
        logger.debug("=" * 50)
        logger.debug("CREATING ALIGNMENT PLOT")
        logger.debug("=" * 50)
        bokeh_plot_file = plot_reads_bokeh(region_data_list, output_config)
        if bokeh_plot_file:
            logger.debug(f"Alignment plot saved to: {bokeh_plot_file}")
            logger.debug(
                f"Use `orographer deploy --outdir {output_config.output_dir}` to serve the results"
            )
        else:
            logger.warning("Failed to create alignment plot")
    else:
        logger.warning("No regions with data to plot.")
