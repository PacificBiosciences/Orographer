#!/usr/bin/env python

import argparse
import logging
import re
import sys

from .__init__ import __version__
from .bam_parser import validate_bam_file
from .deploy import run_deploy
from .orographer import orographer
from .utils import ALLOWED_REGION_TYPES, OutputConfig, parse_coordinate

logger = logging.getLogger(__name__)


def setup_args():
    """Set up argument parser with plot and deploy as subcommands."""
    parser = argparse.ArgumentParser(
        prog="orographer",
        description=(
            "Parse BAM file and genomic coordinates, then log the information "
            "and analyze split alignments."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-v",
        "--version",
        help=f"Installed version ({__version__})",
        action="version",
        version="%(prog)s " + str(__version__),
    )

    # Create subparsers for commands
    subparsers = parser.add_subparsers(
        dest="command", help="Command to run", metavar="COMMAND", required=True
    )

    # Plot command
    plot_parser = subparsers.add_parser(
        "plot",
        help="Generate alignment plot",
        description=("Parse BAM file and genomic coordinates, then generate alignment plot."),
        epilog="""
Examples:
  orographer plot --bam input.bam --coord chr1:1000-2000 --region-type \\
      complex_sv --ref ref.fa --outdir ./output
  orographer plot --bam /path/to/file.bam --coord chrX:50000-60000 \\
      --region-type paraphase --ref ref.fa --outdir ./output
        """,
    )

    plot_parser.add_argument(
        "--bam",
        help="Path to the input BAM file",
        required=True,
    )

    plot_parser.add_argument(
        "--coord",
        help=(
            "Genomic coordinate(s) in format chrom:start-end (e.g., chr1:1000-2000). "
            "Can be specified multiple times for multiple regions."
        ),
        required=True,
        action="append",
    )

    plot_parser.add_argument(
        "--region-type",
        help="Region type to plot",
        choices=ALLOWED_REGION_TYPES,
        required=True,
    )

    plot_parser.add_argument(
        "--ref",
        help="Path to reference FASTA file (required for mismatch visualization)",
        required=True,
    )

    plot_parser.add_argument(
        "--outdir",
        help="Directory path for output HTML and JSON files",
        required=True,
        type=str,
    )

    plot_parser.add_argument(
        "--prefix",
        help=("Alphanumeric prefix (may include underscores) to prepend to output filenames"),
        type=str,
        default=None,
    )

    plot_parser.add_argument(
        "--gtf",
        help=(
            "Path to bgzip compressed GTF/GFF3 annotation file with .tbi index "
            "for gene track visualization (optional). Create index with: "
            "tabix -p gff file.gtf.gz"
        ),
        required=False,
        default=None,
    )

    plot_parser.add_argument(
        "--vcf",
        help=(
            "Path to VCF file (gzipped or uncompressed) for variant visualization "
            "(optional). If .tbi index exists, uses tabix for efficient region access."
        ),
        required=False,
        default=None,
    )

    plot_parser.add_argument(
        "--other-bam",
        help=(
            "Path to an additional BAM file. May be specified up to two times. "
            "Order in plot: first other-bam at top, second below it, primary "
            "--bam at bottom."
        ),
        action="append",
        default=[],
    )

    plot_parser.add_argument(
        "--other-vcf",
        help=(
            "Path to VCF for the corresponding --other-bam (order matches: "
            "first --other-vcf for first --other-bam). May be specified up to "
            "two times."
        ),
        action="append",
        default=[],
    )

    plot_parser.add_argument(
        "--sample-label",
        help="Display label for the primary BAM (--bam).",
        type=str,
        default=None,
    )
    plot_parser.add_argument(
        "--other-sample-label",
        help=(
            "Display label for an --other-bam. Specify once per --other-bam "
            "(order matches). May be specified up to two times."
        ),
        action="append",
        default=[],
    )

    plot_parser.add_argument(
        "--verbose",
        help="Write verbose output to stderr.",
        required=False,
        default=False,
        action="store_true",
    )

    # Deploy command
    deploy_parser = subparsers.add_parser(
        "deploy",
        help="Start HTTP server to serve generated plots",
        description=("Start a simple HTTP server to serve Bokeh plots with external JSON files."),
    )

    deploy_parser.add_argument(
        "--outdir",
        help="Directory path containing HTML and JSON files to serve",
        required=True,
        type=str,
    )

    deploy_parser.add_argument(
        "--port",
        help="Port number to serve on (default: 8000)",
        type=int,
        default=8000,
    )

    # Parse arguments
    args = parser.parse_args()

    # Validate prefix and multi-BAM args (for plot command)
    if args.command == "plot":
        if args.prefix is not None and not re.match(r"^[a-zA-Z0-9_]+$", args.prefix):
            parser.error(
                "--prefix must be alphanumeric with underscores "
                "(letters, numbers, and underscores only)"
            )
        other_bam_list = args.other_bam if args.other_bam else []
        other_vcf_list = args.other_vcf if args.other_vcf else []
        other_sample_label_list = args.other_sample_label if args.other_sample_label else []
        if len(other_bam_list) > 2:
            parser.error("--other-bam may be specified at most two times.")
        if len(other_vcf_list) > 2:
            parser.error("--other-vcf may be specified at most two times.")
        if len(other_sample_label_list) > 2:
            parser.error("--other-sample-label may be specified at most two times.")
        if other_sample_label_list and len(other_sample_label_list) != len(other_bam_list):
            parser.error(
                f"--other-sample-label must be specified once per --other-bam "
                f"(expected {len(other_bam_list)} labels for {len(other_bam_list)} "
                f"other BAMs, got {len(other_sample_label_list)})."
            )

    return args


def setup_logging(verbose: bool):
    """Set up logging configuration to log to standard error."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )
    # Silence NumExpr INFO messages about thread defaults
    logging.getLogger("numexpr").setLevel(logging.WARNING)
    # Silence Bokeh INFO messages about file overwriting
    logging.getLogger("bokeh.io.state").setLevel(logging.WARNING)
    return logging.getLogger(__name__)


def run_plot_command(args):
    """Run the plot command."""
    # Set up logging
    logger = setup_logging(args.verbose)

    # Build BAM and VCF lists (validation already done in setup_args)
    other_bam_list = args.other_bam if args.other_bam else []
    other_vcf_list = args.other_vcf if args.other_vcf else []
    other_sample_label_list = args.other_sample_label if args.other_sample_label else []

    bam_files = [args.bam, *other_bam_list]

    # Validate all BAM files
    for bam_path in bam_files:
        logger.debug(f"Validating BAM file: {bam_path}")
        validate_bam_file(bam_path)
    logger.debug("BAM file validation successful")

    # Parse coordinates
    coords = args.coord if isinstance(args.coord, list) else [args.coord]
    coords = [coord.replace("_", "") for coord in coords]
    for coord in coords:
        chromosome, start, end = parse_coordinate(coord)
        logger.debug(
            f"Coordinate parsing successful: {coord} -> chromosome={chromosome}, "
            f"start={start}, end={end}"
        )

    output_config = OutputConfig(args.outdir, args.prefix)
    orographer(
        args.region_type,
        args.bam,
        coords,
        args.ref,
        output_config,
        args.gtf,
        args.vcf,
        other_bam_list,
        other_vcf_list,
        args.sample_label,
        other_sample_label_list,
    )


def run_deploy_command(args):
    """Run the deploy command."""
    run_deploy(args.outdir, args.port)


def main():
    print(f"\nOrographer v{__version__}", file=sys.stderr)
    args = setup_args()

    if args.command == "plot":
        run_plot_command(args)
    elif args.command == "deploy":
        run_deploy_command(args)
    else:
        # Should not happen since subparsers are required
        logger.error(f"Unknown command: {args.command}")
        sys.exit(1)


if __name__ == "__main__":
    logger.error(
        "You are running this module directly, which should only be done for debugging",
    )
    sys.exit(main() or 0)
