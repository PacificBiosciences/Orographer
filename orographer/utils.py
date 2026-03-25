import logging
import os
import re
from pathlib import Path
from typing import NamedTuple

logger = logging.getLogger(__name__)

COMPLEX_SV_REGION_TYPE = "complex_sv"
PARAPHASE_REGION_TYPE = "paraphase"

ALLOWED_REGION_TYPES = [COMPLEX_SV_REGION_TYPE, PARAPHASE_REGION_TYPE]


class Region(NamedTuple):
    """Genomic region: chromosome, start/end, and original coordinate string."""

    chromosome: str
    start: int
    end: int
    coordinate_str: str


class ProcessingPaths(NamedTuple):
    """File paths for processing one BAM: BAM, reference, GTF, and optional VCF."""

    bam_file: str
    reference_path: str
    gtf_file: str
    vcf_file: str


class OutputConfig(NamedTuple):
    """Output directory and optional filename prefix for plots."""

    output_dir: str
    prefix: str | None


def get_tabix_chromosome_name(tabix_file, target_chrom: str) -> str | None:
    """
    Get the exact chromosome name as it appears in the tabix index.
    Tries common variations (chr1 vs 1, etc.)

    Returns:
        The chromosome name as it appears in the index, or None if not found
    """
    # Get all contigs from the index
    contigs = tabix_file.contigs

    # Try exact match first
    if target_chrom in contigs:
        return target_chrom

    # Try without 'chr' prefix
    chrom_no_chr = target_chrom.replace("chr", "")
    if chrom_no_chr in contigs:
        return chrom_no_chr

    # Try with 'chr' prefix
    chrom_with_chr = f"chr{target_chrom}" if not target_chrom.startswith("chr") else target_chrom
    if chrom_with_chr in contigs:
        return chrom_with_chr

    # Try other variations
    for contig in contigs:
        if contig.replace("chr", "") == target_chrom.replace("chr", ""):
            return contig

    return None


def chromosome_matches(seq_name: str, target_chrom: str) -> bool:
    """Check if sequence name matches target chrom, handling chr prefix."""
    if seq_name == target_chrom:
        return True
    if seq_name == target_chrom.replace("chr", ""):
        return True
    return f"chr{seq_name}" == target_chrom


def validate_output_dir(output_dir):
    """
    Validate that the output directory path is valid and can be created.

    Args:
        output_dir (str): Path to the output directory

    Raises:
        ValueError: If output directory path is invalid or cannot be created
    """
    if not output_dir:
        raise ValueError("Output directory path cannot be empty")

    # Convert to Path object for easier manipulation
    output_path = Path(output_dir)

    # If the directory already exists, check that it's actually a directory
    if output_path.exists():
        if not output_path.is_dir():
            raise ValueError(f"Output path exists but is not a directory: {output_dir}")
        # Check if it's writable
        if not os.access(output_path, os.W_OK):
            raise ValueError(f"Output directory is not writable: {output_dir}")
    else:
        # Directory doesn't exist, try to create it
        try:
            output_path.mkdir(parents=True, exist_ok=True)
        except (OSError, PermissionError) as err:
            raise ValueError(f"Cannot create output directory {output_dir}: {err}") from err


def parse_coordinate(coord_string):
    """
    Parse genomic coordinate string in format 'chrom:start-end'.

    Args:
        coord_string (str): Coordinate string in format 'chrom:start-end'

    Returns:
        tuple: (chromosome, start, end) where start and end are integers

    Raises:
        ValueError: If coordinate format is invalid
    """
    # Pattern to match chrom:start-end format
    pattern = r"^([^:]+):(\d+)-(\d+)$"
    match = re.match(pattern, coord_string)

    if not match:
        raise ValueError(
            f"Invalid coordinate format: {coord_string}. Expected format: chrom:start-end"
        )

    chromosome = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))

    # Validate that start < end
    if start >= end:
        raise ValueError(f"Start position ({start}) must be less than end position ({end})")

    return chromosome, start, end
