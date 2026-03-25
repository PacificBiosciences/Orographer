import gzip
import logging
import os
from dataclasses import dataclass, field

import pysam

from .utils import Region, chromosome_matches, get_tabix_chromosome_name

logger = logging.getLogger(__name__)


@dataclass
class VCFVariant:
    """Represents a variant from a VCF file."""

    chrom: str
    pos: int  # 1-based position
    ref: str  # Reference allele
    alt: str  # Alternate allele
    variant_type: str  # 'SNP', 'INSERTION', 'DELETION'
    alt_base: str | None = None  # For SNPs, the alt base (A, C, G, T)
    haplotypes: list[str] = field(
        default_factory=list
    )  # List of haplotype names with non-ref genotype

    def overlaps_region(self, region_start: int, region_end: int) -> bool:
        """Check if this variant overlaps with the given region."""
        # For deletions, end position is pos + len(ref) - 1
        variant_end = self.pos + len(self.ref) - 1 if self.variant_type == "DELETION" else self.pos
        return self.pos <= region_end and variant_end >= region_start


def _parse_vcf_header(file_path: str) -> list[str]:
    """
    Parse VCF header to extract sample names from the #CHROM line.

    Args:
        file_path: Path to VCF file (can be gzipped or uncompressed)

    Returns:
        List of sample names, empty list if not found
    """
    sample_names = []
    open_func = gzip.open if file_path.endswith(".gz") else open
    mode = "rt" if file_path.endswith(".gz") else "r"

    with open_func(file_path, mode) as f:
        for line in f:
            if line.startswith("#CHROM"):
                # Format: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1 ...
                parts = line.strip().split("\t")
                if len(parts) > 9:
                    sample_names = parts[9:]  # Everything after FORMAT column
                break

    return sample_names


def parse_vcf_file(file_path: str, region: Region) -> list[VCFVariant]:
    """
    Parse VCF file and extract variants overlapping the region.
    Uses tabix if .tbi index is present, otherwise reads the entire file.

    Args:
        file_path: Path to VCF file (can be gzipped or uncompressed)
        region: Genomic region (chromosome, start, end, coordinate_str)

    Returns:
        List of VCFVariant objects for variants overlapping the region
    """
    if not file_path:
        return []

    variants = []
    tbi_path = file_path + ".tbi"
    use_tabix = os.path.exists(tbi_path)

    # Parse header to get sample names
    sample_names = _parse_vcf_header(file_path)

    logger.debug(f"Parsing VCF: {file_path} (tabix={use_tabix}, samples={len(sample_names)})")

    try:
        if use_tabix:
            # Use tabix for indexed VCF
            # Tabix doesn't expose header easily; we parse it separately
            # sample_names should already be parsed from _parse_vcf_header above

            with pysam.TabixFile(file_path) as tabix_file:
                # Get the exact chromosome name as it appears in the index
                tabix_chrom = get_tabix_chromosome_name(tabix_file, region.chromosome)

                if tabix_chrom is None:
                    logger.warning(
                        f"Chromosome {region.chromosome} not in VCF tabix. "
                        f"Contigs: {list(tabix_file.contigs)[:10]}"
                    )
                    return []

                # Tabix uses 0-based start, exclusive end
                tabix_start = region.start - 1
                tabix_end = region.end + 1
                records = tabix_file.fetch(tabix_chrom, tabix_start, tabix_end)
                for line in records:
                    if isinstance(line, bytes):
                        line = line.decode("utf-8")
                    # Tabix already filtered by chrom; pass tabix_chrom
                    # and skip chromosome check in _parse_vcf_line
                    variant = _parse_vcf_line(
                        line,
                        tabix_chrom,
                        region.start,
                        region.end,
                        sample_names,
                        skip_chrom_check=True,
                    )
                    if variant:
                        variants.append(variant)
        else:
            # Read entire file (no index)
            open_func = gzip.open if file_path.endswith(".gz") else open
            mode = "rt" if file_path.endswith(".gz") else "r"

            with open_func(file_path, mode) as f:
                for line in f:
                    if line.startswith("#"):
                        continue

                    variant = _parse_vcf_line(
                        line, region.chromosome, region.start, region.end, sample_names
                    )
                    if variant:
                        variants.append(variant)

        logger.debug(
            f"Found {len(variants)} VCF variants in {region.chromosome}:{region.start}-{region.end}"
        )
        if variants:
            logger.debug(
                f"First few: {[(v.pos, v.ref, v.alt, v.variant_type) for v in variants[:5]]}"
            )
        return variants

    except FileNotFoundError:
        logger.error(f"VCF file not found: {file_path}")
        return []


def _parse_genotype_haplotypes(gt_field: str, sample_name: str) -> list[str]:
    """
    Parse GT (genotype) field and return list of sample names with non-ref alleles.

    Args:
        gt_field: GT field value (e.g., "0|1", "1/0", "1|1", "./.")
        sample_name: Name of the sample

    Returns:
        List containing sample name if it has a non-ref genotype, empty list otherwise
    """
    haplotypes = []

    if not gt_field or gt_field == "." or gt_field == "./." or gt_field == ".|.":
        return haplotypes

    # Check if phased (|) or unphased (/)
    separator = "|" if "|" in gt_field else "/"

    # Split genotype into alleles
    alleles = gt_field.split(separator)

    # Check if any allele is non-reference (> 0)
    has_non_ref = False
    for allele_str in alleles:
        try:
            allele = int(allele_str)
            if allele > 0:
                has_non_ref = True
                break
        except (ValueError, IndexError):
            continue

    # If sample has non-reference genotype, add sample name
    if has_non_ref:
        haplotypes.append(sample_name)

    return haplotypes


def _parse_vcf_line(
    line: str,
    chrom: str,
    region_start: int,
    region_end: int,
    sample_names: list[str],
    skip_chrom_check: bool = False,
) -> VCFVariant | None:
    """
    Parse a single VCF line and return a VCFVariant if it overlaps the region.

    Args:
        line: VCF line to parse
        chrom: Target chromosome name
        region_start: Start of region (1-based)
        region_end: End of region (1-based, inclusive)
        sample_names: List of sample names from VCF header
        skip_chrom_check: If True, skip chromosome matching (e.g., when using tabix)

    Returns:
        VCFVariant object or None if line doesn't overlap or is invalid
    """
    parts = line.strip().split("\t")
    if len(parts) < 5:
        return None

    seq_name = parts[0]
    pos = int(parts[1])  # 1-based
    ref = parts[3].upper()
    alt_string = parts[4].upper()

    # Handle multiple ALT alleles (comma-separated) - take the first one
    alt = alt_string.split(",")[0] if alt_string else ""

    # Skip if not on target chromosome (unless skip_chrom_check is True)
    if not skip_chrom_check and not chromosome_matches(seq_name, chrom):
        return None

    # Handle symbolic alleles (e.g., <DEL>, <INS>, <DUP>)
    is_symbolic = alt.startswith("<") and alt.endswith(">")

    # Skip if doesn't overlap region
    if pos > region_end:
        return None

    # For deletions (including symbolic), check if deletion end overlaps
    if is_symbolic:
        # For symbolic alleles, try to determine type from the symbol
        alt_upper = alt.upper()
        if "DEL" in alt_upper:
            variant_type = "DELETION"
            alt_base = None
            # Symbolic deletions: unknown length, just check position
            if pos < region_start:
                return None
        elif "INS" in alt_upper:
            variant_type = "INSERTION"
            alt_base = None
            if pos < region_start:
                return None
        else:
            # Unknown symbolic variant, skip for now
            return None
    elif len(ref) > len(alt):
        # Deletion
        deletion_end = pos + len(ref) - 1
        if deletion_end < region_start:
            return None
        variant_type = "DELETION"
        alt_base = None
    elif len(ref) < len(alt):
        # Insertion
        if pos < region_start:
            return None
        variant_type = "INSERTION"
        alt_base = None
    elif len(ref) == 1 and len(alt) == 1:
        # SNP
        if pos < region_start:
            return None
        variant_type = "SNP"
        alt_base = alt
    else:
        # Complex variant (MNP or other), treat as SNP if single base
        if len(ref) == 1:
            if pos < region_start:
                return None
            variant_type = "SNP"
            alt_base = alt[0] if alt else None
        else:
            # Skip complex variants for now
            return None

    # Parse genotype information from sample columns
    haplotypes = []
    if len(parts) > 9 and sample_names:
        format_col = parts[8]  # FORMAT column
        format_fields = format_col.split(":")

        # Find GT field index
        gt_index = None
        for idx, fmt_field in enumerate(format_fields):
            if fmt_field == "GT":
                gt_index = idx
                break

        if gt_index is not None:
            # Parse each sample column
            for sample_idx, sample_name in enumerate(sample_names):
                if 9 + sample_idx < len(parts):
                    sample_data = parts[9 + sample_idx]
                    sample_fields = sample_data.split(":")

                    if gt_index < len(sample_fields):
                        gt_field = sample_fields[gt_index]
                        # Get haplotypes with non-ref genotype for this sample
                        sample_haplotypes = _parse_genotype_haplotypes(gt_field, sample_name)
                        haplotypes.extend(sample_haplotypes)

    return VCFVariant(
        chrom=seq_name,
        pos=pos,
        ref=ref,
        alt=alt,
        variant_type=variant_type,
        alt_base=alt_base,
        haplotypes=haplotypes,
    )
