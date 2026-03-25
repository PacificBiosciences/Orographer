"""GTF/GFF3 annotation file parsing for gene track visualization."""

import logging
import os
import re
from dataclasses import dataclass, field

import pysam

from .utils import Region, get_tabix_chromosome_name

logger = logging.getLogger(__name__)


@dataclass
class GeneAnnotation:
    """Represents a gene/transcript with its exons for visualization."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int  # Gene start (min of all exon starts)
    end: int  # Gene end (max of all exon ends)
    strand: str  # '+' or '-'
    exons: list[tuple[int, int, int]] = field(
        default_factory=list
    )  # List of (start, end, exon_number) tuples

    def overlaps_region(self, region_start: int, region_end: int) -> bool:
        """Check if this gene overlaps with the given region."""
        return self.start < region_end and self.end > region_start


def parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GTF attribute column: ``key "value";`` pairs, then unquoted ``key N;`` numerics.
    """
    attrs = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', attr_string):
        key, value = match.groups()
        attrs[key] = value
    for match in re.finditer(r"(\S+)\s+(\d+)\s*;", attr_string):
        key, value = match.group(1), match.group(2)
        attrs.setdefault(key, value)
    return attrs


def parse_gff3_attributes(attr_string: str) -> dict[str, str]:
    """Parse GFF3 attribute string (key=value; format)."""
    attrs = {}
    # GFF3 format: key=value;key2=value2
    for pair in attr_string.split(";"):
        pair = pair.strip()
        if "=" in pair:
            key, value = pair.split("=", 1)
            attrs[key] = value
    return attrs


def detect_annotation_format(file_path: str) -> str:
    """Detect if file is GTF or GFF3 based on extension and content."""
    lower_path = file_path.lower()
    if ".gtf" in lower_path:
        return "gtf"
    elif ".gff3" in lower_path or ".gff" in lower_path:
        return "gff3"
    # Default to GTF
    return "gtf"


def validate_bgzip_index(file_path: str) -> None:
    """
    Validate that the annotation file is bgzip compressed and has a .tbi index.

    Raises:
        ValueError: If file is not bgzip or missing .tbi index
    """
    tbi_path = file_path + ".tbi"
    if not os.path.exists(tbi_path):
        raise ValueError(
            f"Tabix index (.tbi) not found for {file_path}. "
            f"Expected index file: {tbi_path}. "
            f"Create index with: tabix -p gff {file_path}"
        )

    # Try to open with tabix to verify it's properly indexed
    try:
        with pysam.TabixFile(file_path) as _:
            pass
    except OSError as err:
        raise ValueError(f"Failed to open indexed annotation file {file_path}: {err}") from err


def extract_gtf_gene_info(attrs: dict[str, str]) -> tuple[str, str, bool]:
    """
    Extract gene_id and gene_name from GTF attributes.

    Returns:
        tuple: (gene_id, gene_name, should_skip) where should_skip indicates
               if this record should be skipped
    """
    gene_id = attrs.get("gene_id", "")
    gene_name = attrs.get("gene_name", attrs.get("gene_id", "Unknown"))
    return gene_id, gene_name, False


def extract_gff3_gene_info(
    attrs: dict[str, str], feature_type: str, transcript_to_gene: dict[str, str]
) -> tuple[str, str, bool]:
    """
    Extract gene_id and gene_name from GFF3 attributes.
    Also updates transcript_to_gene mapping for mRNA/transcript features.

    Returns:
        tuple: (gene_id, gene_name, should_skip) where should_skip indicates
               if this record should be skipped
    """
    if feature_type == "gene":
        gene_id = attrs.get("ID", attrs.get("gene_id", ""))
        gene_name = attrs.get("Name", attrs.get("gene_name", gene_id))
        return gene_id, gene_name, False

    elif feature_type in ("mrna", "transcript"):
        # Track transcript to gene mapping
        transcript_id = attrs.get("ID", "")
        parent_gene = attrs.get("Parent", "")
        if transcript_id and parent_gene:
            transcript_to_gene[transcript_id] = parent_gene
        return "", "", True  # Skip this record

    elif feature_type == "exon":
        # For exons, find the gene via parent transcript
        parent = attrs.get("Parent", "")
        gene_id = transcript_to_gene.get(parent, parent)
        gene_name = ""  # Will be filled from gene record
        return gene_id, gene_name, False

    else:
        gene_id = attrs.get("gene_id", attrs.get("Parent", ""))
        gene_name = attrs.get("gene_name", "")
        return gene_id, gene_name, False


def update_genes_data(
    genes_data: dict[str, dict],
    gene_id: str,
    gene_name: str,
    feature_type: str,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    attrs: dict[str, str] | None = None,
) -> None:
    """Update genes_data dictionary with gene or exon information."""
    if feature_type == "gene":
        if gene_id not in genes_data:
            genes_data[gene_id] = {
                "gene_name": gene_name,
                "chrom": chrom,
                "strand": strand,
                "exons": [],
                "start": start,
                "end": end,
            }
        else:
            if gene_name and not genes_data[gene_id]["gene_name"]:
                genes_data[gene_id]["gene_name"] = gene_name

    elif feature_type == "exon":
        if gene_id not in genes_data:
            genes_data[gene_id] = {
                "gene_name": gene_name or gene_id,
                "chrom": chrom,
                "strand": strand,
                "exons": [],
                "start": start,
                "end": end,
            }

        exon_number = None
        if attrs:
            raw = attrs.get("exon_number") or attrs.get("Exon_number") or attrs.get("rank") or ""
            exon_number_str = str(raw).strip() if raw is not None else ""
            if exon_number_str:
                try:
                    exon_number = int(exon_number_str)
                except (ValueError, TypeError):
                    exon_number = None

        genes_data[gene_id]["exons"].append((start, end, exon_number))
        genes_data[gene_id]["start"] = min(genes_data[gene_id]["start"], start)
        genes_data[gene_id]["end"] = max(genes_data[gene_id]["end"], end)


def convert_to_annotations(
    genes_data: dict[str, dict], region_start: int, region_end: int
) -> list[GeneAnnotation]:
    """Convert genes_data dictionary to list of GeneAnnotation objects."""
    annotations = []
    for gene_id, data in genes_data.items():
        exons = sorted(data["exons"], key=lambda exon: exon[0])
        if not exons:
            exons = [(data["start"], data["end"], None)]
        exons = [(s, e, n if n is not None else i) for i, (s, e, n) in enumerate(exons, start=1)]

        annotation = GeneAnnotation(
            gene_id=gene_id,
            gene_name=data["gene_name"] or gene_id,
            chrom=data["chrom"],
            start=data["start"],
            end=data["end"],
            strand=data["strand"],
            exons=exons,
        )

        if annotation.overlaps_region(region_start, region_end):
            annotations.append(annotation)

    annotations.sort(key=lambda annotation: annotation.start)
    return annotations


def parse_annotation_file(file_path: str, region: Region) -> list[GeneAnnotation]:
    """
    Parse indexed bgzip GTF/GFF3 and extract gene annotations overlapping the region.
    Uses tabix index for efficient region-based access.

    Args:
        file_path: Path to bgzip compressed GTF or GFF3 file with .tbi index
        region: Genomic region (chromosome, start, end, coordinate_str)

    Returns:
        List of GeneAnnotation objects for genes overlapping the region

    Raises:
        ValueError: If file is not bgzip or missing .tbi index
    """
    if not file_path:
        return []

    validate_bgzip_index(file_path)

    file_format = detect_annotation_format(file_path)
    parse_attrs = parse_gtf_attributes if file_format == "gtf" else parse_gff3_attributes

    genes_data: dict[str, dict] = {}
    transcript_to_gene: dict[str, str] = {}

    logger.debug(f"Parsing indexed {file_format.upper()} file: {file_path}")

    try:
        with pysam.TabixFile(file_path) as tabix_file:
            tabix_chrom = get_tabix_chromosome_name(tabix_file, region.chromosome)

            if tabix_chrom is None:
                logger.warning(
                    f"Chromosome {region.chromosome} not found in tabix index. "
                    f"Available contigs: {list(tabix_file.contigs)[:10]}"
                )
                return []

            # Tabix uses 0-based start, exclusive end (half-open interval)
            # GTF/GFF3 uses 1-based coordinates, inclusive on both ends
            tabix_start = region.start - 1
            tabix_end = region.end + 1

            try:
                records = tabix_file.fetch(tabix_chrom, tabix_start, tabix_end)
            except ValueError as err:
                logger.warning(
                    f"Could not fetch region {tabix_chrom}:{tabix_start}-{tabix_end} "
                    f"from tabix: {err}"
                )
                return []

            for line in records:
                if isinstance(line, bytes):
                    line = line.decode("utf-8")

                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                feature_type = parts[2].lower()
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attr_string = parts[8]

                if end < region.start or start > region.end:
                    continue

                attrs = parse_attrs(attr_string)

                if file_format == "gtf":
                    gene_id, gene_name, should_skip = extract_gtf_gene_info(attrs)
                else:
                    gene_id, gene_name, should_skip = extract_gff3_gene_info(
                        attrs, feature_type, transcript_to_gene
                    )

                if should_skip or not gene_id:
                    continue

                update_genes_data(
                    genes_data,
                    gene_id,
                    gene_name,
                    feature_type,
                    tabix_chrom,
                    start,
                    end,
                    strand,
                    attrs,
                )

        annotations = convert_to_annotations(genes_data, region.start, region.end)
        logger.debug(
            f"Found {len(annotations)} genes overlapping region "
            f"{region.chromosome}:{region.start}-{region.end}"
        )
        return annotations

    except FileNotFoundError:
        logger.error(f"Annotation file not found: {file_path}")
        return []
