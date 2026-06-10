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
    representative_transcript_id: str | None = None
    representative_transcript_name: str | None = None
    representative_selection_method: str = "gene annotation"

    def overlaps_region(self, region_start: int, region_end: int) -> bool:
        """Check if this gene overlaps with the given region."""
        return self.start <= region_end and self.end >= region_start


@dataclass
class TranscriptAnnotation:
    """Represents a transcript model with exon structure for IsoSeq visualization."""

    transcript_id: str
    transcript_name: str
    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: list[tuple[int, int, int]] = field(default_factory=list)

    def overlaps_region(self, region_start: int, region_end: int) -> bool:
        """Check if this transcript overlaps with the given region."""
        return self.start <= region_end and self.end >= region_start

    @property
    def intron_chain(self) -> tuple[tuple[int, int], ...]:
        """Return splice junction chain as (left_exon_end, right_exon_start)."""
        exons = sorted(self.exons, key=lambda exon: exon[0])
        return tuple((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1))


def parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GTF attribute column: ``key "value";`` pairs, then unquoted ``key N;`` numerics.
    """
    attrs = {}
    for match in re.finditer(r'(\S+)\s+"([^"]*)"', attr_string):
        key, value = match.groups()
        if key in attrs:
            attrs[key] = f"{attrs[key]},{value}"
            continue
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


REPRESENTATIVE_TRANSCRIPT_PRIORITY = (
    ("MANE Select", ("maneselect",)),
    ("RefSeq Select", ("refseqselect",)),
    ("Ensembl Canonical", ("ensemblcanonical",)),
    ("APPRIS principal", ("apprisprincipal",)),
    ("GENCODE Primary", ("gencodeprimary",)),
)


def _normalized_annotation_token(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def _representative_methods_from_attrs(attrs: dict[str, str]) -> set[str]:
    """Return representative-transcript methods advertised by annotation attributes."""
    tokens: set[str] = set()
    for key, value in attrs.items():
        tokens.add(_normalized_annotation_token(key))
        for part in str(value).split(","):
            tokens.add(_normalized_annotation_token(part))

    methods = set()
    for method, aliases in REPRESENTATIVE_TRANSCRIPT_PRIORITY:
        if any(alias in token for alias in aliases for token in tokens):
            methods.add(method)
    return methods


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
    transcripts_data: dict[str, dict] = {}
    transcript_to_gene: dict[str, str] = {}
    gene_names: dict[str, str] = {}

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
                    _update_gtf_transcripts_data(
                        transcripts_data, feature_type, attrs, tabix_chrom, start, end, strand
                    )
                else:
                    _update_gff3_transcripts_data(
                        transcripts_data,
                        transcript_to_gene,
                        gene_names,
                        feature_type,
                        attrs,
                        tabix_chrom,
                        start,
                        end,
                        strand,
                    )

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

        annotations = _convert_to_representative_gene_annotations(
            transcripts_data, region.start, region.end
        )
        if not annotations:
            annotations = convert_to_annotations(genes_data, region.start, region.end)
        logger.debug(
            f"Found {len(annotations)} genes overlapping region "
            f"{region.chromosome}:{region.start}-{region.end}"
        )
        return annotations

    except FileNotFoundError:
        logger.error(f"Annotation file not found: {file_path}")
        return []


def _update_gtf_transcripts_data(
    transcripts_data: dict[str, dict],
    feature_type: str,
    attrs: dict[str, str],
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> None:
    """Update transcript_data from a GTF feature."""
    transcript_id = attrs.get("transcript_id", "")
    if not transcript_id:
        return
    gene_id = attrs.get("gene_id", "")
    gene_name = attrs.get("gene_name", gene_id or "Unknown")
    transcript_name = attrs.get("transcript_name", transcript_id)
    selection_methods = _representative_methods_from_attrs(attrs)

    data = transcripts_data.setdefault(
        transcript_id,
        {
            "transcript_name": transcript_name,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "chrom": chrom,
            "strand": strand,
            "exons": [],
            "start": start,
            "end": end,
            "selection_methods": set(),
        },
    )
    data["selection_methods"].update(selection_methods)
    data["start"] = min(data["start"], start)
    data["end"] = max(data["end"], end)
    if gene_id and not data.get("gene_id"):
        data["gene_id"] = gene_id
    if gene_name and data.get("gene_name") in ("", "Unknown"):
        data["gene_name"] = gene_name
    if transcript_name and data.get("transcript_name") == transcript_id:
        data["transcript_name"] = transcript_name

    if feature_type == "exon":
        raw = attrs.get("exon_number") or attrs.get("Exon_number") or ""
        exon_number = None
        if raw:
            try:
                exon_number = int(str(raw).strip())
            except (ValueError, TypeError):
                exon_number = None
        data["exons"].append((start, end, exon_number))


def _update_gff3_transcripts_data(
    transcripts_data: dict[str, dict],
    transcript_to_gene: dict[str, str],
    gene_names: dict[str, str],
    feature_type: str,
    attrs: dict[str, str],
    chrom: str,
    start: int,
    end: int,
    strand: str,
) -> None:
    """Update transcript_data from a GFF3 feature."""
    if feature_type == "gene":
        gene_id = attrs.get("ID", attrs.get("gene_id", ""))
        if gene_id:
            gene_names[gene_id] = attrs.get("Name", attrs.get("gene_name", gene_id))
        return

    if feature_type in ("mrna", "transcript"):
        transcript_id = attrs.get("ID", "")
        gene_id = attrs.get("Parent", "")
        if not transcript_id:
            return
        if gene_id:
            transcript_to_gene[transcript_id] = gene_id
        transcripts_data.setdefault(
            transcript_id,
            {
                "transcript_name": attrs.get("Name", transcript_id),
                "gene_id": gene_id,
                "gene_name": gene_names.get(gene_id, gene_id or "Unknown"),
                "chrom": chrom,
                "strand": strand,
                "exons": [],
                "start": start,
                "end": end,
                "selection_methods": set(),
            },
        )
        transcripts_data[transcript_id]["selection_methods"].update(
            _representative_methods_from_attrs(attrs)
        )
        return

    if feature_type != "exon":
        return

    parent = attrs.get("Parent", "")
    if not parent:
        return
    transcript_id = parent.split(",")[0]
    gene_id = transcript_to_gene.get(transcript_id, "")
    data = transcripts_data.setdefault(
        transcript_id,
        {
            "transcript_name": transcript_id,
            "gene_id": gene_id,
            "gene_name": gene_names.get(gene_id, gene_id or "Unknown"),
            "chrom": chrom,
            "strand": strand,
            "exons": [],
            "start": start,
            "end": end,
            "selection_methods": set(),
        },
    )
    data["selection_methods"].update(_representative_methods_from_attrs(attrs))
    data["start"] = min(data["start"], start)
    data["end"] = max(data["end"], end)
    raw = attrs.get("exon_number") or attrs.get("rank") or ""
    exon_number = None
    if raw:
        try:
            exon_number = int(str(raw).strip())
        except (ValueError, TypeError):
            exon_number = None
    data["exons"].append((start, end, exon_number))


def _convert_to_transcript_annotations(
    transcripts_data: dict[str, dict], region_start: int, region_end: int
) -> list[TranscriptAnnotation]:
    """Convert transcript_data dictionary to TranscriptAnnotation objects."""
    annotations = []
    for transcript_id, data in transcripts_data.items():
        exons = sorted(data["exons"], key=lambda exon: exon[0])
        if not exons:
            continue
        exons = [(s, e, n if n is not None else i) for i, (s, e, n) in enumerate(exons, start=1)]
        annotation = TranscriptAnnotation(
            transcript_id=transcript_id,
            transcript_name=data["transcript_name"] or transcript_id,
            gene_id=data["gene_id"] or data["gene_name"] or "Unknown",
            gene_name=data["gene_name"] or data["gene_id"] or "Unknown",
            chrom=data["chrom"],
            start=min(s for s, _e, _n in exons),
            end=max(e for _s, e, _n in exons),
            strand=data["strand"],
            exons=exons,
        )
        if annotation.overlaps_region(region_start, region_end):
            annotations.append(annotation)
    annotations.sort(
        key=lambda transcript: (
            transcript.gene_name,
            transcript.start,
            transcript.transcript_id,
        )
    )
    return annotations


def _transcript_selection_priority(data: dict) -> tuple[int, str]:
    methods = data.get("selection_methods", set())
    for index, (method, _aliases) in enumerate(REPRESENTATIVE_TRANSCRIPT_PRIORITY):
        if method in methods:
            return index, method
    return len(REPRESENTATIVE_TRANSCRIPT_PRIORITY), "most exons"


def _number_transcript_exons(
    exons: list[tuple[int, int, int | None]],
    strand: str,
) -> list[tuple[int, int, int]]:
    sorted_exons = sorted(exons, key=lambda exon: exon[0])
    exon_count = len(sorted_exons)
    numbered_exons = []
    for index, (start, end, exon_number) in enumerate(sorted_exons):
        if exon_number is not None:
            numbered_exons.append((start, end, exon_number))
            continue
        display_number = exon_count - index if strand == "-" else index + 1
        numbered_exons.append((start, end, display_number))
    return numbered_exons


def _convert_to_representative_gene_annotations(
    transcripts_data: dict[str, dict],
    region_start: int,
    region_end: int,
) -> list[GeneAnnotation]:
    """Select one representative transcript per gene for the collapsed gene track."""
    transcripts_by_gene: dict[str, list[tuple[str, dict]]] = {}
    for transcript_id, data in transcripts_data.items():
        if not data["exons"]:
            continue
        gene_id = data["gene_id"] or data["gene_name"] or transcript_id
        transcripts_by_gene.setdefault(gene_id, []).append((transcript_id, data))

    annotations = []
    for gene_id, transcript_entries in transcripts_by_gene.items():
        selected_transcript_id, selected_data = min(
            transcript_entries,
            key=lambda item: (
                _transcript_selection_priority(item[1])[0],
                -len(item[1]["exons"]),
                item[1]["start"],
                item[0],
            ),
        )
        _priority, selection_method = _transcript_selection_priority(selected_data)
        exons = _number_transcript_exons(selected_data["exons"], selected_data["strand"])
        annotation = GeneAnnotation(
            gene_id=gene_id,
            gene_name=selected_data["gene_name"] or gene_id,
            chrom=selected_data["chrom"],
            start=selected_data["start"],
            end=selected_data["end"],
            strand=selected_data["strand"],
            exons=exons,
            representative_transcript_id=selected_transcript_id,
            representative_transcript_name=(
                selected_data["transcript_name"] or selected_transcript_id
            ),
            representative_selection_method=selection_method,
        )
        if annotation.overlaps_region(region_start, region_end):
            annotations.append(annotation)

    annotations.sort(key=lambda annotation: (annotation.start, annotation.gene_name))
    return annotations


def parse_transcript_annotation_file(file_path: str, region: Region) -> list[TranscriptAnnotation]:
    """Parse indexed bgzip GTF/GFF3 and extract transcript annotations for IsoSeq."""
    if not file_path:
        return []

    validate_bgzip_index(file_path)

    file_format = detect_annotation_format(file_path)
    parse_attrs = parse_gtf_attributes if file_format == "gtf" else parse_gff3_attributes
    transcripts_data: dict[str, dict] = {}
    transcript_to_gene: dict[str, str] = {}
    gene_names: dict[str, str] = {}

    logger.debug(f"Parsing indexed {file_format.upper()} transcripts: {file_path}")

    try:
        with pysam.TabixFile(file_path) as tabix_file:
            tabix_chrom = get_tabix_chromosome_name(tabix_file, region.chromosome)
            if tabix_chrom is None:
                logger.warning(
                    f"Chromosome {region.chromosome} not found in tabix index. "
                    f"Available contigs: {list(tabix_file.contigs)[:10]}"
                )
                return []
            try:
                records = tabix_file.fetch(tabix_chrom, region.start - 1, region.end + 1)
            except ValueError as err:
                logger.warning(
                    f"Could not fetch region {tabix_chrom}:{region.start - 1}-{region.end + 1} "
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
                if feature_type not in ("gene", "transcript", "mrna", "exon"):
                    continue
                start = int(parts[3])
                end = int(parts[4])
                if end < region.start or start > region.end:
                    continue
                attrs = parse_attrs(parts[8])
                if file_format == "gtf":
                    _update_gtf_transcripts_data(
                        transcripts_data, feature_type, attrs, tabix_chrom, start, end, parts[6]
                    )
                else:
                    _update_gff3_transcripts_data(
                        transcripts_data,
                        transcript_to_gene,
                        gene_names,
                        feature_type,
                        attrs,
                        tabix_chrom,
                        start,
                        end,
                        parts[6],
                    )
        annotations = _convert_to_transcript_annotations(transcripts_data, region.start, region.end)
        logger.debug(
            f"Found {len(annotations)} transcripts overlapping region "
            f"{region.chromosome}:{region.start}-{region.end}"
        )
        return annotations
    except FileNotFoundError:
        logger.error(f"Annotation file not found: {file_path}")
        return []
