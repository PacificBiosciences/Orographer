"""
Python implementation of Rust BAM SA parser functions.
Implements data structures and parsing functions for split alignments.
"""

import re
from collections import defaultdict
from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path

import pysam
from intervaltree import Interval, IntervalTree

from .utils import Region

CoverageHaplotype = int | str


@dataclass
class SplitAlignmentSegment:
    """Represents a split alignment segment from SA tag."""

    rname: str
    pos: int  # 0-based position
    is_fwd_strand: bool
    cigar: str
    mapq: int


@dataclass
class ReadMetadata:
    """Metadata for a read including phasing information."""

    readname: str
    phaseset_tag: int | None
    haplotype_tag: int | None
    color_tag: str | None  # YC tag - RGB color string like "255,0,0"


@dataclass
class FwdStrandReadSegment:
    """Forward strand read segment with complete alignment information."""

    fwd_read_start: int
    fwd_read_end: int
    chrom: str
    second_chrom: str
    pos: int
    end: int
    is_fwd_strand: bool
    is_start_softclipped: bool
    is_end_softclipped: bool
    phaseset_tag: int | None
    haplotype_tag: int | None
    color_tag: str | None  # YC tag - RGB color string like "255,0,0"
    spans: bool
    from_primary_bam_record: bool
    readname: str
    mismatches: list[tuple[int, str]] = None  # List of (ref_position, alt_base) tuples
    insertions: list[tuple[int, str]] = None  # List of (ref_position, inserted_bases) tuples
    deletions: list[tuple[int, int]] = None  # List of (ref_start, ref_end) tuples
    aligned_blocks: list[tuple[int, int]] = None  # Reference blocks, 0-based half-open
    alignment_order: int = 0  # 1-based order within the read (set in order_alignments)
    inclusion_reason: str | None = None

    def __post_init__(self):
        if self.mismatches is None:
            self.mismatches = []
        if self.insertions is None:
            self.insertions = []
        if self.deletions is None:
            self.deletions = []
        if self.aligned_blocks is None:
            self.aligned_blocks = [(self.pos, self.end)]


@dataclass(slots=True)
class AlignmentSummary:
    """Cheap discovery-time record for one BAM alignment. No FASTA, no variant extraction."""

    read_name: str
    chromosome: str
    start: int  # reference_start (0-based)
    end: int  # reference_end
    is_supplementary: bool
    has_sa: bool
    sa_tag: str | None  # raw SA string; doubles as SA parse-cache key
    max_indel: int = 0  # longest I or D op in cigartuples; 0 if unset or unmapped
    max_softclip: int = 0  # longest S op in cigartuples; 0 if unset or unmapped


@dataclass
class ExcludeRegions:
    """Container for regions to exclude from analysis."""

    regions: dict[str, set]  # Simplified from Rust BTreeSet


@dataclass
class SegmentProcessingContext:
    """Context for processing alignment segments."""

    fwd_read_split_segments: list[FwdStrandReadSegment]
    exclude_regions: ExcludeRegions
    excluded: bool


# Constants
MIN_MAPQ = 5
MIN_SOFTCLIPS = 10


def parse_sa_segment(seg: str) -> SplitAlignmentSegment:
    """Parse a single SA segment string into a SplitAlignmentSegment."""
    sa_fields = seg.split(",")
    if len(sa_fields) != 6:
        raise ValueError(f"Unexpected segment in bam SA tag: {seg}")

    rname = sa_fields[0]
    pos = int(sa_fields[1]) - 1  # Convert to 0-based
    is_fwd_strand = sa_fields[2] == "+"
    cigar = sa_fields[3]
    mapq = int(sa_fields[4])

    return SplitAlignmentSegment(
        rname=rname, pos=pos, is_fwd_strand=is_fwd_strand, cigar=cigar, mapq=mapq
    )


def parse_sa_aux_val(sa_aux_val: str) -> list[SplitAlignmentSegment]:
    """Parse SA auxiliary value into list of SplitAlignmentSegment objects."""
    if not sa_aux_val:
        return []

    segments = sa_aux_val.split(";")
    return [parse_sa_segment(seg) for seg in segments if seg]


def get_optional_string_aux_tag(record: pysam.AlignedSegment, aux_tag: str) -> str | None:
    """Retrieve a string aux tag from BAM record."""
    try:
        return record.get_tag(aux_tag)
    except KeyError:
        return None


def get_optional_int_aux_tag(record: pysam.AlignedSegment, aux_tag: str) -> int | None:
    """Retrieve an integer aux tag from BAM record."""
    try:
        return record.get_tag(aux_tag)
    except KeyError:
        return None


def extract_variants(
    record: pysam.AlignedSegment, ref_fasta: pysam.FastaFile = None
) -> tuple[list[tuple[int, str]], list[tuple[int, str]], list[tuple[int, int]]]:
    """
    Extract mismatches, insertions, and deletions from an alignment.

    Args:
        record: pysam.AlignedSegment object
        ref_fasta: pysam.FastaFile object for fetching reference bases

    Returns:
        Tuple of (mismatches, insertions, deletions) where:
        - mismatches: List of (ref_position, alt_base) tuples
        - insertions: List of (ref_position, inserted_bases)
        - deletions: List of (ref_start, ref_end) tuples (half-open interval)
    """
    mismatches = []
    insertions = []
    deletions = []

    if ref_fasta is None:
        return mismatches, insertions, deletions

    query_seq = record.query_sequence
    if query_seq is None:
        return mismatches, insertions, deletions

    chrom = record.reference_name
    if chrom is None:
        return mismatches, insertions, deletions

    ref_start = record.reference_start
    ref_end = record.reference_end
    if ref_start is None or ref_end is None:
        return mismatches, insertions, deletions

    # Fetch the entire reference region at once (much faster than base-by-base)
    ref_seq = ref_fasta.fetch(chrom, ref_start, ref_end).upper()

    # Get aligned pairs (query_pos, ref_pos); no with_seq to avoid MD tag
    aligned_pairs = record.get_aligned_pairs(matches_only=False)
    insertion_query_positions = insertion_query_positions_from_cigar(record.cigartuples)

    # Track current insertion/deletion runs
    current_insertion_start_ref = None
    current_insertion_bases = []
    current_deletion_start = None
    last_ref_pos = None

    for query_pos, ref_pos in aligned_pairs:
        # Handle insertions (query has base, reference doesn't)
        if ref_pos is None and query_pos is not None and query_pos in insertion_query_positions:
            if current_insertion_start_ref is None:
                current_insertion_start_ref = (
                    last_ref_pos if last_ref_pos is not None else ref_start
                )
            current_insertion_bases.append(query_seq[query_pos].upper())
            continue
        else:
            # End of insertion run - save it
            if current_insertion_bases:
                insertions.append((current_insertion_start_ref, "".join(current_insertion_bases)))
                current_insertion_bases = []
                current_insertion_start_ref = None

        # Handle deletions (reference has base, query doesn't)
        if query_pos is None and ref_pos is not None:
            if current_deletion_start is None:
                current_deletion_start = ref_pos
            last_ref_pos = ref_pos
            continue
        else:
            # End of deletion run - save it
            if current_deletion_start is not None:
                deletions.append((current_deletion_start, last_ref_pos + 1))
                current_deletion_start = None

        # Both have values - check for mismatch
        if query_pos is not None and ref_pos is not None:
            last_ref_pos = ref_pos

            read_base = query_seq[query_pos].upper()

            # Get reference base by indexing into the fetched sequence
            ref_offset = ref_pos - ref_start
            if ref_offset < 0 or ref_offset >= len(ref_seq):
                continue
            ref_base = ref_seq[ref_offset]

            # Check for mismatch (skip if bases match or if ref_base is N)
            if read_base != ref_base and ref_base != "N":
                mismatches.append((ref_pos, read_base))

    # Handle any trailing insertion
    if current_insertion_bases:
        insertions.append((current_insertion_start_ref, "".join(current_insertion_bases)))

    # Handle any trailing deletion
    if current_deletion_start is not None and last_ref_pos is not None:
        deletions.append((current_deletion_start, last_ref_pos + 1))

    return mismatches, insertions, deletions


def insertion_query_positions_from_cigar(cigartuples: list[tuple[int, int]] | None) -> set[int]:
    """Return query positions that belong to true CIGAR insertion operations."""
    insertion_positions: set[int] = set()
    query_pos = 0
    for op, length in cigartuples or []:
        if op == 1:
            insertion_positions.update(range(query_pos, query_pos + length))
        if op in (0, 1, 4, 7, 8):
            query_pos += length
    return insertion_positions


def get_cigarseg_ref_offset(cigar_op: tuple) -> int:
    """Get reference offset for a CIGAR operation."""
    op, length = cigar_op
    if isinstance(op, int):
        # Convert numeric code to string operation
        op_map = {
            0: "M",
            1: "I",
            2: "D",
            3: "N",
            4: "S",
            5: "H",
            6: "P",
            7: "=",
            8: "X",
        }
        op = op_map.get(op, "M")

    if op in ["D", "N", "X", "=", "M"]:
        return length
    return 0


def get_cigarseg_complete_read_offset(cigar_op: tuple) -> int:
    """Get read offset for a CIGAR operation."""
    op, length = cigar_op
    if isinstance(op, int):
        # Convert numeric code to string operation
        op_map = {
            0: "M",
            1: "I",
            2: "D",
            3: "N",
            4: "S",
            5: "H",
            6: "P",
            7: "=",
            8: "X",
        }
        op = op_map.get(op, "M")

    if op in ["H", "I", "S", "X", "=", "M"]:
        return length
    return 0


def update_ref_and_complete_read_pos(
    cigar_op: tuple[str, int], ref_pos: int, read_pos: int
) -> tuple[int, int]:
    """Update reference and read positions based on CIGAR operation."""
    ref_offset = get_cigarseg_ref_offset(cigar_op)
    read_offset = get_cigarseg_complete_read_offset(cigar_op)
    return ref_pos + ref_offset, read_pos + read_offset


def get_complete_read_clip_positions(cigar: list[tuple]) -> tuple[int, int, int]:
    """Get complete read clip positions from CIGAR string."""
    ref_pos = 0
    read_pos = 0
    left_clip_size = 0
    right_clip_size = 0
    left_clip = True

    # Convert numeric CIGAR codes to string operations if needed
    cigar_ops = []
    for op, length in cigar:
        if isinstance(op, int):
            # Convert numeric code to string operation
            op_map = {
                0: "M",
                1: "I",
                2: "D",
                3: "N",
                4: "S",
                5: "H",
                6: "P",
                7: "=",
                8: "X",
            }
            op = op_map.get(op, "M")
        cigar_ops.append((op, length))

    for op, length in cigar_ops:
        if op in ["H", "S"]:
            if left_clip:
                left_clip_size += length
            else:
                right_clip_size += length
        else:
            left_clip = False

        ref_pos, read_pos = update_ref_and_complete_read_pos((op, length), ref_pos, read_pos)

    return left_clip_size, read_pos - right_clip_size, read_pos


def get_sa_softclip_status(cigar: list[tuple], min_softclips: int) -> tuple[bool, bool]:
    """Check if start and end of CIGAR string are softclipped."""
    if not cigar:
        return False, False

    # Convert numeric codes to string operations if needed
    first_op, first_length = cigar[0]
    last_op, last_length = cigar[-1]

    if isinstance(first_op, int):
        op_map = {
            0: "M",
            1: "I",
            2: "D",
            3: "N",
            4: "S",
            5: "H",
            6: "P",
            7: "=",
            8: "X",
        }
        first_op = op_map.get(first_op, "M")
    if isinstance(last_op, int):
        op_map = {
            0: "M",
            1: "I",
            2: "D",
            3: "N",
            4: "S",
            5: "H",
            6: "P",
            7: "=",
            8: "X",
        }
        last_op = op_map.get(last_op, "M")

    is_start_softclipped = first_op == "S" and first_length > min_softclips
    is_end_softclipped = last_op == "S" and last_length > min_softclips

    return is_start_softclipped, is_end_softclipped


def get_fwd_strand_read_pos(
    read_start: int, read_end: int, read_size: int, is_fwd_strand: bool
) -> tuple[int, int]:
    """Get forward strand read positions."""
    if is_fwd_strand:
        return read_start, read_end
    else:
        return read_size - read_end, read_size - read_start


def get_seq_len_from_cigar(cigar: list[tuple]) -> int:
    """Get sequence length from CIGAR string."""
    aligned_len = 0
    for op, length in cigar:
        if isinstance(op, int):
            # Convert numeric code to string operation
            op_map = {
                0: "M",
                1: "I",
                2: "D",
                3: "N",
                4: "S",
                5: "H",
                6: "P",
                7: "=",
                8: "X",
            }
            op = op_map.get(op, "M")

        if op in ["M", "=", "D", "N", "X"]:
            aligned_len += length
    return aligned_len


def process_primary_alignment(
    record: pysam.AlignedSegment,
    chrom: str,
    read_metadata: ReadMetadata,
    context: SegmentProcessingContext,
    ref_fasta: pysam.FastaFile = None,
) -> int:
    """Process the primary alignment record and create a FwdStrandReadSegment."""
    cigar = list(record.cigartuples)
    read_start, read_end, read_size = get_complete_read_clip_positions(cigar)
    fwd_read_start, fwd_read_end = get_fwd_strand_read_pos(
        read_start, read_end, read_size, not record.is_reverse
    )

    aligned_len = get_seq_len_from_cigar(cigar)
    is_start_softclipped, is_end_softclipped = get_sa_softclip_status(cigar, MIN_SOFTCLIPS)

    # Extract variants (mismatches, insertions, deletions) from the alignment
    mismatches, insertions, deletions = extract_variants(record, ref_fasta)

    alignment = FwdStrandReadSegment(
        fwd_read_start=fwd_read_start,
        fwd_read_end=fwd_read_end,
        chrom=chrom,
        second_chrom=chrom,
        pos=record.reference_start,
        end=record.reference_start + aligned_len,
        is_fwd_strand=not record.is_reverse,
        is_start_softclipped=is_start_softclipped,
        is_end_softclipped=is_end_softclipped,
        phaseset_tag=read_metadata.phaseset_tag,
        haplotype_tag=read_metadata.haplotype_tag,
        color_tag=read_metadata.color_tag,
        from_primary_bam_record=True,
        readname=read_metadata.readname,
        spans=True,
        mismatches=mismatches,
        insertions=insertions,
        deletions=deletions,
        aligned_blocks=record.get_blocks(),
    )

    context.fwd_read_split_segments.append(alignment)
    return read_size


def process_sa_segments(
    sa_segments: list[SplitAlignmentSegment],
    primary_read_size: int,
    read_metadata: ReadMetadata,
    context: SegmentProcessingContext,
) -> None:
    """Process supplementary alignment segments from the SA tag."""
    for sa_segment in sa_segments:
        if sa_segment.mapq < MIN_MAPQ:
            context.excluded = True
            continue

        # Parse CIGAR string using regex for better reliability
        cigar = []
        cigar_str = sa_segment.cigar

        # Use regex to find all number+letter pairs
        cigar_pairs = re.findall(r"(\d+)([MIDNSHPX=])", cigar_str)
        for length_str, op in cigar_pairs:
            cigar.append((op, int(length_str)))

        # If CIGAR parsing failed, skip this segment
        if not cigar:
            continue

        read_start, read_end, read_size = get_complete_read_clip_positions(cigar)
        if read_size != primary_read_size:
            continue  # Skip if read sizes don't match

        fwd_read_start, fwd_read_end = get_fwd_strand_read_pos(
            read_start, read_end, read_size, sa_segment.is_fwd_strand
        )

        aligned_len = get_seq_len_from_cigar(cigar)
        if not cigar:
            continue

        is_start_softclipped, is_end_softclipped = get_sa_softclip_status(cigar, MIN_SOFTCLIPS)

        # SA tag does not carry HP/PS; segment gets None so we set it when we see the
        # supplementary BAM record (which can have its own haplotype assignment).
        alignment = FwdStrandReadSegment(
            fwd_read_start=fwd_read_start,
            fwd_read_end=fwd_read_end,
            chrom=sa_segment.rname,
            second_chrom=sa_segment.rname,
            pos=sa_segment.pos,
            end=sa_segment.pos + aligned_len,
            is_fwd_strand=sa_segment.is_fwd_strand,
            is_start_softclipped=is_start_softclipped,
            is_end_softclipped=is_end_softclipped,
            phaseset_tag=None,
            haplotype_tag=None,
            color_tag=read_metadata.color_tag,
            from_primary_bam_record=False,
            readname=read_metadata.readname,
            spans=True,
        )

        context.fwd_read_split_segments.append(alignment)


def get_fwd_read_split_segments(
    record: pysam.AlignedSegment,
    chrom: str,
    exclude_regions: ExcludeRegions,
    ref_fasta: pysam.FastaFile = None,
) -> list[FwdStrandReadSegment]:
    """Main function to get forward read split segments from a BAM record."""
    fwd_read_split_segments = []
    excluded = False

    readname = record.query_name

    # Check minimum mapping quality
    primary_mapq = record.mapping_quality
    if primary_mapq < MIN_MAPQ:
        return fwd_read_split_segments

    # Check for SA tag
    sa_aux_val = get_optional_string_aux_tag(record, "SA")
    if sa_aux_val:
        sa_segments = parse_sa_aux_val(sa_aux_val)
        phaseset_tag = get_optional_int_aux_tag(record, "PS")
        haplotype_tag = get_optional_int_aux_tag(record, "HP")
        color_tag = get_optional_string_aux_tag(record, "YC")

        # Create read metadata
        read_metadata = ReadMetadata(
            readname=readname,
            phaseset_tag=phaseset_tag,
            haplotype_tag=haplotype_tag,
            color_tag=color_tag,
        )

        # Create processing context
        context = SegmentProcessingContext(
            fwd_read_split_segments=fwd_read_split_segments,
            exclude_regions=exclude_regions,
            excluded=excluded,
        )

        # Process primary alignment (with reference for mismatch extraction)
        primary_read_size = process_primary_alignment(
            record, chrom, read_metadata, context, ref_fasta
        )

        # Process SA segments (no mismatch extraction for supplementary alignments)
        process_sa_segments(sa_segments, primary_read_size, read_metadata, context)

    return fwd_read_split_segments


def convert_alignment_to_segment(
    record: pysam.AlignedSegment, chrom: str, ref_fasta: pysam.FastaFile = None
) -> FwdStrandReadSegment | None:
    """
    Convert a single pysam alignment to a FwdStrandReadSegment.

    Args:
        record: pysam.AlignedSegment object
        chrom: Chromosome name
        ref_fasta: pysam.FastaFile for mismatch extraction

    Returns:
        FwdStrandReadSegment or None if alignment doesn't meet quality criteria
    """
    # Check minimum mapping quality
    if record.mapping_quality < MIN_MAPQ:
        return None

    readname = record.query_name
    phaseset_tag = get_optional_int_aux_tag(record, "PS")
    haplotype_tag = get_optional_int_aux_tag(record, "HP")
    color_tag = get_optional_string_aux_tag(record, "YC")

    cigar = list(record.cigartuples)
    read_start, read_end, read_size = get_complete_read_clip_positions(cigar)
    fwd_read_start, fwd_read_end = get_fwd_strand_read_pos(
        read_start, read_end, read_size, not record.is_reverse
    )

    aligned_len = get_seq_len_from_cigar(cigar)
    is_start_softclipped, is_end_softclipped = get_sa_softclip_status(cigar, MIN_SOFTCLIPS)

    # Extract variants (mismatches, insertions, deletions) from the alignment
    mismatches, insertions, deletions = extract_variants(record, ref_fasta)

    return FwdStrandReadSegment(
        fwd_read_start=fwd_read_start,
        fwd_read_end=fwd_read_end,
        chrom=chrom,
        second_chrom=chrom,
        pos=record.reference_start,
        end=record.reference_start + aligned_len,
        is_fwd_strand=not record.is_reverse,
        is_start_softclipped=is_start_softclipped,
        is_end_softclipped=is_end_softclipped,
        phaseset_tag=phaseset_tag,
        haplotype_tag=haplotype_tag,
        color_tag=color_tag,
        from_primary_bam_record=True,
        readname=readname,
        spans=True,
        mismatches=mismatches,
        insertions=insertions,
        deletions=deletions,
        aligned_blocks=record.get_blocks(),
    )


def fetch_all_alignments(
    bam_path,
    region: Region,
    only_split=False,
    reference_path=None,
    read_name_filter=None,
    record_callback=None,
):
    """
    Fetch alignments from BAM file for the specified region.

    Args:
        bam_path (str): Path to the BAM file
        region (Region): Genomic region (chromosome, start, end, coordinate_str)
        only_split (bool): If True, only fetch reads with SA tags (split alignments).
                          If False, fetch all reads (split and non-split).
        reference_path (str): Path to reference FASTA file for mismatch extraction.
        read_name_filter (set[str] | None): Optional raw read names to include.
        record_callback: Optional callable invoked once per fetched, non-secondary,
                         de-duplicated record before segment conversion.

    Returns:
        dict: Read name -> list of FwdStrandReadSegment
    """
    try:
        bam_file = pysam.AlignmentFile(bam_path, "rb")

        # Open reference FASTA for mismatch extraction
        ref_fasta = None
        if reference_path:
            ref_fasta = pysam.FastaFile(reference_path)

        start_0based = region.start - 1
        alignments = bam_file.fetch(region.chromosome, start_0based, region.end)

        # Initialize exclude regions (empty for now)
        exclude_regions = ExcludeRegions(regions={})

        # Segments by read name (one list per read; haplotype split at display)
        segments_by_read = defaultdict(list)

        # Track (read_name, ref_start) so each alignment added once per haplotype
        processed_alignments = set()

        for alignment in alignments:
            # Skip secondary alignments early (fast check)
            if alignment.is_secondary:
                continue

            # Early exit for only_split mode (fast check)
            if only_split and not alignment.has_tag("SA"):
                continue

            read_name = alignment.query_name
            if read_name_filter is not None and read_name not in read_name_filter:
                continue
            ref_start = alignment.reference_start

            alignment_key = (read_name, ref_start)
            if alignment_key in processed_alignments:
                continue
            processed_alignments.add(alignment_key)
            if record_callback is not None:
                record_callback(alignment)

            # Supplementary alignments can have their own HP/PS; update existing segment
            # (from primary's SA parse) or add one if we saw this supp first.
            if alignment.is_supplementary:
                supp_hp = get_optional_int_aux_tag(alignment, "HP")
                supp_ps = get_optional_int_aux_tag(alignment, "PS")
                existing = [
                    s for s in segments_by_read[read_name] if s.pos == alignment.reference_start
                ]
                if existing:
                    for s in existing:
                        s.haplotype_tag = supp_hp
                        s.phaseset_tag = supp_ps
                        s.aligned_blocks = alignment.get_blocks() or [(s.pos, s.end)]
                    if ref_fasta:
                        mismatches, insertions, deletions = extract_variants(alignment, ref_fasta)
                        for s in existing:
                            s.mismatches = list(mismatches)
                            s.insertions = list(insertions)
                            s.deletions = list(deletions)
                else:
                    # No existing segment (e.g. supplementary seen before primary). Add segment
                    # with variant annotations so non-primary alignments get SNV/INDEL drawn.
                    segment = convert_alignment_to_segment(
                        alignment, alignment.reference_name, ref_fasta
                    )
                    if segment:
                        segment.from_primary_bam_record = False
                        segments_by_read[read_name].append(segment)
                continue

            # Process alignment (store under read_name, one list per read)
            if alignment.has_tag("SA"):
                split_segments = get_fwd_read_split_segments(
                    alignment, region.chromosome, exclude_regions, ref_fasta
                )
                if split_segments:
                    existing_positions = {s.pos for s in segments_by_read[read_name]}
                    for seg in split_segments:
                        if seg.from_primary_bam_record:
                            segments_by_read[read_name].append(seg)
                        elif seg.pos not in existing_positions:
                            segments_by_read[read_name].append(seg)
                            existing_positions.add(seg.pos)
            elif not only_split:
                segment = convert_alignment_to_segment(alignment, region.chromosome, ref_fasta)
                if segment:
                    segments_by_read[read_name].append(segment)

        bam_file.close()
        if ref_fasta:
            ref_fasta.close()

        # Order per read, then group by haplotype; each segment under its haplotype_tag.
        order_alignments(segments_by_read)
        segments_by_haplotype = {}
        for rname, segs in segments_by_read.items():
            for seg in segs:
                hp = seg.haplotype_tag if seg.haplotype_tag is not None else 0
                key = rname if hp == 0 else f"{rname}_HP{hp}"
                segments_by_haplotype.setdefault(key, []).append(seg)
        for key in segments_by_haplotype:
            segments_by_haplotype[key].sort(key=lambda segment: segment.alignment_order)
        return dict(segments_by_haplotype)
    except FileNotFoundError as err:
        raise FileNotFoundError(f"BAM file not found or not readable: {bam_path}") from err
    except ValueError as err:
        raise ValueError(f"Invalid coordinates: {err}") from err
    except OSError as err:
        raise RuntimeError(f"Error reading BAM file: {err}") from err


def build_alignments_from_records(
    records,
    region: Region,
    reference_path=None,
    read_name_filter=None,
):
    """Build display segments from already-fetched BAM records.

    This mirrors fetch_all_alignments() without performing a second indexed BAM
    fetch, which is useful when a caller needs full variant-bearing segments for
    a sampled subset after a lightweight first pass.
    """
    ref_fasta = pysam.FastaFile(reference_path) if reference_path else None
    try:
        exclude_regions = ExcludeRegions(regions={})
        segments_by_read = defaultdict(list)
        processed_alignments = set()
        for alignment in records:
            if alignment.is_secondary:
                continue
            read_name = alignment.query_name
            if read_name_filter is not None and read_name not in read_name_filter:
                continue
            ref_start = alignment.reference_start
            alignment_key = (read_name, ref_start)
            if alignment_key in processed_alignments:
                continue
            processed_alignments.add(alignment_key)

            if alignment.is_supplementary:
                supp_hp = get_optional_int_aux_tag(alignment, "HP")
                supp_ps = get_optional_int_aux_tag(alignment, "PS")
                existing = [
                    s for s in segments_by_read[read_name] if s.pos == alignment.reference_start
                ]
                if existing:
                    for s in existing:
                        s.haplotype_tag = supp_hp
                        s.phaseset_tag = supp_ps
                        s.aligned_blocks = alignment.get_blocks() or [(s.pos, s.end)]
                    if ref_fasta:
                        mismatches, insertions, deletions = extract_variants(alignment, ref_fasta)
                        for s in existing:
                            s.mismatches = list(mismatches)
                            s.insertions = list(insertions)
                            s.deletions = list(deletions)
                else:
                    segment = convert_alignment_to_segment(
                        alignment, alignment.reference_name, ref_fasta
                    )
                    if segment:
                        segment.from_primary_bam_record = False
                        segments_by_read[read_name].append(segment)
                continue

            if alignment.has_tag("SA"):
                split_segments = get_fwd_read_split_segments(
                    alignment, region.chromosome, exclude_regions, ref_fasta
                )
                if split_segments:
                    existing_positions = {s.pos for s in segments_by_read[read_name]}
                    for seg in split_segments:
                        if seg.from_primary_bam_record:
                            segments_by_read[read_name].append(seg)
                        elif seg.pos not in existing_positions:
                            segments_by_read[read_name].append(seg)
                            existing_positions.add(seg.pos)
            else:
                segment = convert_alignment_to_segment(alignment, region.chromosome, ref_fasta)
                if segment:
                    segments_by_read[read_name].append(segment)

        order_alignments(segments_by_read)
        segments_by_haplotype = {}
        for rname, segs in segments_by_read.items():
            for seg in segs:
                hp = seg.haplotype_tag if seg.haplotype_tag is not None else 0
                key = rname if hp == 0 else f"{rname}_HP{hp}"
                segments_by_haplotype.setdefault(key, []).append(seg)
        for key in segments_by_haplotype:
            segments_by_haplotype[key].sort(key=lambda segment: segment.alignment_order)
        return dict(segments_by_haplotype)
    finally:
        if ref_fasta:
            ref_fasta.close()


def validate_bam_file(bam_path):
    """
    Validate that the BAM file exists and is readable.

    Args:
        bam_path (str): Path to the BAM file

    Returns:
        bool: True if file exists and is readable

    Raises:
        FileNotFoundError: If BAM file doesn't exist
    """
    bam_file = Path(bam_path)
    if not bam_file.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    if not bam_file.is_file():
        raise FileNotFoundError(f"Path is not a file: {bam_path}")

    return True


def order_split_alignments(alignments):
    """
    Robustly order split alignments (primary + supplementary) of a read in query space.

    Parameters:
        alignments (list): pysam.AlignedSegment objects with same query_name.

    Returns:
        List of alignments ordered by their position on the original read.
    """
    if not alignments:
        return []

    # Prefer primary alignment as anchor
    primary = next(
        (aln for aln in alignments if not aln.is_secondary and not aln.is_supplementary),
        None,
    )

    # Try SA tag-based ordering if available
    if primary and primary.has_tag("SA"):
        sa_tag = primary.get_tag("SA").strip(";")
        sa_entries = sa_tag.split(";")

        # Build ordered keys from SA tag (skip malformed / incomplete entries)
        sa_keys = []
        for entry in sa_entries:
            entry = entry.strip()
            if not entry:
                continue
            parts = entry.split(",")
            if len(parts) < 6:
                continue
            try:
                rname, pos_s, strand = parts[0], parts[1], parts[2]
                pos = int(pos_s)
            except (ValueError, IndexError):
                continue
            sa_keys.append((rname, pos, strand))

        if sa_keys:

            def sa_sort_key(aln):
                key = (
                    aln.reference_name,
                    aln.reference_start,
                    "-" if aln.is_reverse else "+",
                )
                try:
                    return sa_keys.index(key)
                except ValueError:
                    return len(sa_keys)  # Not found → put last

            return sorted(alignments, key=sa_sort_key)

    # Fallback: Use query coordinate intervals (with overlap detection)
    def read_span(aln):
        return min(aln.query_alignment_start, aln.query_alignment_end), max(
            aln.query_alignment_start, aln.query_alignment_end
        )

    spans = [(read_span(aln), aln) for aln in alignments]

    # Build interval tree for overlap detection
    tree = IntervalTree()
    for (start, end), aln in spans:
        tree.add(Interval(start, end, aln))

    # Sort alignments in left-to-right read order, breaking ties by:
    # - lower start
    # - longer alignment
    # - primary > supplementary > secondary
    def sort_key(aln):
        start, end = read_span(aln)
        return (
            start,
            -(end - start),
            (
                0
                if not aln.is_secondary and not aln.is_supplementary
                else (1 if aln.is_supplementary else 2)
            ),
        )

    # Sort alignments by that key
    return sorted(alignments, key=sort_key)


def fetch_alignment_summaries(
    bam_file: pysam.AlignmentFile,
    region: Region,
    mapq_threshold: int = MIN_MAPQ,
) -> list[AlignmentSummary]:
    """
    Cheap discovery-time fetch: returns AlignmentSummary objects without opening
    the reference FASTA, extracting variants, or building FwdStrandReadSegment objects.

    Args:
        bam_file: Open pysam.AlignmentFile (caller owns the handle).
        region: Genomic region to fetch (start is 1-based inclusive, converted internally).
        mapq_threshold: Minimum mapping quality; alignments below this are skipped.

    Returns:
        List of AlignmentSummary, one per qualifying alignment (secondaries excluded).
    """
    summaries: list[AlignmentSummary] = []
    start_0based = region.start - 1
    for record in bam_file.fetch(region.chromosome, start_0based, region.end):
        if record.is_secondary:
            continue
        if record.mapping_quality < mapq_threshold:
            continue
        if record.reference_name is None or record.reference_start is None:
            continue
        has_sa = record.has_tag("SA")
        sa_tag: str | None = record.get_tag("SA") if has_sa else None
        cigar = record.cigartuples or []
        max_indel = max((length for op, length in cigar if op in (1, 2)), default=0)
        max_softclip = max((length for op, length in cigar if op == 4), default=0)
        summaries.append(
            AlignmentSummary(
                read_name=record.query_name,
                chromosome=record.reference_name,
                start=record.reference_start,
                end=record.reference_end or (record.reference_start + 1),
                is_supplementary=record.is_supplementary,
                has_sa=has_sa,
                sa_tag=sa_tag,
                max_indel=max_indel,
                max_softclip=max_softclip,
            )
        )
    return summaries


def _merged_coverage_blocks(record, min_gap: int):
    """Yield (start, end) coverage blocks from record.get_blocks(), merging
    adjacent blocks whose intervening gap is smaller than *min_gap* bp.

    Small indels (< min_gap) are treated as covered; large gaps are preserved
    as dips in the depth track.
    """
    blocks = record.get_blocks()
    if not blocks:
        return
    cur_start, cur_end = blocks[0]
    for block_start, block_end in blocks[1:]:
        if block_start - cur_end < min_gap:
            cur_end = block_end
        else:
            yield cur_start, cur_end
            cur_start, cur_end = block_start, block_end
    yield cur_start, cur_end


def compute_coverage(
    bam_path: str,
    region: Region,
    haplotypes_to_track: set[CoverageHaplotype] | None = None,
    include_observed_haplotypes: bool = False,
    bin_size: int = 10,
    mapq_threshold: int = MIN_MAPQ,
    min_cigar_gap: int = 20,
) -> dict[CoverageHaplotype, tuple[list[int], list[int]]]:
    """
    Compute binned coverage for a region using a difference-array approach.

    Adjacent CIGAR-aligned blocks are merged when the gap between them is smaller
    than *min_cigar_gap* bp, smoothing over small indels while preserving large
    structural gaps as dips in the depth track.
    Key -1 is the total-all-reads series. Other keys are HP tag values.
    Records without HP are included in total but not emitted as a separate series
    unless 0 is explicitly included in haplotypes_to_track.

    Args:
        bam_path: Path to BAM file.
        region: Genomic region (start is 1-based inclusive).
        haplotypes_to_track: HP tag values to emit separate series for.
            Pass None to emit total only.
        include_observed_haplotypes: Emit separate series for each observed HP tag.
        bin_size: Resolution in base pairs.
        mapq_threshold: Minimum mapping quality.

    Returns:
        dict mapping haplotype value (-1=total, HP value otherwise) to x/y depth series.
    """
    region_len = region.end - region.start + 1
    n_bins = (region_len + bin_size - 1) // bin_size

    # Difference arrays: index = bin index; +1 at coverage start, -1 at end
    total_diff: list[int] = [0] * (n_bins + 1)
    hap_diffs: dict[CoverageHaplotype, list[int]] = {}
    if haplotypes_to_track:
        for hp in haplotypes_to_track:
            hap_diffs[hp] = [0] * (n_bins + 1)

    region_start_0 = region.start - 1  # convert to 0-based for pysam
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        for record in bam.fetch(region.chromosome, region_start_0, region.end):
            if record.is_secondary or record.mapping_quality < mapq_threshold:
                continue
            hp_tag: CoverageHaplotype | None = None
            with suppress(KeyError):
                hp_tag = record.get_tag("HP")
            if include_observed_haplotypes and hp_tag is not None and hp_tag not in hap_diffs:
                hap_diffs[hp_tag] = [0] * (n_bins + 1)

            for block_start, block_end in _merged_coverage_blocks(record, min_cigar_gap):
                cov_start = max(block_start, region_start_0)
                cov_end = min(block_end, region.end)
                if cov_start < cov_end:
                    bin_s = (cov_start - region_start_0) // bin_size
                    bin_e = (cov_end - region_start_0 - 1) // bin_size + 1
                    bin_s = max(0, min(bin_s, n_bins))
                    bin_e = max(0, min(bin_e, n_bins))
                    if bin_s < bin_e:
                        total_diff[bin_s] += 1
                        total_diff[bin_e] -= 1
                        hp_key = hp_tag if hp_tag is not None else 0
                        if hp_key in hap_diffs:
                            hap_diffs[hp_key][bin_s] += 1
                            hap_diffs[hp_key][bin_e] -= 1
        bam.close()
    except (FileNotFoundError, OSError, ValueError):
        pass

    # Prefix-sum to get depths; build x/y lists
    result: dict[CoverageHaplotype, tuple[list[int], list[int]]] = {}
    x_vals = list(range(region.start, region.end + 1, bin_size))

    depth = 0
    y_total = []
    for i in range(n_bins):
        depth += total_diff[i]
        y_total.append(max(0, depth))
    result[-1] = (x_vals[: len(y_total)], y_total)

    for hp, diff in hap_diffs.items():
        depth = 0
        y_hp = []
        for i in range(n_bins):
            depth += diff[i]
            y_hp.append(max(0, depth))
        result[hp] = (x_vals[: len(y_hp)], y_hp)

    return result


def order_alignments(segments_by_read):
    """
    Sort the split segments for each read by their fwd_read_start (ascending),
    and assign alignment_order 1..n so numbering is stable for the plot.
    Modifies the input dictionary in place.
    """
    for _read_name, segments in segments_by_read.items():
        segments.sort(key=lambda seg: seg.fwd_read_start)
        for i, seg in enumerate(segments, start=1):
            seg.alignment_order = i


def fetch_reference_sequence(reference_path: str, region: Region) -> str | None:
    """Return the reference sequence string for region, or None on any error.

    Uses 0-based half-open coordinates internally (pysam convention).
    Region coordinates are 1-based inclusive.
    """
    try:
        with pysam.FastaFile(reference_path) as fasta:
            if region.chromosome not in fasta.references:
                return None
            # region.start/end are 1-based inclusive; pysam.fetch is 0-based half-open
            seq = fasta.fetch(region.chromosome, region.start - 1, region.end)
            return seq.upper() if seq else None
    except Exception:
        return None
