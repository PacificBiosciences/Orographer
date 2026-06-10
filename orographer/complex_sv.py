"""
Two-phase complex SV workflow: BFS-based connected-locus discovery, then
combined rendering + coverage accumulation per BAM sample.
"""

import logging
import math
import re
import statistics
import time
from collections import defaultdict
from contextlib import suppress
from dataclasses import dataclass

import pysam

from .bam_parser import (
    MIN_MAPQ,
    AlignmentSummary,
    CoverageHaplotype,
    ExcludeRegions,
    FwdStrandReadSegment,
    SplitAlignmentSegment,
    _merged_coverage_blocks,
    convert_alignment_to_segment,
    extract_variants,
    fetch_alignment_summaries,
    get_fwd_read_split_segments,
    get_optional_int_aux_tag,
    order_alignments,
    parse_sa_aux_val,
)
from .gtf_parser import parse_annotation_file
from .utils import Region, chromosome_sort_key, parse_coordinate
from .vcf_parser import parse_vcf_file

logger = logging.getLogger(__name__)

# ---- Constants ---------------------------------------------------------------

DISCOVERY_MARGIN_FRACTION = 0.10
COVERAGE_BIN_SIZE = 10
MAX_DISCOVERY_INTERVALS = 500
MAX_CONNECTED_READS = 10_000
MAX_TOTAL_RENDER_BP = 25_000_000
MAX_BFS_ROUNDS = 50
COMPLEX_SV_MIN_INDEL_BP = 20  # reads with no SA tag and no I/D op exceeding this are coverage-only
DEFAULT_REGION_CONNECTION_THRESHOLD = 2
FINAL_REGION_MERGE_DISTANCE_BP = 200_000
INCLUSION_REASON_CHIMERIC = "Chimeric"
INCLUSION_REASON_LARGE_INSERTION = "INS"
INCLUSION_REASON_LARGE_DELETION = "DEL"
INCLUSION_REASON_LARGE_SOFTCLIP = "Large soft clip"


# ---- Interval free functions -------------------------------------------------


def expand_interval(
    chrom: str,
    start: int,
    end: int,
    fraction: float = DISCOVERY_MARGIN_FRACTION,
) -> tuple[str, int, int]:
    """Expand a 1-based inclusive interval by *fraction* on each side."""
    margin = math.ceil((end - start + 1) * fraction)
    return chrom, max(1, start - margin), end + margin


def merge_sorted_intervals(
    intervals: list[tuple[str, int, int]],
    max_gap_bp: int = 0,
) -> list[tuple[str, int, int]]:
    """Merge intervals across any chromosomes when their same-chrom gap is small enough.

    Groups by chromosome, sorts by start, merges in one pass per chromosome.
    Returns a chromosome-sorted list of non-overlapping intervals.
    """
    if not intervals:
        return []

    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in intervals:
        by_chrom.setdefault(chrom, []).append((start, end))

    merged: list[tuple[str, int, int]] = []
    for chrom in sorted(by_chrom, key=chromosome_sort_key):
        spans = sorted(by_chrom[chrom])
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start <= cur_end + max_gap_bp:
                cur_end = max(cur_end, end)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((chrom, cur_start, cur_end))

    return merged


def subtract_visited_spans(
    candidates: list[tuple[str, int, int]],
    visited_by_chrom: dict[str, list[tuple[int, int]]],
) -> list[tuple[str, int, int]]:
    """Return the portions of *candidates* not already covered by *visited_by_chrom*.

    visited_by_chrom maps chromosome → sorted list of non-overlapping (start, end).
    May split a candidate into sub-intervals around visited gaps.
    """
    import bisect

    result: list[tuple[str, int, int]] = []
    for chrom, cand_start, cand_end in candidates:
        visited = visited_by_chrom.get(chrom, [])
        if not visited:
            result.append((chrom, cand_start, cand_end))
            continue

        starts = [v[0] for v in visited]
        # Find first visited span that might overlap [cand_start, cand_end]
        idx = bisect.bisect_right(starts, cand_end)
        # Walk back to first span whose end exceeds cand_start
        while idx > 0 and visited[idx - 1][1] > cand_start:
            idx -= 1

        pos = cand_start
        for i in range(idx, len(visited)):
            v_start, v_end = visited[i]
            if v_start >= cand_end:
                break
            if v_start > pos:
                result.append((chrom, pos, v_start))
            pos = max(pos, v_end)
        if pos < cand_end:
            result.append((chrom, pos, cand_end))

    return result


def _add_visited_span(
    visited_by_chrom: dict[str, list[tuple[int, int]]],
    chrom: str,
    start: int,
    end: int,
) -> None:
    """Insert (start, end) into visited_by_chrom[chrom], merging in place."""
    if chrom not in visited_by_chrom:
        visited_by_chrom[chrom] = [(start, end)]
        return
    spans = visited_by_chrom[chrom] + [(start, end)]
    spans.sort()
    merged = [spans[0]]
    for s, e in spans[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    visited_by_chrom[chrom] = merged


def _cigar_ref_len(cigar_str: str) -> int:
    """Compute the reference-consuming length of a CIGAR string."""
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_str):
        if op in ("M", "D", "N", "=", "X"):
            total += int(length)
    return total


def inclusion_reason_for_record(record) -> str:
    """Return why a BAM record is rendered in the complex-SV read track."""
    if record.is_supplementary or record.has_tag("SA"):
        return INCLUSION_REASON_CHIMERIC

    indel_reason = large_indel_inclusion_reason(record)
    if indel_reason:
        return indel_reason

    cigar = record.cigartuples or []
    max_softclip = max((length for op, length in cigar if op == 4), default=0)
    if max_softclip > COMPLEX_SV_MIN_INDEL_BP:
        return INCLUSION_REASON_LARGE_SOFTCLIP

    return ""


def large_indel_inclusion_reason(record) -> str:
    """Return a specific large insertion/deletion reason with its 1-based POS."""
    cigar = record.cigartuples or []
    ref_pos = record.reference_start
    if ref_pos is None:
        return ""

    chrom = record.reference_name or "unknown"
    best_reason = ""
    best_length = 0
    for op, length in cigar:
        if op == 1 and length > COMPLEX_SV_MIN_INDEL_BP and length > best_length:
            ins_pos = max(record.reference_start + 1, ref_pos)
            best_length = length
            best_reason = f"{length:,} bp {INCLUSION_REASON_LARGE_INSERTION} at {chrom}:{ins_pos}"
        elif op == 2:
            if length > COMPLEX_SV_MIN_INDEL_BP and length > best_length:
                best_length = length
                best_reason = (
                    f"{length:,} bp {INCLUSION_REASON_LARGE_DELETION} at {chrom}:{ref_pos + 1}"
                )

        if op in (0, 2, 3, 7, 8):
            ref_pos += length

    return best_reason


@dataclass
class SATargetObservation:
    """One read-supported SA target interval before target-locus merging."""

    read_name: str
    chrom: str
    start: int
    end: int


@dataclass
class SATargetLocus:
    """Merged SA target locus with unique read support."""

    chrom: str
    start: int
    end: int
    read_names: set[str]


def sa_target_observations_for_summary(
    summary: AlignmentSummary,
    mapq_threshold: int,
    sa_cache: dict[str, list[SplitAlignmentSegment]],
) -> list[SATargetObservation]:
    """Return qualifying SA target observations for one alignment summary."""
    if not summary.has_sa or not summary.sa_tag:
        return []
    if summary.sa_tag not in sa_cache:
        sa_cache[summary.sa_tag] = parse_sa_aux_val(summary.sa_tag)

    observations = []
    for sa_seg in sa_cache[summary.sa_tag]:
        if sa_seg.mapq < mapq_threshold:
            continue
        ref_len = _cigar_ref_len(sa_seg.cigar)
        # sa_seg.pos is 0-based, so convert to 1-based inclusive coordinates.
        sa_start = sa_seg.pos + 1
        sa_end = sa_seg.pos + max(ref_len, 1)
        observations.append(
            SATargetObservation(
                read_name=summary.read_name,
                chrom=sa_seg.rname,
                start=sa_start,
                end=sa_end,
            )
        )
    return observations


def merge_sa_target_observations(
    observations: list[SATargetObservation],
) -> list[SATargetLocus]:
    """Expand and merge nearby SA targets, preserving unique read support."""
    if not observations:
        return []

    expanded = []
    for observation in observations:
        chrom, start, end = expand_interval(
            observation.chrom,
            observation.start,
            observation.end,
        )
        expanded.append((chrom, start, end, observation.read_name))

    merged: list[SATargetLocus] = []

    def sort_key(item: tuple[str, int, int, str]) -> tuple[tuple[int, int, str], int, int]:
        return chromosome_sort_key(item[0]), item[1], item[2]

    for chrom, start, end, read_name in sorted(expanded, key=sort_key):
        if merged and merged[-1].chrom == chrom and start <= merged[-1].end:
            merged[-1].end = max(merged[-1].end, end)
            merged[-1].read_names.add(read_name)
            continue
        merged.append(SATargetLocus(chrom, start, end, {read_name}))
    return merged


def filter_supported_sa_target_loci(
    observations: list[SATargetObservation],
    region_connection_threshold: int,
) -> list[tuple[str, int, int]]:
    """Return merged SA target loci meeting the unique-read support threshold."""
    target_loci = merge_sa_target_observations(observations)
    retained = [
        (locus.chrom, locus.start, locus.end)
        for locus in target_loci
        if len(locus.read_names) >= region_connection_threshold
    ]
    filtered_count = len(target_loci) - len(retained)
    logger.debug(
        "SA target filter: %d observation(s), %d retained locus/loci, %d filtered.",
        len(observations),
        len(retained),
        filtered_count,
    )
    for locus in target_loci:
        support_count = len(locus.read_names)
        if support_count >= region_connection_threshold:
            continue
        logger.debug(
            "Filtered SA target %s:%d-%d with %d unique supporting read(s); threshold is %d.",
            locus.chrom,
            locus.start,
            locus.end,
            support_count,
            region_connection_threshold,
        )
    return retained


# ---- BFS discovery -----------------------------------------------------------


def discover_connected_reads(
    bam_file: pysam.AlignmentFile,
    initial_frontier: list[tuple[str, int, int]],
    mapq_threshold: int = MIN_MAPQ,
    region_connection_threshold: int = DEFAULT_REGION_CONNECTION_THRESHOLD,
) -> tuple[set[str], dict[str, list[AlignmentSummary]], list[tuple[str, int, int]]]:
    """Round-based BFS: follow SA-tag links from *initial_frontier*.

    Returns:
        discovered_reads:   all read names encountered
        summaries_by_read:  AlignmentSummary lists keyed by read name
        visited_intervals:  every 1-based interval that was fetched (already expanded)
    """
    if region_connection_threshold < 1:
        raise ValueError("region_connection_threshold must be at least 1")

    visited_by_chrom: dict[str, list[tuple[int, int]]] = {}
    discovered_reads: set[str] = set()
    expanded_reads: set[str] = set()
    summaries_by_read: dict[str, list[AlignmentSummary]] = defaultdict(list)
    sa_cache: dict[str, list[SplitAlignmentSegment]] = {}
    all_visited: list[tuple[str, int, int]] = []

    frontier: list[tuple[str, int, int]] = list(initial_frontier)

    for round_num in range(MAX_BFS_ROUNDS):
        if not frontier:
            logger.debug("BFS complete after %d round(s).", round_num)
            break

        # Expand, merge, subtract visited
        expanded = [expand_interval(*iv) for iv in frontier]
        merged = merge_sorted_intervals(expanded)
        new_intervals = subtract_visited_spans(merged, visited_by_chrom)

        if not new_intervals:
            logger.debug("BFS round %d: no new intervals after subtraction.", round_num)
            break

        # Guardrail: interval count
        remaining = MAX_DISCOVERY_INTERVALS - len(all_visited)
        if len(new_intervals) > remaining:
            logger.warning(
                "Discovery interval limit (%d) reached at round %d. Truncating BFS.",
                MAX_DISCOVERY_INTERVALS,
                round_num,
            )
            new_intervals = new_intervals[:remaining]

        logger.debug(
            "BFS round %d: fetching %d interval(s) across %d chromosome(s).",
            round_num,
            len(new_intervals),
            len({c for c, _, _ in new_intervals}),
        )

        sa_target_observations: list[SATargetObservation] = []

        for chrom, start, end in new_intervals:
            _add_visited_span(visited_by_chrom, chrom, start, end)
            all_visited.append((chrom, start, end))

            region = Region(chrom, start, end, "")
            for summary in fetch_alignment_summaries(bam_file, region, mapq_threshold):
                qualifies = (
                    summary.has_sa
                    or summary.max_indel > COMPLEX_SV_MIN_INDEL_BP
                    or summary.max_softclip > COMPLEX_SV_MIN_INDEL_BP
                )
                if qualifies:
                    summaries_by_read[summary.read_name].append(summary)
                    discovered_reads.add(summary.read_name)

                if summary.has_sa and summary.read_name not in expanded_reads:
                    expanded_reads.add(summary.read_name)
                    sa_target_observations.extend(
                        sa_target_observations_for_summary(
                            summary,
                            mapq_threshold,
                            sa_cache,
                        )
                    )

        if len(discovered_reads) >= MAX_CONNECTED_READS:
            logger.warning(
                "Connected read limit (%d) reached at round %d. Truncating BFS.",
                MAX_CONNECTED_READS,
                round_num,
            )
            break

        frontier = filter_supported_sa_target_loci(
            sa_target_observations,
            region_connection_threshold,
        )
    else:
        logger.warning("BFS reached max rounds (%d). Stopping early.", MAX_BFS_ROUNDS)

    logger.debug(
        "BFS finished: %d read(s), %d interval(s) visited.",
        len(discovered_reads),
        len(all_visited),
    )
    return discovered_reads, dict(summaries_by_read), all_visited


def discover_reads_in_regions(
    bam_file: pysam.AlignmentFile,
    regions: list[tuple[str, int, int]],
    mapq_threshold: int = MIN_MAPQ,
) -> tuple[set[str], dict[str, list[AlignmentSummary]], list[tuple[str, int, int]]]:
    """Discover qualifying reads inside explicit regions without following SA targets."""
    discovered_reads: set[str] = set()
    summaries_by_read: dict[str, list[AlignmentSummary]] = defaultdict(list)

    for chrom, start, end in regions:
        region = Region(chrom, start, end, "")
        for summary in fetch_alignment_summaries(bam_file, region, mapq_threshold):
            qualifies = (
                summary.has_sa
                or summary.max_indel > COMPLEX_SV_MIN_INDEL_BP
                or summary.max_softclip > COMPLEX_SV_MIN_INDEL_BP
            )
            if not qualifies:
                continue
            summaries_by_read[summary.read_name].append(summary)
            discovered_reads.add(summary.read_name)

    logger.debug(
        "Explicit-region discovery finished: %d read(s), %d region(s) visited.",
        len(discovered_reads),
        len(regions),
    )
    return discovered_reads, dict(summaries_by_read), list(regions)


# ---- Final region building ---------------------------------------------------


def _build_final_regions(
    summaries_by_read: dict[str, list[AlignmentSummary]],
    visited_intervals: list[tuple[str, int, int]],
) -> list[Region]:
    """Derive final visualization Regions from alignment extents of discovered reads.

    Uses only alignment extents from summaries_by_read so final regions correspond
    to actual read positions, not intermediate BFS fetch intervals (which can produce
    empty-looking panels with no connected-read segments). Falls back to
    visited_intervals only when no summaries are available.
    """
    candidates: list[tuple[str, int, int]] = []
    for summaries in summaries_by_read.values():
        for s in summaries:
            start_1 = s.start + 1  # 0-based → 1-based inclusive
            end_1 = s.end  # 0-based exclusive = 1-based inclusive
            if start_1 <= end_1:
                candidates.append((s.chromosome, start_1, end_1))

    if not candidates:
        candidates = list(visited_intervals)

    expanded = [expand_interval(*iv) for iv in candidates]
    merged = merge_sorted_intervals(
        expanded,
        max_gap_bp=FINAL_REGION_MERGE_DISTANCE_BP,
    )
    # Sort by natural chromosome order, then start for deterministic display order.
    merged.sort(key=lambda t: (chromosome_sort_key(t[0]), t[1]))
    return [Region(chrom, start, end, "") for chrom, start, end in merged]


def _truncate_regions(regions: list[Region], max_bp: int) -> list[Region]:
    """Drop trailing regions once the cumulative span would exceed *max_bp*."""
    kept: list[Region] = []
    total = 0
    for r in regions:
        span = r.end - r.start + 1
        if total + span > max_bp:
            logger.warning(
                "Dropping region %s:%d-%d: total render span limit reached.",
                r.chromosome,
                r.start,
                r.end,
            )
            continue
        kept.append(r)
        total += span
    return kept


# ---- Combined rendering + coverage pass --------------------------------------
@dataclass
class CoverageState:
    """Mutable binned coverage diff arrays for one rendered region."""

    region: Region
    bin_size: int
    total_diff: list[int]
    hap_diffs: dict[CoverageHaplotype, list[int]]
    include_observed_haplotypes: bool = False


@dataclass
class SegmentRenderState:
    """Mutable connected-read segment extraction state for one BAM and region."""

    segments_by_read: dict[str, list[FwdStrandReadSegment]]
    processed_alignments: set[tuple[str, int]]
    exclude_regions: ExcludeRegions


def create_coverage_state(
    region: Region,
    haplotypes_to_track: set[CoverageHaplotype] | None,
    bin_size: int,
    include_observed_haplotypes: bool = False,
) -> CoverageState:
    """Create empty diff arrays for binned coverage accumulation."""
    region_len = region.end - region.start + 1
    n_bins = (region_len + bin_size - 1) // bin_size
    hap_diffs: dict[CoverageHaplotype, list[int]] = {}
    if haplotypes_to_track:
        for hp in haplotypes_to_track:
            hap_diffs[hp] = [0] * (n_bins + 1)
    return CoverageState(
        region=region,
        bin_size=bin_size,
        total_diff=[0] * (n_bins + 1),
        hap_diffs=hap_diffs,
        include_observed_haplotypes=include_observed_haplotypes,
    )


def add_record_coverage(record, coverage_state: CoverageState, min_cigar_gap: int) -> None:
    """Accumulate one BAM record into total and requested haplotype coverage diffs."""
    region = coverage_state.region
    region_start_0 = region.start - 1
    n_bins = len(coverage_state.total_diff) - 1
    hp_tag: CoverageHaplotype | None = None
    with suppress(KeyError):
        hp_tag = record.get_tag("HP")
    if (
        coverage_state.include_observed_haplotypes
        and hp_tag is not None
        and hp_tag not in coverage_state.hap_diffs
    ):
        coverage_state.hap_diffs[hp_tag] = [0] * (n_bins + 1)

    for block_start, block_end in _merged_coverage_blocks(record, min_cigar_gap):
        cov_start = max(block_start, region_start_0)
        cov_end = min(block_end, region.end)
        if cov_start >= cov_end:
            continue

        bin_s = (cov_start - region_start_0) // coverage_state.bin_size
        bin_e = (cov_end - region_start_0 - 1) // coverage_state.bin_size + 1
        bin_s = max(0, min(bin_s, n_bins))
        bin_e = max(0, min(bin_e, n_bins))
        if bin_s >= bin_e:
            continue

        coverage_state.total_diff[bin_s] += 1
        coverage_state.total_diff[bin_e] -= 1
        hp_key = hp_tag if hp_tag is not None else 0
        if hp_key in coverage_state.hap_diffs:
            coverage_state.hap_diffs[hp_key][bin_s] += 1
            coverage_state.hap_diffs[hp_key][bin_e] -= 1


def build_coverage_tracks(
    coverage_state: CoverageState,
) -> dict[CoverageHaplotype, tuple[list[int], list[int]]]:
    """Convert accumulated coverage diffs into plottable depth tracks."""
    region = coverage_state.region
    x_vals = list(range(region.start, region.end + 1, coverage_state.bin_size))
    coverage_tracks: dict[CoverageHaplotype, tuple[list[int], list[int]]] = {}

    depth = 0
    y_total = []
    for diff_value in coverage_state.total_diff[:-1]:
        depth += diff_value
        y_total.append(max(0, depth))
    coverage_tracks[-1] = (x_vals[: len(y_total)], y_total)

    for hp, diff in coverage_state.hap_diffs.items():
        depth = 0
        y_hp = []
        for diff_value in diff[:-1]:
            depth += diff_value
            y_hp.append(max(0, depth))
        coverage_tracks[hp] = (x_vals[: len(y_hp)], y_hp)

    return coverage_tracks


def create_segment_render_state() -> SegmentRenderState:
    """Create empty state for extracting connected-read segments."""
    return SegmentRenderState(
        segments_by_read=defaultdict(list),
        processed_alignments=set(),
        exclude_regions=ExcludeRegions(regions={}),
    )


def add_connected_read_segments(
    record,
    segment_state: SegmentRenderState,
    region: Region,
    ref_fasta,
) -> None:
    """Extract or enrich connected-read segments from one qualifying BAM record."""
    read_name = record.query_name
    inclusion_reason = inclusion_reason_for_record(record)
    aln_key = (read_name, record.reference_start)
    if aln_key in segment_state.processed_alignments:
        return
    segment_state.processed_alignments.add(aln_key)

    if record.is_supplementary:
        add_supplementary_record_segments(record, segment_state, ref_fasta, inclusion_reason)
        return

    if record.has_tag("SA"):
        add_split_read_segments(record, segment_state, region, ref_fasta, inclusion_reason)
        return

    seg = convert_alignment_to_segment(record, region.chromosome, ref_fasta)
    if seg:
        seg.inclusion_reason = inclusion_reason
        segment_state.segments_by_read[read_name].append(seg)


def add_supplementary_record_segments(
    record,
    segment_state: SegmentRenderState,
    ref_fasta,
    inclusion_reason: str,
) -> None:
    """Add a supplementary record, or enrich the matching SA-derived segment."""
    read_name = record.query_name
    supp_hp = get_optional_int_aux_tag(record, "HP")
    supp_ps = get_optional_int_aux_tag(record, "PS")
    existing = [
        seg
        for seg in segment_state.segments_by_read[read_name]
        if seg.pos == record.reference_start
    ]
    if existing:
        for seg in existing:
            seg.haplotype_tag = supp_hp
            seg.phaseset_tag = supp_ps
            seg.inclusion_reason = inclusion_reason
        if ref_fasta:
            mismatches, insertions, deletions = extract_variants(record, ref_fasta)
            for seg in existing:
                seg.mismatches = list(mismatches)
                seg.insertions = list(insertions)
                seg.deletions = list(deletions)
        return

    seg = convert_alignment_to_segment(record, record.reference_name, ref_fasta)
    if seg:
        seg.from_primary_bam_record = False
        seg.inclusion_reason = inclusion_reason
        segment_state.segments_by_read[read_name].append(seg)


def add_split_read_segments(
    record,
    segment_state: SegmentRenderState,
    region: Region,
    ref_fasta,
    inclusion_reason: str,
) -> None:
    """Add primary and SA-derived split-read segments from one primary record."""
    split_segs = get_fwd_read_split_segments(
        record,
        region.chromosome,
        segment_state.exclude_regions,
        ref_fasta,
    )
    if not split_segs:
        return

    existing_pos = {seg.pos for seg in segment_state.segments_by_read[record.query_name]}
    for seg in split_segs:
        seg.inclusion_reason = inclusion_reason
        if seg.from_primary_bam_record:
            segment_state.segments_by_read[record.query_name].append(seg)
        elif seg.pos not in existing_pos:
            segment_state.segments_by_read[record.query_name].append(seg)
            existing_pos.add(seg.pos)


def group_segments_by_haplotype(
    segments_by_read: dict[str, list[FwdStrandReadSegment]],
) -> dict[str, list[FwdStrandReadSegment]]:
    """Group ordered read segments into read or read_HP keys used by plotting."""
    order_alignments(segments_by_read)
    segments_by_haplotype: dict[str, list[FwdStrandReadSegment]] = {}
    for read_name, segs in segments_by_read.items():
        for seg in segs:
            hp = seg.haplotype_tag if seg.haplotype_tag is not None else 0
            hkey = read_name if hp == 0 else f"{read_name}_HP{hp}"
            segments_by_haplotype.setdefault(hkey, []).append(seg)
    for hkey in segments_by_haplotype:
        segments_by_haplotype[hkey].sort(key=lambda seg: seg.alignment_order)
    return segments_by_haplotype


def collect_insertion_summary(
    segments_by_read: dict[str, list[FwdStrandReadSegment]],
    region: Region,
) -> dict[int, list[dict]]:
    """Summarize large discovered-read insertions by haplotype and reference site."""
    insertion_by_hp: dict[int, dict[int, dict]] = {}
    for read_name, segs in segments_by_read.items():
        seen: set[tuple[int, int]] = set()
        for seg in segs:
            if not seg.insertions:
                continue
            hp = seg.haplotype_tag if seg.haplotype_tag is not None else 0
            for ref_pos_0, inserted_bases in seg.insertions:
                if len(inserted_bases) <= COMPLEX_SV_MIN_INDEL_BP:
                    continue
                pos_1 = ref_pos_0 + 1
                if not (region.start <= pos_1 <= region.end):
                    continue
                key = (hp, pos_1)
                if key in seen:
                    continue
                seen.add(key)
                site = insertion_by_hp.setdefault(hp, {}).setdefault(
                    pos_1, {"sizes": [], "read_names": []}
                )
                site["sizes"].append(len(inserted_bases))
                site["read_names"].append(read_name)

    insertion_summary: dict[int, list[dict]] = {}
    for hp, sites in insertion_by_hp.items():
        summary = []
        for pos, data in sorted(sites.items()):
            all_names = data["read_names"]
            summary.append(
                {
                    "pos": pos,
                    "count": len(all_names),
                    "median_size": round(statistics.median(data["sizes"])),
                    "read_names": all_names[:5],
                    "total_count": len(all_names),
                    "chrom": region.chromosome,
                }
            )
        insertion_summary[hp] = summary
    return insertion_summary


def render_region_for_bam(
    bam_path: str,
    region: Region,
    discovered_reads: set[str],
    reference_path: str | None,
    haplotypes_to_track: set[CoverageHaplotype] | None,
    bin_size: int = COVERAGE_BIN_SIZE,
    mapq_threshold: int = MIN_MAPQ,
) -> tuple[
    dict[str, list[FwdStrandReadSegment]],
    dict[CoverageHaplotype, tuple[list[int], list[int]]],
    dict[int, list[dict]],
]:
    """Single BAM pass: extract read segments for connected reads and accumulate coverage.

    Returns:
        segments_by_haplotype_key: same shape as fetch_all_alignments output
        coverage_tracks: haplotype value → (x_positions, y_depths)
        insertion_summary: haplotype int → list of site dicts
            each dict: {pos, count, median_size, read_names (first 5), total_count, chrom}
    """
    region_start_0 = region.start - 1  # 0-based for pysam
    coverage_state = create_coverage_state(
        region,
        haplotypes_to_track,
        bin_size,
        include_observed_haplotypes=True,
    )
    segment_state = create_segment_render_state()

    ref_fasta = None
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        if reference_path:
            ref_fasta = pysam.FastaFile(reference_path)

        for record in bam.fetch(region.chromosome, region_start_0, region.end):
            if record.is_secondary:
                continue
            if record.mapping_quality < mapq_threshold:
                continue

            add_record_coverage(record, coverage_state, COMPLEX_SV_MIN_INDEL_BP)

            if record.query_name not in discovered_reads:
                continue

            add_connected_read_segments(record, segment_state, region, ref_fasta)

        bam.close()
        if ref_fasta:
            ref_fasta.close()

    except (FileNotFoundError, OSError, ValueError) as exc:
        logger.error("Error reading BAM %s: %s", bam_path, exc)
        return {}, {}, {}

    segments_by_haplotype = group_segments_by_haplotype(segment_state.segments_by_read)
    coverage_tracks = build_coverage_tracks(coverage_state)
    insertion_summary = collect_insertion_summary(segment_state.segments_by_read, region)
    return dict(segments_by_haplotype), coverage_tracks, insertion_summary


# ---- Top-level orchestrator --------------------------------------------------


def process_complex_sv(
    coordinate_strs: list[str],
    bam_files: list[str],
    vcf_files: list[str | None],
    reference_path: str | None,
    gtf_file: str | None,
    sample_labels: list[str | None],
    region_type: str,
    region_connection_threshold: int = DEFAULT_REGION_CONNECTION_THRESHOLD,
    no_expansion: bool = False,
) -> list[dict]:
    """Two-phase complex SV workflow: discovery then render+coverage.

    Discovery runs once against the primary BAM. All other BAMs are rendered
    against the shared set of final visualization regions.

    Returns region_data entries (one per final discovered region) in the same
    shape expected by plot_reads_bokeh:
        {"region": Region, "gene_annotations": [...], "bam_rows": [...]}
    """
    t0 = time.monotonic()

    # Parse seeds (1-based inclusive)
    seeds: list[tuple[str, int, int]] = [parse_coordinate(c) for c in coordinate_strs]

    if no_expansion:
        logger.debug("Complex SV: expansion disabled for %d explicit seed region(s).", len(seeds))
    else:
        # Unified initial BFS frontier: expand all seeds, merge
        initial_frontier = merge_sorted_intervals(
            [expand_interval(chrom, start, end) for chrom, start, end in seeds]
        )
        logger.debug(
            "Complex SV: %d seed(s) → %d initial frontier interval(s)",
            len(seeds),
            len(initial_frontier),
        )

    # Phase 1: BFS discovery with primary BAM
    primary_bam_path = bam_files[0]
    try:
        primary_bam = pysam.AlignmentFile(primary_bam_path, "rb")
    except (FileNotFoundError, OSError) as exc:
        raise RuntimeError(f"Cannot open primary BAM: {primary_bam_path}") from exc

    try:
        if no_expansion:
            discovered_reads, summaries_by_read, visited_intervals = discover_reads_in_regions(
                primary_bam,
                seeds,
                mapq_threshold=0,
            )
        else:
            discovered_reads, summaries_by_read, visited_intervals = discover_connected_reads(
                primary_bam,
                initial_frontier,
                mapq_threshold=0,
                region_connection_threshold=region_connection_threshold,
            )
    finally:
        primary_bam.close()

    t_discovery = time.monotonic() - t0
    logger.debug(
        "Discovery: %.2fs — %d read(s), %d interval(s) visited",
        t_discovery,
        len(discovered_reads),
        len(visited_intervals),
    )

    if no_expansion:
        final_regions = [
            Region(chrom, start, end, f"{chrom}:{start}-{end}") for chrom, start, end in seeds
        ]
    elif not discovered_reads:
        logger.warning("No connected reads found; falling back to seed regions.")
        final_regions: list[Region] = [
            Region(chrom, start, end, f"{chrom}:{start}-{end}") for chrom, start, end in seeds
        ]
    else:
        final_regions = _build_final_regions(summaries_by_read, visited_intervals)

    logger.debug("Final visualization regions: %d", len(final_regions))

    total_bp = sum(r.end - r.start + 1 for r in final_regions)
    if total_bp > MAX_TOTAL_RENDER_BP:
        logger.warning(
            "Total render span %d bp exceeds limit %d. Truncating.",
            total_bp,
            MAX_TOTAL_RENDER_BP,
        )
        final_regions = _truncate_regions(final_regions, MAX_TOTAL_RENDER_BP)

    # Standard PacBio HP values are reserved; named HP tags are added as observed.
    haplotypes_to_track: set[CoverageHaplotype] = {0, 1, 2}

    # Phase 2: Render + coverage for each final region, each BAM
    region_data_list: list[dict] = []

    for region in final_regions:
        gene_annotations = parse_annotation_file(gtf_file, region)
        rows_in_bam_order: list[dict] = []

        for bam_index, bam_path in enumerate(bam_files):
            vcf_path = vcf_files[bam_index] if bam_index < len(vcf_files) else None
            vcf_variants = parse_vcf_file(vcf_path, region)

            segments_by_read, coverage_tracks, insertion_summary = render_region_for_bam(
                bam_path,
                region,
                discovered_reads,
                reference_path,
                haplotypes_to_track,
                mapq_threshold=0,
            )

            rows_in_bam_order.append(
                {
                    "segments_by_read": segments_by_read,
                    "vcf_variants": vcf_variants,
                    "coverage_tracks": coverage_tracks,
                    "insertion_summary": insertion_summary,
                    "region_type": region_type,
                    "sample_label": (
                        sample_labels[bam_index] if bam_index < len(sample_labels) else None
                    ),
                    "sample_index": bam_index,
                }
            )

        # Display order: [other1, other2, primary] (primary at bottom)
        if len(rows_in_bam_order) == 3:
            bam_rows = [rows_in_bam_order[1], rows_in_bam_order[2], rows_in_bam_order[0]]
        elif len(rows_in_bam_order) == 2:
            bam_rows = [rows_in_bam_order[1], rows_in_bam_order[0]]
        else:
            bam_rows = [rows_in_bam_order[0]]

        region_data_list.append(
            {
                "region": region,
                "gene_annotations": gene_annotations,
                "bam_rows": bam_rows,
            }
        )

    t_render = time.monotonic() - t0 - t_discovery
    logger.debug(
        "Rendering: %.2fs across %d region(s) for %d BAM(s)",
        t_render,
        len(final_regions),
        len(bam_files),
    )

    return region_data_list
