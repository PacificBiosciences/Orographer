"""IsoSeq transcript-evidence workflow and assignment helpers."""

from __future__ import annotations

import hashlib
import logging
import os
import pickle
import tempfile
import time
from contextlib import suppress
from dataclasses import dataclass

import pysam

from .bam_parser import (
    MIN_MAPQ,
    MIN_SOFTCLIPS,
    get_complete_read_clip_positions,
    get_fwd_strand_read_pos,
    get_optional_int_aux_tag,
    get_optional_string_aux_tag,
)
from .gtf_parser import TranscriptAnnotation, parse_transcript_annotation_file
from .utils import Region, parse_coordinate
from .vcf_parser import parse_vcf_file

logger = logging.getLogger(__name__)
ISOSEQ_TRANSCRIPT_CACHE_VERSION = "v1"
CIGAR_DEL = 2
CIGAR_REF_SKIP = 3
CIGAR_MATCH_OPS = {0, 7, 8}


@dataclass(frozen=True, slots=True)
class IsoSeqReadModel:
    """Splice-chain and endpoint summary for transcript assignment."""

    read_name: str
    start: int
    end: int
    intron_chain: tuple[tuple[int, int], ...]


@dataclass(frozen=True, slots=True)
class IsoSeqAssignment:
    """Read assignment result."""

    transcript_id: str | None
    gene_id: str | None
    gene_name: str | None


@dataclass(slots=True)
class IsoSeqReadSegment:
    """Minimal read segment needed for IsoSeq read-page chunk generation."""

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
    color_tag: str | None
    spans: bool
    from_primary_bam_record: bool
    readname: str
    aligned_blocks: list[tuple[int, int]]
    alignment_order: int = 1


def _pg_texts(bam_path: str) -> list[str]:
    """Return lower-cased @PG records flattened into searchable strings."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        pg_records = bam.header.to_dict().get("PG", [])
    texts = []
    for record in pg_records:
        parts = [str(value) for value in record.values() if value is not None]
        texts.append(" ".join(parts).lower())
    return texts


def missing_isoseq_header_requirements(bam_path: str) -> list[str]:
    """Return missing IsoSeq provenance requirements from a BAM header."""
    texts = _pg_texts(bam_path)
    joined = "\n".join(texts)

    has_lima = "lima" in joined
    has_refine = any("isoseq" in text and "refine" in text for text in texts)
    has_pbmm2_isoseq = any(
        "pbmm2" in text and "isoseq" in text and ("preset" in text or "--preset" in text)
        for text in texts
    )

    missing = []
    if not has_lima:
        missing.append("lima")
    if not has_refine:
        missing.append("isoseq refine")
    if not has_pbmm2_isoseq:
        missing.append("pbmm2 align --preset ISOSEQ")
    return missing


def validate_isoseq_bam_header(bam_path: str, *, force: bool = False) -> None:
    """Validate IsoSeq pipeline provenance in BAM @PG headers."""
    missing = missing_isoseq_header_requirements(bam_path)
    if not missing:
        return

    message = f"IsoSeq input validation failed for {bam_path}: missing {', '.join(missing)}"
    if force:
        logger.warning("%s; continuing because --force-isoseq was provided.", message)
        return
    raise ValueError(message)


def _read_model_from_record(record: pysam.AlignedSegment) -> IsoSeqReadModel | None:
    """Build an assignment model from one primary alignment record."""
    if record.is_unmapped or record.is_secondary or record.is_supplementary:
        return None
    blocks = _reference_blocks_split_on_refskip(
        record.cigartuples or [],
        record.reference_start,
    )
    if not blocks:
        return None
    introns = tuple((blocks[i][1], blocks[i + 1][0] + 1) for i in range(len(blocks) - 1))
    return IsoSeqReadModel(
        read_name=record.query_name,
        start=blocks[0][0] + 1,
        end=blocks[-1][1],
        intron_chain=introns,
    )


def collect_read_models(bam_path: str, region: Region) -> dict[str, IsoSeqReadModel]:
    """Collect read splice-chain models from primary alignments in a region."""
    models = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam.fetch(region.chromosome, region.start - 1, region.end):
            model = _read_model_from_record(record)
            if model is None:
                continue
            previous = models.get(model.read_name)
            if previous is None:
                models[model.read_name] = model
                continue
            # Keep the primary alignment with the largest visible span if duplicates exist.
            if (model.end - model.start) > (previous.end - previous.start):
                models[model.read_name] = model
    return models


def _reference_blocks_split_on_refskip(
    cigartuples: list[tuple[int, int]],
    reference_start: int,
) -> list[tuple[int, int]]:
    """Return display blocks split only by CIGAR reference skips."""
    blocks = []
    block_start = None
    ref_pos = reference_start
    for op, length in cigartuples:
        if op in CIGAR_MATCH_OPS or op == CIGAR_DEL:
            if block_start is None:
                block_start = ref_pos
            ref_pos += length
            continue
        if op == CIGAR_REF_SKIP:
            if block_start is not None and block_start != ref_pos:
                blocks.append((block_start, ref_pos))
            ref_pos += length
            block_start = None
    if block_start is not None and block_start != ref_pos:
        blocks.append((block_start, ref_pos))
    return blocks


def _display_read_name(read_name: str, haplotype_tag: int | None) -> str:
    """Return the same haplotype-suffixed display name used by read plotting."""
    if haplotype_tag is None or haplotype_tag == 0:
        return read_name
    return f"{read_name}_HP{haplotype_tag}"


def _softclip_status(cigartuples: list[tuple[int, int]]) -> tuple[bool, bool]:
    if not cigartuples:
        return False, False
    first_op, first_length = cigartuples[0]
    last_op, last_length = cigartuples[-1]
    return (
        first_op == 4 and first_length > MIN_SOFTCLIPS,
        last_op == 4 and last_length > MIN_SOFTCLIPS,
    )


def _isoseq_segment_from_record(record: pysam.AlignedSegment) -> IsoSeqReadSegment | None:
    """Create a minimal plotting segment from one primary IsoSeq alignment."""
    if record.is_unmapped or record.is_secondary or record.is_supplementary:
        return None
    if record.mapping_quality < MIN_MAPQ:
        return None
    cigartuples = list(record.cigartuples or [])
    blocks = _reference_blocks_split_on_refskip(cigartuples, record.reference_start)
    if not blocks:
        return None
    read_start, read_end, read_size = get_complete_read_clip_positions(cigartuples)
    fwd_read_start, fwd_read_end = get_fwd_strand_read_pos(
        read_start,
        read_end,
        read_size,
        not record.is_reverse,
    )
    is_start_softclipped, is_end_softclipped = _softclip_status(cigartuples)
    return IsoSeqReadSegment(
        fwd_read_start=fwd_read_start,
        fwd_read_end=fwd_read_end,
        chrom=record.reference_name,
        second_chrom=record.reference_name,
        pos=record.reference_start,
        end=record.reference_end,
        is_fwd_strand=not record.is_reverse,
        is_start_softclipped=is_start_softclipped,
        is_end_softclipped=is_end_softclipped,
        phaseset_tag=get_optional_int_aux_tag(record, "PS"),
        haplotype_tag=get_optional_int_aux_tag(record, "HP"),
        color_tag=get_optional_string_aux_tag(record, "YC"),
        spans=True,
        from_primary_bam_record=True,
        readname=record.query_name,
        aligned_blocks=blocks,
    )


def collect_isoseq_segments_and_models(
    bam_path: str,
    region: Region,
) -> tuple[dict[str, list[IsoSeqReadSegment]], dict[str, IsoSeqReadModel]]:
    """Collect lightweight IsoSeq display segments and assignment models in one BAM pass."""
    segments_by_read: dict[str, list[IsoSeqReadSegment]] = {}
    models: dict[str, IsoSeqReadModel] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam.fetch(region.chromosome, region.start - 1, region.end):
            segment = _isoseq_segment_from_record(record)
            if segment is None:
                continue
            display_name = _display_read_name(segment.readname, segment.haplotype_tag)
            previous_segments = segments_by_read.get(display_name)
            if previous_segments:
                previous = previous_segments[0]
                if (segment.end - segment.pos) <= (previous.end - previous.pos):
                    continue
            segments_by_read[display_name] = [segment]

            model = _read_model_from_record(record)
            if model is None:
                continue
            previous_model = models.get(model.read_name)
            if previous_model is None:
                models[model.read_name] = model
                continue
            if (model.end - model.start) > (previous_model.end - previous_model.start):
                models[model.read_name] = model
    return segments_by_read, models


def _assignment_sort_key(
    transcript: TranscriptAnnotation, read_model: IsoSeqReadModel
) -> tuple[int, str, str, str, int]:
    endpoint_score = abs(read_model.start - transcript.start) + abs(read_model.end - transcript.end)
    return (
        endpoint_score,
        transcript.gene_name,
        transcript.gene_id,
        transcript.transcript_id,
        transcript.start,
    )


def _read_utrs_within_transcript(
    read_model: IsoSeqReadModel,
    transcript: TranscriptAnnotation,
) -> bool:
    return transcript.start <= read_model.start and read_model.end <= transcript.end


def assign_read_to_transcript(
    read_model: IsoSeqReadModel, transcripts: list[TranscriptAnnotation]
) -> TranscriptAnnotation | None:
    """Assign one read to the best full-splice-matching transcript."""
    candidates = [
        transcript
        for transcript in transcripts
        if transcript.intron_chain == read_model.intron_chain
        and _read_utrs_within_transcript(read_model, transcript)
    ]
    if not candidates:
        return None
    return sorted(
        candidates,
        key=lambda transcript: _assignment_sort_key(transcript, read_model),
    )[0]


def index_transcripts_by_intron_chain(
    transcripts: list[TranscriptAnnotation],
) -> dict[tuple[tuple[int, int], ...], list[TranscriptAnnotation]]:
    """Group transcript candidates by exact intron chain for fast assignment."""
    indexed: dict[tuple[tuple[int, int], ...], list[TranscriptAnnotation]] = {}
    for transcript in transcripts:
        indexed.setdefault(transcript.intron_chain, []).append(transcript)
    for intron_chain, candidates in indexed.items():
        indexed[intron_chain] = sorted(
            candidates,
            key=lambda transcript: (
                transcript.gene_name,
                transcript.gene_id,
                transcript.transcript_id,
                transcript.start,
            ),
        )
    return indexed


def assign_read_to_transcript_from_index(
    read_model: IsoSeqReadModel,
    transcript_index: dict[tuple[tuple[int, int], ...], list[TranscriptAnnotation]],
) -> TranscriptAnnotation | None:
    """Assign one read using a pre-indexed exact splice-chain candidate map."""
    candidates = [
        transcript
        for transcript in transcript_index.get(read_model.intron_chain, [])
        if _read_utrs_within_transcript(read_model, transcript)
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda transcript: _assignment_sort_key(transcript, read_model))


def assign_reads_to_transcripts(
    read_models: dict[str, IsoSeqReadModel], transcripts: list[TranscriptAnnotation]
) -> dict[str, IsoSeqAssignment]:
    """Assign read names to transcript IDs using exact splice chains and endpoint score."""
    transcript_index = index_transcripts_by_intron_chain(transcripts)
    return assign_reads_to_transcripts_from_index(read_models, transcript_index)


def assign_reads_to_transcripts_from_index(
    read_models: dict[str, IsoSeqReadModel],
    transcript_index: dict[tuple[tuple[int, int], ...], list[TranscriptAnnotation]],
) -> dict[str, IsoSeqAssignment]:
    """Assign read names using a pre-indexed transcript candidate map."""
    assignments = {}
    for read_name, read_model in read_models.items():
        transcript = assign_read_to_transcript_from_index(read_model, transcript_index)
        if transcript is None:
            assignments[read_name] = IsoSeqAssignment(None, None, None)
        else:
            assignments[read_name] = IsoSeqAssignment(
                transcript.transcript_id,
                transcript.gene_id,
                transcript.gene_name,
            )
    return assignments


def parse_transcripts_cached(
    gtf_file: str,
    region: Region,
    cache: dict[tuple[str, str, int, int], list[TranscriptAnnotation]],
) -> list[TranscriptAnnotation]:
    """Parse transcript annotations once per exact region during one process call."""
    cache_key = (gtf_file, region.chromosome, region.start, region.end)
    if cache_key not in cache:
        start_time = time.perf_counter()
        cache[cache_key] = _parse_transcripts_disk_cached(gtf_file, region)
        logger.debug(
            "IsoSeq transcript parsing for %s returned %d transcripts in %.3fs.",
            region.coordinate_str,
            len(cache[cache_key]),
            time.perf_counter() - start_time,
        )
    return cache[cache_key]


def _transcript_disk_cache_path(gtf_file: str, region: Region) -> str:
    stat = os.stat(gtf_file)
    cache_root = os.environ.get(
        "OROGRAPHER_ISOSEQ_CACHE_DIR",
        os.path.join(tempfile.gettempdir(), "orographer_isoseq_cache"),
    )
    key = "|".join(
        [
            ISOSEQ_TRANSCRIPT_CACHE_VERSION,
            os.path.abspath(gtf_file),
            str(stat.st_mtime_ns),
            str(stat.st_size),
            region.chromosome,
            str(region.start),
            str(region.end),
        ]
    )
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()
    return os.path.join(cache_root, f"{digest}.pkl")


def _parse_transcripts_disk_cached(gtf_file: str, region: Region) -> list[TranscriptAnnotation]:
    try:
        cache_path = _transcript_disk_cache_path(gtf_file, region)
    except FileNotFoundError:
        return parse_transcript_annotation_file(gtf_file, region)

    try:
        with open(cache_path, "rb") as handle:
            cached = pickle.load(handle)
        if isinstance(cached, list):
            return cached
    except (FileNotFoundError, OSError, pickle.PickleError, EOFError) as exc:
        logger.debug("IsoSeq transcript cache read failed for %s: %s", cache_path, exc)

    transcripts = parse_transcript_annotation_file(gtf_file, region)
    cache_dir = os.path.dirname(cache_path)
    try:
        os.makedirs(cache_dir, exist_ok=True)
        fd, tmp_path = tempfile.mkstemp(
            prefix=".tmp_transcripts_",
            suffix=".pkl",
            dir=cache_dir,
        )
        with os.fdopen(fd, "wb") as handle:
            pickle.dump(transcripts, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_path, cache_path)
    except OSError:
        with suppress(NameError, OSError):
            os.remove(tmp_path)
    return transcripts


def _display_rows(rows_in_bam_order: list[dict]) -> list[dict]:
    """Return rows in Orographer display order: other BAMs above primary."""
    if len(rows_in_bam_order) == 3:
        return [rows_in_bam_order[1], rows_in_bam_order[2], rows_in_bam_order[0]]
    if len(rows_in_bam_order) == 2:
        return [rows_in_bam_order[1], rows_in_bam_order[0]]
    return [rows_in_bam_order[0]]


def _gene_transcript_counts(transcripts: list[TranscriptAnnotation]) -> dict[str, int]:
    counts = {}
    for transcript in transcripts:
        counts[transcript.gene_id] = counts.get(transcript.gene_id, 0) + 1
    return counts


def build_isoseq_groups(
    segments_by_read: dict[str, list],
    assignments: dict[str, IsoSeqAssignment],
    transcripts: list[TranscriptAnnotation],
) -> list[dict]:
    """Build transcript-group metadata for plotting."""
    gene_counts = _gene_transcript_counts(transcripts)
    assigned_by_transcript = {transcript.transcript_id: [] for transcript in transcripts}
    unassigned = []
    display_read_start = {
        read_name: segments[0].pos for read_name, segments in segments_by_read.items() if segments
    }

    for display_read_name, segments in segments_by_read.items():
        raw_read_name = segments[0].readname if segments else display_read_name
        assignment = assignments.get(raw_read_name, IsoSeqAssignment(None, None, None))
        if assignment.transcript_id in assigned_by_transcript:
            assigned_by_transcript[assignment.transcript_id].append(display_read_name)
        else:
            unassigned.append(display_read_name)

    groups = []
    for transcript in transcripts:
        read_names = sorted(
            assigned_by_transcript.get(transcript.transcript_id, []),
            key=lambda name: display_read_start.get(name, 0),
        )
        groups.append(
            {
                "group_id": transcript.transcript_id,
                "transcript": transcript,
                "read_names": read_names,
                "assigned_read_count": len(read_names),
                "gene_transcript_count": gene_counts.get(transcript.gene_id, 1),
                "is_unassigned": False,
            }
        )

    groups.sort(
        key=lambda group: (
            group["transcript"].gene_name,
            group["transcript"].start,
            group["transcript"].transcript_id,
        )
    )
    groups.append(
        {
            "group_id": "UNASSIGNED",
            "transcript": None,
            "read_names": sorted(
                unassigned,
                key=lambda name: display_read_start.get(name, 0),
            ),
            "assigned_read_count": len(unassigned),
            "gene_transcript_count": 0,
            "is_unassigned": True,
        }
    )
    return groups


def process_isoseq(
    coordinate_strs: list[str],
    bam_files: list[str],
    vcf_files: list[str | None],
    reference_path: str | None,
    gtf_file: str,
    sample_labels: list[str | None],
    region_type: str,
    force_isoseq: bool = False,
) -> list[dict]:
    """Build IsoSeq region data in the shape expected by plot_reads_bokeh."""
    total_start_time = time.perf_counter()
    for bam_path in bam_files:
        validation_start_time = time.perf_counter()
        validate_isoseq_bam_header(bam_path, force=force_isoseq)
        logger.debug(
            "IsoSeq header validation for %s completed in %.3fs.",
            bam_path,
            time.perf_counter() - validation_start_time,
        )

    region_data_list = []
    transcript_cache: dict[tuple[str, str, int, int], list[TranscriptAnnotation]] = {}
    for coordinate_str in coordinate_strs:
        region_start_time = time.perf_counter()
        chromosome, start, end = parse_coordinate(coordinate_str)
        region = Region(chromosome, start, end, coordinate_str)
        transcripts = parse_transcripts_cached(gtf_file, region, transcript_cache)
        transcript_index = index_transcripts_by_intron_chain(transcripts)
        rows_in_bam_order = []

        for bam_index, bam_path in enumerate(bam_files):
            bam_start_time = time.perf_counter()
            vcf_path = vcf_files[bam_index] if bam_index < len(vcf_files) else None
            segments_by_read, read_models = collect_isoseq_segments_and_models(
                bam_path,
                region,
            )
            assignments = assign_reads_to_transcripts_from_index(
                read_models,
                transcript_index,
            )
            logger.debug(
                "IsoSeq BAM read/assignment pass for %s %s collected %d display reads, "
                "%d read models, and %d assignments in %.3fs.",
                bam_path,
                region.coordinate_str,
                len(segments_by_read),
                len(read_models),
                len(assignments),
                time.perf_counter() - bam_start_time,
            )
            isoseq_groups = build_isoseq_groups(segments_by_read, assignments, transcripts)
            rows_in_bam_order.append(
                {
                    "bam_path": bam_path,
                    "reference_path": reference_path,
                    "region_chromosome": region.chromosome,
                    "region_coordinate_str": region.coordinate_str,
                    "segments_by_read": segments_by_read,
                    "vcf_variants": parse_vcf_file(vcf_path, region),
                    "coverage_tracks": None,
                    "insertion_summary": {},
                    "region_type": region_type,
                    "sample_label": (
                        sample_labels[bam_index] if bam_index < len(sample_labels) else None
                    ),
                    "sample_index": bam_index,
                    "isoseq_groups": isoseq_groups,
                    "transcript_annotations": transcripts,
                }
            )

        region_data_list.append(
            {
                "region": region,
                "gene_annotations": [],
                "transcript_annotations": transcripts,
                "isoseq": True,
                "bam_rows": _display_rows(rows_in_bam_order),
            }
        )
        logger.debug(
            "IsoSeq region %s processing completed in %.3fs.",
            region.coordinate_str,
            time.perf_counter() - region_start_time,
        )

    logger.debug(
        "IsoSeq processing for %d region(s) and %d BAM(s) completed in %.3fs.",
        len(coordinate_strs),
        len(bam_files),
        time.perf_counter() - total_start_time,
    )
    return region_data_list
