from pathlib import Path

import pysam
import pytest

from orographer.gtf_parser import parse_transcript_annotation_file
from orographer.isoseq import (
    IsoSeqReadModel,
    _read_model_from_record,
    _reference_blocks_split_on_refskip,
    assign_reads_to_transcripts,
    assign_reads_to_transcripts_from_index,
    index_transcripts_by_intron_chain,
    missing_isoseq_header_requirements,
    parse_transcripts_cached,
    process_isoseq,
    validate_isoseq_bam_header,
)
from orographer.plot_bokeh.plot_bokeh import (
    _isoseq_coverage_block_cache,
    _isoseq_coverage_for_cached_reads,
)
from orographer.utils import Region


def _write_bam_with_pg(path: Path, pg_records: list[dict]) -> None:
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "PG": pg_records,
    }
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    pysam.index(str(path))


def _write_transcript_gtf(path: Path) -> None:
    plain = path.with_suffix("")
    attrs_t1 = 'gene_id "G1"; gene_name "Gene1"; transcript_id "T1"; transcript_name "Tx1";'
    attrs_t2 = 'gene_id "G1"; gene_name "Gene1"; transcript_id "T2"; transcript_name "Tx2";'
    attrs_t3 = 'gene_id "G2"; gene_name "Gene2"; transcript_id "T3"; transcript_name "Tx3";'
    lines = [
        f"chr1\ttest\ttranscript\t100\t300\t.\t+\t.\t{attrs_t1}",
        f'chr1\ttest\texon\t100\t150\t.\t+\t.\t{attrs_t1} exon_number "1";',
        f"chr1\ttest\ttranscript\t110\t310\t.\t+\t.\t{attrs_t2}",
        f'chr1\ttest\texon\t110\t150\t.\t+\t.\t{attrs_t2} exon_number "1";',
        f'chr1\ttest\texon\t250\t300\t.\t+\t.\t{attrs_t1} exon_number "2";',
        f'chr1\ttest\texon\t250\t310\t.\t+\t.\t{attrs_t2} exon_number "2";',
        f"chr1\ttest\ttranscript\t500\t560\t.\t-\t.\t{attrs_t3}",
        f'chr1\ttest\texon\t500\t560\t.\t-\t.\t{attrs_t3} exon_number "1";',
    ]
    plain.write_text("\n".join(lines) + "\n", encoding="utf-8")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="gff", force=True)
    plain.unlink()


def _aligned_record(cigar: list[tuple[int, int]]) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment()
    record.query_name = "read_1"
    record.flag = 0
    record.reference_id = 0
    record.reference_start = 99
    record.mapping_quality = 60
    record.cigartuples = cigar
    query_length = sum(length for op, length in cigar if op in {0, 1, 4, 7, 8})
    record.query_sequence = "A" * query_length
    return record


def test_isoseq_header_validation_detects_required_pg_records(tmp_path: Path):
    bam = tmp_path / "ok.bam"
    _write_bam_with_pg(
        bam,
        [
            {"ID": "lima", "PN": "lima", "CL": "lima input.bam output.bam"},
            {"ID": "refine", "PN": "isoseq", "CL": "isoseq refine in.bam primers.fasta out.bam"},
            {
                "ID": "pbmm2",
                "PN": "pbmm2",
                "CL": "pbmm2 align ref.fa in.bam out.bam --preset ISOSEQ",
            },
        ],
    )

    assert missing_isoseq_header_requirements(str(bam)) == []
    validate_isoseq_bam_header(str(bam))


def test_isoseq_header_validation_reports_multiple_missing_items(tmp_path: Path):
    bam = tmp_path / "bad.bam"
    _write_bam_with_pg(
        bam,
        [{"ID": "pbmm2", "PN": "pbmm2", "CL": "pbmm2 align ref.fa in.bam out.bam"}],
    )

    assert missing_isoseq_header_requirements(str(bam)) == [
        "lima",
        "isoseq refine",
        "pbmm2 align --preset ISOSEQ",
    ]
    with pytest.raises(ValueError, match="missing lima, isoseq refine, pbmm2 align"):
        validate_isoseq_bam_header(str(bam))
    validate_isoseq_bam_header(str(bam), force=True)


def test_isoseq_reference_blocks_split_on_splice_skips_not_indels():
    blocks = _reference_blocks_split_on_refskip(
        [
            (0, 40),
            (1, 3),
            (0, 20),
            (2, 5),
            (0, 10),
            (3, 100),
            (0, 30),
        ],
        99,
    )

    assert blocks == [(99, 174), (274, 304)]


def test_isoseq_read_model_introns_ignore_deletions():
    record = _aligned_record(
        [
            (0, 40),
            (2, 5),
            (0, 30),
            (3, 100),
            (0, 30),
        ]
    )

    model = _read_model_from_record(record)

    assert model is not None
    assert model.start == 100
    assert model.end == 304
    assert model.intron_chain == ((174, 275),)


def test_parse_transcript_annotation_file_preserves_transcripts(tmp_path: Path):
    gtf = tmp_path / "tx.gtf.gz"
    _write_transcript_gtf(gtf)

    transcripts = parse_transcript_annotation_file(str(gtf), Region("chr1", 1, 800, "chr1:1-800"))

    by_id = {transcript.transcript_id: transcript for transcript in transcripts}
    assert set(by_id) == {"T1", "T2", "T3"}
    assert by_id["T1"].transcript_name == "Tx1"
    assert by_id["T1"].gene_name == "Gene1"
    assert by_id["T1"].strand == "+"
    assert by_id["T1"].exons == [(100, 150, 1), (250, 300, 2)]
    assert by_id["T1"].intron_chain == ((150, 250),)


def test_assign_reads_to_transcripts_exact_chain_and_endpoint_tie_break(tmp_path: Path):
    gtf = tmp_path / "tx.gtf.gz"
    _write_transcript_gtf(gtf)
    transcripts = parse_transcript_annotation_file(str(gtf), Region("chr1", 1, 800, "chr1:1-800"))
    read_models = {
        "best_t1": IsoSeqReadModel("best_t1", 100, 300, ((150, 250),)),
        "best_t2": IsoSeqReadModel("best_t2", 110, 310, ((150, 250),)),
        "shorter_t2": IsoSeqReadModel("shorter_t2", 112, 308, ((150, 250),)),
        "longer_five_prime": IsoSeqReadModel(
            "longer_five_prime",
            95,
            300,
            ((150, 250),),
        ),
        "longer_three_prime": IsoSeqReadModel(
            "longer_three_prime",
            110,
            315,
            ((150, 250),),
        ),
        "longer_both_ends": IsoSeqReadModel(
            "longer_both_ends",
            105,
            305,
            ((150, 250),),
        ),
        "single": IsoSeqReadModel("single", 500, 560, ()),
        "unassigned": IsoSeqReadModel("unassigned", 100, 300, ((151, 250),)),
    }

    assignments = assign_reads_to_transcripts(read_models, transcripts)

    assert assignments["best_t1"].transcript_id == "T1"
    assert assignments["best_t2"].transcript_id == "T2"
    assert assignments["shorter_t2"].transcript_id == "T2"
    assert assignments["longer_five_prime"].transcript_id is None
    assert assignments["longer_three_prime"].transcript_id is None
    assert assignments["longer_both_ends"].transcript_id is None
    assert assignments["single"].transcript_id == "T3"
    assert assignments["unassigned"].transcript_id is None


def test_indexed_assignment_matches_public_assignment(tmp_path: Path):
    gtf = tmp_path / "tx.gtf.gz"
    _write_transcript_gtf(gtf)
    transcripts = parse_transcript_annotation_file(str(gtf), Region("chr1", 1, 800, "chr1:1-800"))
    read_models = {
        "best_t1": IsoSeqReadModel("best_t1", 100, 300, ((150, 250),)),
        "best_t2": IsoSeqReadModel("best_t2", 110, 310, ((150, 250),)),
        "single": IsoSeqReadModel("single", 500, 560, ()),
        "unassigned": IsoSeqReadModel("unassigned", 100, 300, ((151, 250),)),
    }

    public_assignments = assign_reads_to_transcripts(read_models, transcripts)
    indexed_assignments = assign_reads_to_transcripts_from_index(
        read_models,
        index_transcripts_by_intron_chain(transcripts),
    )

    assert indexed_assignments == public_assignments


def test_parse_transcripts_cached_reuses_exact_region(monkeypatch):
    calls = []
    region = Region("chr1", 1, 800, "chr1:1-800")
    transcript = object()

    def fake_parse(gtf_file, parse_region):
        calls.append((gtf_file, parse_region))
        return [transcript]

    monkeypatch.setattr("orographer.isoseq.parse_transcript_annotation_file", fake_parse)
    cache = {}

    first = parse_transcripts_cached("annotation.gtf.gz", region, cache)
    second = parse_transcripts_cached("annotation.gtf.gz", region, cache)

    assert first == [transcript]
    assert second is first
    assert len(calls) == 1


class _Segment:
    def __init__(self, readname, pos, end):
        self.readname = readname
        self.pos = pos
        self.end = end


def test_isoseq_coverage_counts_all_assigned_reads():
    segments_by_read = {f"read_{i}": [_Segment(f"read_{i}", 100 + i, 240)] for i in range(125)}
    coverage_blocks = _isoseq_coverage_block_cache(
        segments_by_read,
        90,
        240,
    )
    coverage = _isoseq_coverage_for_cached_reads(
        list(segments_by_read),
        coverage_blocks,
        90,
        240,
    )

    assert max(coverage["y"]) == 125


def test_process_isoseq_uses_one_lightweight_bam_pass_per_region(monkeypatch):
    collect_calls = []
    transcript = type(
        "Transcript",
        (),
        {
            "transcript_id": "T1",
            "gene_id": "G1",
            "gene_name": "Gene1",
            "transcript_name": "Tx1",
            "chrom": "chr1",
            "start": 100,
            "end": 200,
            "strand": "+",
            "exons": [(100, 200, 1)],
            "intron_chain": (),
        },
    )()
    segment = _Segment("read_1", 100, 200)

    def fake_collect_segments_and_models(bam_path, region):
        collect_calls.append((bam_path, region))
        read_models = {"read_1": IsoSeqReadModel("read_1", 100, 200, ())}
        return {"read_1": [segment]}, read_models

    monkeypatch.setattr("orographer.isoseq.validate_isoseq_bam_header", lambda *_args, **_kw: None)
    monkeypatch.setattr(
        "orographer.isoseq.parse_transcript_annotation_file",
        lambda *_args: [transcript],
    )
    monkeypatch.setattr(
        "orographer.isoseq.collect_isoseq_segments_and_models",
        fake_collect_segments_and_models,
    )

    result = process_isoseq(
        ["chr1:1-300"],
        ["sample.bam"],
        [None],
        "ref.fa",
        "annotation.gtf.gz",
        [None],
        "isoseq",
    )

    assert len(collect_calls) == 1
    assert collect_calls[0][0] == "sample.bam"
    assert result[0]["bam_rows"][0]["isoseq_groups"][0]["assigned_read_count"] == 1
