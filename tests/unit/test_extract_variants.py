from __future__ import annotations

from pathlib import Path

import pysam

from orographer.bam_parser import extract_variants


def _write_fasta(tmp_path: Path, ref_name: str, ref_seq: str) -> Path:
    fasta_path = tmp_path / f"{ref_name}.fa"
    fasta_path.write_text(f">{ref_name}\n{ref_seq}\n", encoding="utf-8")
    pysam.faidx(str(fasta_path))
    return fasta_path


def _write_single_record_bam(
    tmp_path: Path,
    ref_name: str,
    ref_seq: str,
    query_name: str,
    query_seq: str,
    reference_start: int,
    cigartuples: list[tuple[int, int]],
) -> Path:
    """
    Write a minimal BAM containing a single aligned read.

    `cigartuples` uses pysam CIGAR codes:
      - 0: M
      - 1: I
      - 2: D
    """

    bam_path = tmp_path / "test.bam"
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"LN": len(ref_seq), "SN": ref_name}],
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        seg = pysam.AlignedSegment()
        seg.query_name = query_name
        seg.query_sequence = query_seq
        seg.flag = 0
        seg.mapping_quality = 60
        seg.reference_id = 0
        seg.reference_start = reference_start
        seg.cigartuples = cigartuples
        bam.write(seg)

    return bam_path


def _load_first_record(bam_path: Path) -> pysam.AlignedSegment:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for rec in bam.fetch(until_eof=True):
            return rec
    raise AssertionError("No BAM records found")


def test_extract_variants_snps(tmp_path: Path):
    # Reference: positions 0..3 are A A A A
    # Read:      A A A T (SNP at ref_pos=3 -> alt_base=T)
    ref_name = "chr1"
    ref_seq = "AAAAACCCCC"

    fasta_path = _write_fasta(tmp_path, ref_name, ref_seq)
    bam_path = _write_single_record_bam(
        tmp_path=tmp_path,
        ref_name=ref_name,
        ref_seq=ref_seq,
        query_name="q1",
        query_seq="AAAT",
        reference_start=0,
        cigartuples=[(0, 4)],  # 4M
    )
    rec = _load_first_record(bam_path)
    ref_fa = pysam.FastaFile(str(fasta_path))

    mismatches, insertions, deletions = extract_variants(rec, ref_fa)

    assert insertions == []
    assert deletions == []
    assert mismatches == [(3, "T")]


def test_extract_variants_insertion(tmp_path: Path):
    # ref:  positions 0-1 are AA, positions 2-4 are AAA
    # read: AA + T + AAA
    ref_name = "chr1"
    ref_seq = "AAAAA"  # length 5

    fasta_path = _write_fasta(tmp_path, ref_name, ref_seq)
    bam_path = _write_single_record_bam(
        tmp_path=tmp_path,
        ref_name=ref_name,
        ref_seq=ref_seq,
        query_name="q1",
        query_seq="AATAAA",
        reference_start=0,
        cigartuples=[(0, 2), (1, 1), (0, 3)],  # 2M1I3M
    )
    rec = _load_first_record(bam_path)
    ref_fa = pysam.FastaFile(str(fasta_path))

    mismatches, insertions, deletions = extract_variants(rec, ref_fa)

    assert mismatches == []
    assert deletions == []
    # insertion_start_ref uses last_ref_pos from aligned_pairs; after 2M it's ref_pos=1
    assert insertions == [(1, "T")]


def test_extract_variants_deletion(tmp_path: Path):
    # ref:   AAAA AA
    # read:    AAAA    (deletion of 2 bases in the middle)
    # CIGAR: 2M2D2M deletes ref_pos=2 and 3 (0-based), so expected deletion tuple is [2,4)
    ref_name = "chr1"
    ref_seq = "AAAAAA"  # length 6

    fasta_path = _write_fasta(tmp_path, ref_name, ref_seq)
    bam_path = _write_single_record_bam(
        tmp_path=tmp_path,
        ref_name=ref_name,
        ref_seq=ref_seq,
        query_name="q1",
        query_seq="AAAA",  # 2M + 2M
        reference_start=0,
        cigartuples=[(0, 2), (2, 2), (0, 2)],  # 2M2D2M
    )
    rec = _load_first_record(bam_path)
    ref_fa = pysam.FastaFile(str(fasta_path))

    mismatches, insertions, deletions = extract_variants(rec, ref_fa)

    assert mismatches == []
    assert insertions == []
    assert deletions == [(2, 4)]
