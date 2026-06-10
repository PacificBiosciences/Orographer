"""Unit tests for compute_coverage in bam_parser.py."""

import os
import re
from pathlib import Path

import pysam

from orographer.bam_parser import compute_coverage
from orographer.utils import Region


def _make_simple_bam(tmp_path: Path, alignments: list[dict]) -> str:
    """
    Build a minimal BAM file with the given alignments.

    Each alignment dict may have:
        chrom, start (0-based), cigar, mapq, hp (int HP tag, optional)
    """
    bam_path = str(tmp_path / "test.bam")
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10_000_000}],
        }
    )
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for aln in alignments:
            rec = pysam.AlignedSegment(header)
            rec.query_name = aln.get("name", "read1")
            rec.query_sequence = "A" * _cigar_query_len(aln["cigar"])
            rec.flag = 0
            rec.reference_id = 0
            rec.reference_start = aln["start"]
            rec.mapping_quality = aln.get("mapq", 60)
            rec.cigarstring = aln["cigar"]
            rec.query_qualities = pysam.qualitystring_to_array("I" * _cigar_query_len(aln["cigar"]))
            if "hp" in aln:
                rec.set_tag("HP", aln["hp"])
            bam.write(rec)
    pysam.sort("-o", bam_path + ".sorted", bam_path)
    os.rename(bam_path + ".sorted", bam_path)
    pysam.index(bam_path)
    return bam_path


def _cigar_query_len(cigar_str: str) -> int:
    """Length of query sequence for a CIGAR string."""
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_str):
        if op in ("M", "I", "S", "=", "X"):
            total += int(length)
    return max(total, 1)


class TestComputeCoverageBasic:
    def test_empty_region_returns_total_key(self, tmp_path):
        bam_path = _make_simple_bam(tmp_path, [])
        region = Region("chr1", 1001, 2000, "chr1:1001-2000")
        result = compute_coverage(bam_path, region)
        assert -1 in result
        _, y = result[-1]
        assert all(d == 0 for d in y)

    def test_single_read_total_coverage(self, tmp_path):
        # Read spans 1000..1099 (0-based) = positions 1001-1100 (1-based)
        bam_path = _make_simple_bam(
            tmp_path,
            [{"name": "r1", "chrom": "chr1", "start": 1000, "cigar": "100M", "mapq": 60}],
        )
        region = Region("chr1", 1001, 1100, "chr1:1001-1100")
        result = compute_coverage(bam_path, region)
        assert -1 in result
        _, y = result[-1]
        assert max(y) >= 1

    def test_large_deletion_creates_dip(self, tmp_path):
        # 50M 50D 50M: the 50 bp deletion gap exceeds min_cigar_gap (20), so depth dips.
        bam_path = _make_simple_bam(
            tmp_path,
            [{"name": "r1", "chrom": "chr1", "start": 1000, "cigar": "50M50D50M", "mapq": 60}],
        )
        region = Region("chr1", 1051, 1100, "chr1:1051-1100")
        result = compute_coverage(bam_path, region)
        _, y = result[-1]
        assert max(y) == 0

    def test_small_deletion_filled(self, tmp_path):
        # 50M 10D 50M: the 10 bp deletion gap is below min_cigar_gap, so it's merged/filled.
        bam_path = _make_simple_bam(
            tmp_path,
            [{"name": "r1", "chrom": "chr1", "start": 1000, "cigar": "50M10D50M", "mapq": 60}],
        )
        # Deletion region: 0-based [1050, 1060), which is 1-based [1051, 1060].
        region = Region("chr1", 1051, 1060, "chr1:1051-1060")
        result = compute_coverage(bam_path, region)
        _, y = result[-1]
        assert max(y) >= 1

    def test_haplotype_keys_present(self, tmp_path):
        bam_path = _make_simple_bam(
            tmp_path,
            [
                {
                    "name": "r1",
                    "chrom": "chr1",
                    "start": 1000,
                    "cigar": "100M",
                    "mapq": 60,
                    "hp": 1,
                },
                {
                    "name": "r2",
                    "chrom": "chr1",
                    "start": 1000,
                    "cigar": "100M",
                    "mapq": 60,
                    "hp": 2,
                },
            ],
        )
        region = Region("chr1", 1001, 1100, "chr1:1001-1100")
        result = compute_coverage(bam_path, region, haplotypes_to_track={1, 2})
        assert 1 in result
        assert 2 in result
        _, y1 = result[1]
        _, y2 = result[2]
        assert max(y1) >= 1
        assert max(y2) >= 1

    def test_observed_string_haplotypes_are_tracked(self, tmp_path):
        bam_path = _make_simple_bam(
            tmp_path,
            [
                {
                    "name": "r1",
                    "chrom": "chr1",
                    "start": 1000,
                    "cigar": "100M",
                    "mapq": 60,
                    "hp": "hba_hba1hap1",
                },
                {
                    "name": "r2",
                    "chrom": "chr1",
                    "start": 1000,
                    "cigar": "100M",
                    "mapq": 60,
                    "hp": "hba_hba2hap2",
                },
            ],
        )
        region = Region("chr1", 1001, 1100, "chr1:1001-1100")

        result = compute_coverage(bam_path, region, include_observed_haplotypes=True)

        assert "hba_hba1hap1" in result
        assert "hba_hba2hap2" in result
        assert max(result["hba_hba1hap1"][1]) >= 1
        assert max(result["hba_hba2hap2"][1]) >= 1

    def test_haplotype_none_emits_only_total(self, tmp_path):
        bam_path = _make_simple_bam(
            tmp_path,
            [{"name": "r1", "chrom": "chr1", "start": 1000, "cigar": "100M", "mapq": 60, "hp": 1}],
        )
        region = Region("chr1", 1001, 1100, "chr1:1001-1100")
        result = compute_coverage(bam_path, region, haplotypes_to_track=None)
        assert -1 in result
        assert 1 not in result

    def test_x_positions_match_bin_size(self, tmp_path):
        bam_path = _make_simple_bam(tmp_path, [])
        region = Region("chr1", 1001, 1100, "chr1:1001-1100")
        result = compute_coverage(bam_path, region, bin_size=10)
        x, _ = result[-1]
        # x_vals = range(1001, 1101, 10) => 10 values
        assert len(x) == 10
        assert x[0] == 1001

    def test_two_reads_stack(self, tmp_path):
        bam_path = _make_simple_bam(
            tmp_path,
            [
                {"name": "r1", "chrom": "chr1", "start": 1000, "cigar": "100M", "mapq": 60},
                {"name": "r2", "chrom": "chr1", "start": 1000, "cigar": "100M", "mapq": 60},
            ],
        )
        region = Region("chr1", 1001, 1100, "chr1:1001-1100")
        result = compute_coverage(bam_path, region)
        _, y = result[-1]
        assert max(y) >= 2
