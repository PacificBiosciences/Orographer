from pathlib import Path

import pysam
import pytest

from orographer.bam_parser import fetch_all_alignments, validate_bam_file
from orographer.utils import Region

# HG002.paraphase.small @ chr16:171801-175500 (hba): stable lower bounds for regression.
_MIN_READS_CHR16_PARAPHASE_NON_SPLIT = 60
_MIN_READS_CHR16_PARAPHASE_ONLY_SPLIT = 55


def _bed_first_coord_to_cli_region(bed_path: Path, chrom: str):
    """
    Convert BED (0-based start, end) into orographer CLI coordinate (1-based inclusive).
    The e2e script uses: coord="${chrom}:$((start + 1))-${end}"
    """
    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] != chrom:
                continue
            bed_chrom, start0, end, _name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            cli_start1 = start0 + 1
            coordinate_str = f"{bed_chrom}:{cli_start1}-{end}"
            return Region(bed_chrom, cli_start1, end, coordinate_str)
    raise AssertionError(f"No BED line for chrom={chrom} in {bed_path}")


def _infer_haplotype_from_key(key: str) -> int:
    # plot_bokeh/data.py uses 0/unassigned semantics; bam_parser keys look like:
    #   rname                 -> haplotype_tag is None/0
    #   rname_HP{hp}          -> haplotype_tag is hp (int)
    if "_HP" not in key:
        return 0
    # Split on last occurrence to be robust to rname containing "_HP"
    _prefix, hp_str = key.rsplit("_HP", 1)
    try:
        return int(hp_str)
    except ValueError:
        # If the key format is unexpected, defaulting to 0 makes the assertion clearer elsewhere.
        return 0


def _infer_haplotype_from_segment(seg) -> int:
    if seg.haplotype_tag is None:
        return 0
    try:
        return int(seg.haplotype_tag)
    except (TypeError, ValueError):
        return 0


def _assert_segment_invariants(seg) -> None:
    assert seg.spans is True
    assert seg.readname  # query_name
    assert seg.chrom
    assert seg.second_chrom
    assert seg.pos <= seg.end
    assert seg.fwd_read_start <= seg.fwd_read_end
    assert seg.alignment_order >= 1

    # These are populated during conversion; they should never be None.
    assert isinstance(seg.mismatches, list)
    assert isinstance(seg.insertions, list)
    assert isinstance(seg.deletions, list)


def _assert_segment_overlaps_query_region(seg, region: Region) -> None:
    """
    Region is 1-based inclusive; BAM segment uses 0-based ref start and exclusive end.
    Only enforced when the segment is on the queried chromosome (skip SA hits elsewhere).
    """
    if seg.chrom != region.chromosome:
        return
    # 1-based inclusive [start, end] -> 0-based half-open [start-1, end]
    r0 = region.start - 1
    r1_excl = region.end
    assert seg.pos < r1_excl and seg.end > r0, (
        f"Segment {seg.chrom}:{seg.pos}-{seg.end} should overlap "
        f"{region.chromosome}:{region.start}-{region.end}"
    )


def _assert_segments_dict_invariants(segments_by_hap: dict, region: Region | None = None) -> None:
    assert isinstance(segments_by_hap, dict)
    assert segments_by_hap  # expect at least one read/track for these fixtures

    for key, seg_list in segments_by_hap.items():
        assert key  # read key
        assert seg_list, f"Empty segment list for key={key}"

        # alignment_order should be sorted and stable within each list
        assert seg_list == sorted(seg_list, key=lambda s: s.alignment_order)
        assert all(seg.alignment_order >= 1 for seg in seg_list)

        expected_hp = _infer_haplotype_from_key(key)
        for seg in seg_list:
            _assert_segment_invariants(seg)
            if region is not None:
                _assert_segment_overlaps_query_region(seg, region)
            inferred_hp = _infer_haplotype_from_segment(seg)
            assert inferred_hp == expected_hp, (
                f"Haplotype key mismatch for key={key}: inferred_hp={inferred_hp}"
            )


def test_validate_bam_file(tmp_path: Path):
    missing = tmp_path / "missing.bam"
    with pytest.raises(FileNotFoundError):
        validate_bam_file(str(missing))


def test_fetch_all_alignments_non_split_smoke():
    bam_path = Path("tests/data/inputs/paraviewer/HG002.paraphase.small.bam")
    ref_path = Path("tests/data/inputs/hg38_small.fa.gz")
    bed_path = Path("tests/data/inputs/paraviewer/HG002_regions.bed")

    region = _bed_first_coord_to_cli_region(bed_path, chrom="chr16")

    segments = fetch_all_alignments(
        bam_path=str(bam_path),
        region=region,
        only_split=False,
        reference_path=str(ref_path),
    )

    assert isinstance(segments, dict)
    assert len(segments) >= _MIN_READS_CHR16_PARAPHASE_NON_SPLIT, (
        "Expected a stable minimum number of reads in chr16 hba for this fixture"
    )
    assert sum(len(v) for v in segments.values()) >= _MIN_READS_CHR16_PARAPHASE_NON_SPLIT
    _assert_segments_dict_invariants(segments, region=region)


def test_fetch_all_alignments_only_split_smoke():
    bam_path = Path("tests/data/inputs/paraviewer/HG002.paraphase.small.bam")
    ref_path = Path("tests/data/inputs/hg38_small.fa.gz")
    bed_path = Path("tests/data/inputs/paraviewer/HG002_regions.bed")

    # Same region as non-split: hba on chr16 has many SA-tagged reads in this BAM.
    region = _bed_first_coord_to_cli_region(bed_path, chrom="chr16")

    split_segments = fetch_all_alignments(
        bam_path=str(bam_path),
        region=region,
        only_split=True,
        reference_path=str(ref_path),
    )

    assert isinstance(split_segments, dict)
    assert len(split_segments) >= _MIN_READS_CHR16_PARAPHASE_ONLY_SPLIT, (
        "Fixture HG002.paraphase.small must yield SA-tagged reads in chr16 hba"
    )
    assert sum(len(v) for v in split_segments.values()) >= _MIN_READS_CHR16_PARAPHASE_ONLY_SPLIT
    _assert_segments_dict_invariants(split_segments, region=region)

    for seg_list in split_segments.values():
        assert seg_list == sorted(seg_list, key=lambda s: s.alignment_order)


def _write_fasta(tmp_path: Path, ref_name: str, ref_seq: str) -> str:
    fasta_path = tmp_path / f"{ref_name}.fa"
    fasta_path.write_text(f">{ref_name}\n{ref_seq}\n", encoding="utf-8")
    pysam.faidx(str(fasta_path))
    return str(fasta_path)


def _write_sa_tagged_bam(
    tmp_path: Path,
    ref_name: str,
    ref_seq: str,
    query_name: str,
    query_seq: str,
    reference_start0: int,
    primary_cigartuples: list[tuple[int, int]],
    sa_rname: str,
    sa_pos_1based: int,
    sa_strand: str,
    sa_cigar: str,
    sa_mapq: int,
    sa_nm: int,
) -> str:
    """
    Create a BAM with a single primary record containing an `SA` aux tag.

    The SA tag format used by `bam_parser.parse_sa_segment` is:
      rname,pos,strand,cigar,mapq,nm;
    where `pos` is 1-based in the tag.
    """
    bam_path = tmp_path / "sa_tagged.bam"
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
        seg.reference_start = reference_start0
        seg.cigartuples = primary_cigartuples

        # Attach SA tag to the primary record.
        sa_val = f"{sa_rname},{sa_pos_1based},{sa_strand},{sa_cigar},{sa_mapq},{sa_nm};"
        seg.set_tag("SA", sa_val)

        bam.write(seg)

    # `fetch_all_alignments()` uses `AlignmentFile.fetch()`, which requires an index.
    pysam.index(str(bam_path))

    return str(bam_path)


def test_fetch_all_alignments_only_split_sa_present(tmp_path: Path):
    """
    Deterministically validate `only_split=True` for a BAM where the primary record
    has an `SA` aux tag: it should yield at least the primary segment + SA-derived
    segment(s).
    """
    ref_name = "chr1"
    ref_seq = "AAAAACCCCC"  # length 10
    fasta_path = _write_fasta(tmp_path, ref_name, ref_seq)

    query_name = "read1"
    query_seq = "AAAA"  # length 4
    bam_path = _write_sa_tagged_bam(
        tmp_path=tmp_path,
        ref_name=ref_name,
        ref_seq=ref_seq,
        query_name=query_name,
        query_seq=query_seq,
        reference_start0=0,  # 0-based, maps to chr1:1
        primary_cigartuples=[(0, 4)],  # 4M
        sa_rname=ref_name,
        sa_pos_1based=1,
        sa_strand="+",
        sa_cigar="4M",
        sa_mapq=60,
        sa_nm=0,
    )

    region = Region(ref_name, start=1, end=10, coordinate_str="chr1:1-10")

    segments = fetch_all_alignments(
        bam_path=bam_path,
        region=region,
        only_split=True,
        reference_path=fasta_path,
    )

    assert isinstance(segments, dict)
    assert len(segments) == 1

    seg_list = segments[query_name]
    assert len(seg_list) >= 2, "Expected primary + at least one SA-derived segment"

    # Primary segment must be present and ordered first.
    assert seg_list[0].from_primary_bam_record is True
    assert seg_list[0].readname == query_name

    # All segments should have alignment_order assigned and valid.
    assert seg_list == sorted(seg_list, key=lambda s: s.alignment_order)
    assert all(s.alignment_order >= 1 for s in seg_list)
    for seg in seg_list:
        _assert_segment_overlaps_query_region(seg, region)
