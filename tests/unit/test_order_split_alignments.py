"""Unit tests for order_split_alignments: SA ordering, malformed SA, fallback tie-breaks."""

from unittest.mock import MagicMock

from orographer.bam_parser import order_split_alignments


def _aln(
    ref_name: str,
    ref_start: int,
    *,
    strand: str = "+",
    q_start: int = 0,
    q_end: int = 50,
    primary: bool = False,
    supplementary: bool = False,
    secondary: bool = False,
    sa_tag: str | None = None,
):
    """Minimal mock pysam.AlignedSegment for ordering tests."""
    m = MagicMock()
    m.reference_name = ref_name
    m.reference_start = ref_start
    m.is_reverse = strand == "-"
    m.is_secondary = secondary
    m.is_supplementary = supplementary
    m.query_alignment_start = q_start
    m.query_alignment_end = q_end

    if sa_tag is not None:

        def has_tag(name):
            return name == "SA"

        def get_tag(name):
            if name == "SA":
                return sa_tag
            raise KeyError(name)

        m.has_tag = has_tag
        m.get_tag = get_tag
    else:
        m.has_tag = lambda _n: False

    if primary:
        m.is_secondary = False
        m.is_supplementary = False
    return m


def test_order_split_alignments_sa_order_matches_tag():
    """Supplementary order follows SA list; primary (not in SA) sorts last."""
    sa = "chr1,5000,+,50M,60,0;chr1,6000,+,50M,60,0"
    primary = _aln("chr1", 1000, primary=True, sa_tag=sa)
    supp_a = _aln("chr1", 5000, supplementary=True, q_start=0, q_end=40)
    supp_b = _aln("chr1", 6000, supplementary=True, q_start=40, q_end=80)
    ordered = order_split_alignments([primary, supp_b, supp_a])
    assert [a.reference_start for a in ordered] == [5000, 6000, 1000]


def test_order_split_alignments_skips_malformed_sa_entries():
    """Malformed SA fragments are ignored; valid entries still define order."""
    sa = "chr1,notint,+,x;chr1,7000,+,50M,60,0"
    primary = _aln("chr1", 2000, primary=True, sa_tag=sa)
    supp = _aln("chr1", 7000, supplementary=True)
    ordered = order_split_alignments([primary, supp])
    assert ordered[0].reference_start == 7000
    assert ordered[-1].reference_start == 2000


def test_order_split_alignments_all_sa_malformed_falls_back_to_query_span():
    """No valid SA keys → interval fallback by query coordinates."""
    sa = "short;also,too,few,fields"
    primary = _aln("chr1", 1000, primary=True, q_start=0, q_end=30, sa_tag=sa)
    supp = _aln("chr1", 5000, supplementary=True, q_start=30, q_end=80)
    ordered = order_split_alignments([supp, primary])
    assert [a.query_alignment_start for a in ordered] == [0, 30]


def test_order_split_alignments_fallback_tie_longer_span_first():
    """Same query start: longer alignment sorts before shorter (sort key tie-break)."""
    a = _aln("chr1", 100, q_start=0, q_end=20, primary=True)
    b = _aln("chr1", 200, supplementary=True, q_start=0, q_end=50)
    ordered = order_split_alignments([a, b])
    # Same start 0; longer span (50) before shorter (20)
    assert ordered[0].query_alignment_end == 50
    assert ordered[1].query_alignment_end == 20


def test_order_split_alignments_fallback_primary_before_supplementary_on_tie():
    """Equal query span: primary before supplementary."""
    a = _aln("chr1", 100, q_start=0, q_end=40, primary=True)
    b = _aln("chr1", 200, supplementary=True, q_start=0, q_end=40)
    ordered = order_split_alignments([b, a])
    assert not ordered[0].is_supplementary
    assert ordered[1].is_supplementary


def test_order_split_alignments_empty():
    assert order_split_alignments([]) == []


def test_order_split_alignments_single():
    a = _aln("chr1", 100, primary=True)
    assert order_split_alignments([a]) == [a]
