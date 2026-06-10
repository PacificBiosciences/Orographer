from pathlib import Path

import pysam

from orographer.gtf_parser import GeneAnnotation, TranscriptAnnotation, parse_annotation_file
from orographer.utils import Region
from tests.helpers.synth_vcf_gtf import write_bgzip_tabix_gtf


def _write_indexed_gtf(path_gtf_gz: Path, lines: list[str]) -> None:
    tmp_plain = path_gtf_gz.with_suffix("")
    sorted_lines = sorted(lines, key=lambda line: (line.split("\t")[0], int(line.split("\t")[3])))
    tmp_plain.write_text("\n".join(sorted_lines) + "\n", encoding="utf-8")
    pysam.tabix_compress(str(tmp_plain), str(path_gtf_gz), force=True)
    pysam.tabix_index(str(path_gtf_gz), preset="gff", force=True)
    tmp_plain.unlink()


def test_annotations_overlap_region_at_inclusive_boundaries() -> None:
    gene = GeneAnnotation(
        gene_id="G1",
        gene_name="Gene1",
        chrom="chr1",
        start=100,
        end=200,
        strand="+",
    )
    transcript = TranscriptAnnotation(
        transcript_id="T1",
        transcript_name="Tx1",
        gene_id="G1",
        gene_name="Gene1",
        chrom="chr1",
        start=100,
        end=200,
        strand="+",
    )

    assert gene.overlaps_region(200, 250)
    assert gene.overlaps_region(50, 100)
    assert transcript.overlaps_region(200, 250)
    assert transcript.overlaps_region(50, 100)


def test_parse_annotation_file_returns_overlapping_gene(tmp_path: Path):
    gtf_gz = tmp_path / "test.gtf.gz"
    write_bgzip_tabix_gtf(gtf_gz)

    # Region overlapping gene/exon
    region = Region(chromosome="chr1", start=150, end=160, coordinate_str="chr1:150-160")
    annotations = parse_annotation_file(str(gtf_gz), region)

    assert len(annotations) >= 1
    gene_ids = {a.gene_id for a in annotations}
    assert "G1" in gene_ids

    g1 = next(a for a in annotations if a.gene_id == "G1")
    assert g1.gene_name == "Gene1"
    assert g1.chrom == "chr1"
    assert g1.start <= 150
    assert g1.end >= 160
    assert len(g1.exons) >= 1


def test_parse_annotation_file_prefers_representative_transcript_tag(tmp_path: Path):
    gtf_gz = tmp_path / "representative.gtf.gz"
    _write_indexed_gtf(
        gtf_gz,
        [
            'chr1\ttest\tgene\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "Gene1";',
            'chr1\ttest\ttranscript\t100\t300\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_primary"; transcript_name "Primary"; tag "GENCODE_Primary";',
            'chr1\ttest\texon\t100\t150\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_primary"; exon_number "1";',
            'chr1\ttest\texon\t250\t300\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_primary"; exon_number "2";',
            'chr1\ttest\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_long"; transcript_name "Long";',
            'chr1\ttest\texon\t100\t150\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_long"; exon_number "1";',
            'chr1\ttest\texon\t250\t300\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_long"; exon_number "2";',
            'chr1\ttest\texon\t400\t500\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_long"; exon_number "3";',
        ],
    )
    region = Region(chromosome="chr1", start=120, end=450, coordinate_str="chr1:120-450")

    annotations = parse_annotation_file(str(gtf_gz), region)

    assert len(annotations) == 1
    gene = annotations[0]
    assert gene.representative_transcript_id == "T_primary"
    assert gene.representative_selection_method == "GENCODE Primary"
    assert gene.exons == [(100, 150, 1), (250, 300, 2)]


def test_parse_annotation_file_falls_back_to_most_exons(tmp_path: Path):
    gtf_gz = tmp_path / "fallback.gtf.gz"
    _write_indexed_gtf(
        gtf_gz,
        [
            'chr1\ttest\ttranscript\t100\t300\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_short";',
            'chr1\ttest\texon\t100\t150\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_short";',
            'chr1\ttest\ttranscript\t100\t500\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_long";',
            'chr1\ttest\texon\t100\t150\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_long";',
            'chr1\ttest\texon\t250\t300\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T_long";',
        ],
    )
    region = Region(chromosome="chr1", start=120, end=450, coordinate_str="chr1:120-450")

    annotations = parse_annotation_file(str(gtf_gz), region)

    assert len(annotations) == 1
    gene = annotations[0]
    assert gene.representative_transcript_id == "T_long"
    assert gene.representative_selection_method == "most exons"
    assert gene.exons == [(100, 150, 1), (250, 300, 2)]


def test_parse_annotation_file_preserves_representative_transcript_span(tmp_path: Path):
    gtf_gz = tmp_path / "transcript_span.gtf.gz"
    _write_indexed_gtf(
        gtf_gz,
        [
            'chr1\ttest\ttranscript\t100\t400\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T1";',
            'chr1\ttest\texon\t200\t250\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'transcript_id "T1"; exon_number "2";',
        ],
    )
    region = Region(chromosome="chr1", start=190, end=260, coordinate_str="chr1:190-260")

    annotations = parse_annotation_file(str(gtf_gz), region)

    assert len(annotations) == 1
    assert annotations[0].start == 100
    assert annotations[0].end == 400
    assert annotations[0].exons == [(200, 250, 2)]


def test_missing_index_raises(tmp_path: Path):
    plain = tmp_path / "test.gtf.gz"
    # Create a placeholder so the parser hits validate_bgzip_index, but no .tbi exists.
    plain.write_text('chr1\t.\texon\t1\t2\t.\t+\t.\tgene_id "G1";\n', encoding="utf-8")

    region = Region(chromosome="chr1", start=1, end=2, coordinate_str="chr1:1-2")

    # validate_bgzip_index raises ValueError if `.tbi` is missing
    import pytest

    with pytest.raises(ValueError):
        parse_annotation_file(str(plain), region)
