from pathlib import Path

from orographer.gtf_parser import parse_annotation_file
from orographer.utils import Region
from tests.helpers.synth_vcf_gtf import write_bgzip_tabix_gtf


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


def test_missing_index_raises(tmp_path: Path):
    plain = tmp_path / "test.gtf.gz"
    # Create a placeholder so the parser hits validate_bgzip_index, but no .tbi exists.
    plain.write_text('chr1\t.\texon\t1\t2\t.\t+\t.\tgene_id "G1";\n', encoding="utf-8")

    region = Region(chromosome="chr1", start=1, end=2, coordinate_str="chr1:1-2")

    # validate_bgzip_index raises ValueError if `.tbi` is missing
    import pytest

    with pytest.raises(ValueError):
        parse_annotation_file(str(plain), region)
