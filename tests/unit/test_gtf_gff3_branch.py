"""GFF3 branch, detect_annotation_format, and GFF3 attribute parsing."""

from pathlib import Path

from orographer.gtf_parser import (
    detect_annotation_format,
    extract_gff3_gene_info,
    parse_annotation_file,
    parse_gff3_attributes,
)
from orographer.utils import Region
from tests.helpers.synth_vcf_gtf import write_bgzip_tabix_gff3


def test_detect_annotation_format_gtf():
    assert detect_annotation_format("/data/annotations.gtf.gz") == "gtf"
    assert detect_annotation_format("foo.GTF.gz") == "gtf"


def test_detect_annotation_format_gff3_and_gff():
    assert detect_annotation_format("genes.gff3.gz") == "gff3"
    assert detect_annotation_format("/path/to/ref.gff.gz") == "gff3"


def test_detect_annotation_format_defaults_to_gtf():
    assert detect_annotation_format("compressed.gz") == "gtf"
    assert detect_annotation_format("annotation.bgz") == "gtf"


def test_parse_gff3_attributes():
    attrs = parse_gff3_attributes("ID=G1;Name=Gene1;Parent=T1")
    assert attrs["ID"] == "G1"
    assert attrs["Name"] == "Gene1"
    assert attrs["Parent"] == "T1"


def test_extract_gff3_gene_info_gene_and_transcript_and_exon():
    t2g: dict[str, str] = {}
    gid, gname, skip = extract_gff3_gene_info({"ID": "G1", "Name": "N1"}, "gene", t2g)
    assert gid == "G1" and gname == "N1" and not skip

    _, _, skip = extract_gff3_gene_info({"ID": "T1", "Parent": "G1"}, "mrna", t2g)
    assert skip and t2g["T1"] == "G1"

    gid, gname, skip = extract_gff3_gene_info({"Parent": "T1"}, "exon", t2g)
    assert gid == "G1" and gname == "" and not skip


def test_parse_annotation_file_gff3_overlapping_gene(tmp_path: Path):
    path = tmp_path / "annot.gff3.gz"
    write_bgzip_tabix_gff3(path, gene_id="GFF3G1", gene_name_attr="GffThreeGene")

    region = Region(chromosome="chr1", start=150, end=160, coordinate_str="chr1:150-160")
    annotations = parse_annotation_file(str(path), region)

    assert len(annotations) >= 1
    g = next(a for a in annotations if a.gene_id == "GFF3G1")
    assert g.gene_name == "GffThreeGene"
    assert g.chrom == "chr1"
    assert len(g.exons) >= 1


def test_parse_annotation_file_gff3_gene_name_falls_back_to_id(tmp_path: Path):
    path = tmp_path / "nofname.gff3.gz"
    write_bgzip_tabix_gff3(path, gene_id="G2_ONLY_ID", gene_name_attr=None)

    region = Region(chromosome="chr1", start=150, end=160, coordinate_str="chr1:150-160")
    annotations = parse_annotation_file(str(path), region)

    g = next(a for a in annotations if a.gene_id == "G2_ONLY_ID")
    assert g.gene_name == "G2_ONLY_ID"
