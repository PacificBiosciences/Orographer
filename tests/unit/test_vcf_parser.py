from pathlib import Path

from orographer.utils import Region
from orographer.vcf_parser import parse_vcf_file
from tests.helpers.synth_vcf_gtf import write_simple_uncompressed_vcf


def test_parse_vcf_file_classifies_variants_and_haplotypes(tmp_path: Path):
    vcf_path = tmp_path / "test.vcf"
    write_simple_uncompressed_vcf(vcf_path)

    # Include positions 10..20
    region = Region(chromosome="chr1", start=9, end=20, coordinate_str="chr1:9-20")

    variants = parse_vcf_file(str(vcf_path), region)

    # We expect 3 variants inside the region (pos 10, 15, 18)
    assert len(variants) == 3

    types_by_pos = {v.pos: v.variant_type for v in variants}
    assert types_by_pos[10] == "SNP"
    assert types_by_pos[15] == "INSERTION"
    assert types_by_pos[18] == "DELETION"

    # Haplotypes: sampleA should appear for non-ref genotypes we set
    for v in variants:
        assert "sampleA" in v.haplotypes


def test_parse_vcf_file_outside_region_returns_empty(tmp_path: Path):
    vcf_path = tmp_path / "test.vcf"
    write_simple_uncompressed_vcf(vcf_path)

    region = Region(chromosome="chr1", start=50, end=60, coordinate_str="chr1:50-60")
    variants = parse_vcf_file(str(vcf_path), region)
    assert variants == []
