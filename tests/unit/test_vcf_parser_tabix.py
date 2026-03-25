import json
from pathlib import Path

import pysam

from orographer.utils import Region, get_tabix_chromosome_name
from orographer.vcf_parser import parse_vcf_file


def _write_and_index_vcf_gz(tmp_path: Path) -> str:
    """
    Write a tiny VCF and create a bgzipped + tabix-indexed `.vcf.gz` using pysam.

    The VCF includes:
    - a SNP at chr1:10
    - an insertion at chr1:15
    - a deletion at chr1:18
    - one variant outside the test region at chr1:100
    """
    vcf_plain = tmp_path / "test.vcf"
    vcf_gz = tmp_path / "test.vcf.gz"

    header = "\n".join(
        [
            "##fileformat=VCFv4.2",
            '##source="orographer-tests"',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA",
        ]
    )

    lines = [
        # SNP: ref=A alt=T -> SNP
        "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1",
        # Insertion: ref=A alt=AT -> insertion (len(ref)<len(alt))
        "chr1\t15\t.\tA\tAT\t.\t.\t.\tGT\t1/1",
        # Deletion: ref=AT alt=A -> deletion (len(ref)>len(alt))
        "chr1\t18\t.\tAT\tA\t.\t.\t.\tGT\t0/1",
        # Outside region
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1",
    ]

    vcf_plain.write_text(header + "\n" + "\n".join(lines) + "\n", encoding="utf-8")

    # bgzip compress + tabix index
    pysam.tabix_compress(str(vcf_plain), str(vcf_gz), force=True)
    pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)

    return str(vcf_gz)


def test_parse_vcf_file_uses_tabix_branch(tmp_path: Path):
    vcf_gz = _write_and_index_vcf_gz(tmp_path)

    region = Region(chromosome="chr1", start=9, end=20, coordinate_str="chr1:9-20")
    variants = parse_vcf_file(vcf_gz, region)

    assert len(variants) == 3

    types_by_pos = {v.pos: v.variant_type for v in variants}
    assert types_by_pos[10] == "SNP"
    assert types_by_pos[15] == "INSERTION"
    assert types_by_pos[18] == "DELETION"

    for v in variants:
        assert v.haplotypes == ["sampleA"] or "sampleA" in v.haplotypes

    # Quick sanity: variants can be JSON-serialized (ensures no weird objects leak).
    json.dumps([v.__dict__ for v in variants])


def _write_indexed_vcf_numeric_chrom(tmp_path: Path) -> str:
    """VCF.gz with CHROM column `1` (no chr prefix), tabix-indexed."""
    vcf_plain = tmp_path / "numeric.vcf"
    vcf_gz = tmp_path / "numeric.vcf.gz"
    body = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##contig=<ID=1,length=249250621>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA",
            "1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1",
            "1\t15\t.\tA\tAT\t.\t.\t.\tGT\t1/1",
            "1\t100\t.\tG\tC\t.\t.\t.\tGT\t0/1",
        ]
    )
    vcf_plain.write_text(body + "\n", encoding="utf-8")
    pysam.tabix_compress(str(vcf_plain), str(vcf_gz), force=True)
    pysam.tabix_index(str(vcf_gz), preset="vcf", force=True)
    return str(vcf_gz)


def test_get_tabix_chromosome_name_maps_chr_prefix_to_numeric_contig(tmp_path: Path):
    vcf_gz = _write_indexed_vcf_numeric_chrom(tmp_path)
    with pysam.TabixFile(vcf_gz) as tabix_file:
        assert tabix_file.contigs == ["1"]
        assert get_tabix_chromosome_name(tabix_file, "chr1") == "1"
        assert get_tabix_chromosome_name(tabix_file, "1") == "1"


def test_parse_vcf_file_tabix_chr1_region_on_numeric_chrom_vcf(tmp_path: Path):
    """Region uses chr1; index contig is 1 — must resolve via get_tabix_chromosome_name."""
    vcf_gz = _write_indexed_vcf_numeric_chrom(tmp_path)
    region = Region(chromosome="chr1", start=9, end=20, coordinate_str="chr1:9-20")
    variants = parse_vcf_file(vcf_gz, region)
    assert {v.pos for v in variants} == {10, 15}


def test_parse_vcf_file_tabix_numeric_region_on_chr_prefixed_vcf(tmp_path: Path):
    """VCF rows use chr1; query with chromosome 1 — alias back to chr1."""
    vcf_gz = _write_and_index_vcf_gz(tmp_path)
    with pysam.TabixFile(vcf_gz) as tabix_file:
        assert get_tabix_chromosome_name(tabix_file, "1") == "chr1"
    region = Region(chromosome="1", start=9, end=20, coordinate_str="1:9-20")
    variants = parse_vcf_file(vcf_gz, region)
    assert {v.pos for v in variants} == {10, 15, 18}
