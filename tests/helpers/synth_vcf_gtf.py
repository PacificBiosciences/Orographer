import contextlib
from pathlib import Path

import pysam


def write_simple_uncompressed_vcf(path: Path) -> None:
    """
    Write a minimal uncompressed VCF with:
    - one SNP in range
    - one insertion in range
    - one deletion in range
    - one variant outside range
    """
    header = "\n".join(
        [
            "##fileformat=VCFv4.2",
            '##source="orographer-tests"',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA",
        ]
    )

    # Coordinates below assume region start/end are chosen to include pos 10..20
    lines = [
        # SNP: ref=A alt=T
        "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1",
        # Insertion: ref=A alt=AT (len(ref) < len(alt))
        "chr1\t15\t.\tA\tAT\t.\t.\t.\tGT\t1/1",
        # Deletion: ref=AT alt=A (len(ref) > len(alt))
        "chr1\t18\t.\tAT\tA\t.\t.\t.\tGT\t0/1",
        # Outside region
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1",
    ]

    path.write_text(header + "\n" + "\n".join(lines) + "\n", encoding="utf-8")


def write_bgzip_tabix_gtf(path_gtf_gz: Path) -> None:
    """
    Write a tiny GTF file (bgzipped + tabix-indexed).

    `path_gtf_gz` should end with `.gtf.gz`.
    """
    if path_gtf_gz.suffixes[-2:] != [".gtf", ".gz"] and not str(path_gtf_gz).endswith(".gtf.gz"):
        raise ValueError("Expected a .gtf.gz path for indexed GTF creation.")

    tmp_plain = path_gtf_gz.with_suffix("")  # strip .gz

    # One gene (chr1:100-200) with one exon.
    # Attributes follow the parser's expectations: key "value";
    gtf_text = "\n".join(
        [
            'chr1\torographer-tests\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "Gene1";',
            'chr1\torographer-tests\texon\t120\t180\t.\t+\t.\tgene_id "G1"; gene_name "Gene1"; '
            'exon_number "1";',
        ]
    )
    tmp_plain.write_text(gtf_text + "\n", encoding="utf-8")

    # Create bgzip-compressed file + tabix index.
    # pysam.tabix_compress expects a plain input file and writes to output.
    pysam.tabix_compress(str(tmp_plain), str(path_gtf_gz), force=True)
    pysam.tabix_index(str(path_gtf_gz), preset="gff", force=True)

    # Cleanup plain file.
    with contextlib.suppress(FileNotFoundError):
        tmp_plain.unlink()


def write_bgzip_tabix_gff3(path_gff3_gz: Path, *, gene_id: str, gene_name_attr: str | None) -> None:
    """
    Write a minimal GFF3 (gene + mRNA + exon), bgzip + tabix-indexed.

    If gene_name_attr is None, the gene line omits Name= so the parser falls back to ID.
    """
    if not str(path_gff3_gz).endswith(".gff3.gz"):
        raise ValueError("Expected a .gff3.gz path for indexed GFF3 creation.")

    tmp_plain = path_gff3_gz.with_suffix("")  # .gff3
    tid = f"T_{gene_id}"
    gene_attrs = f"ID={gene_id};"
    if gene_name_attr is not None:
        gene_attrs += f"Name={gene_name_attr};"
    gff3_text = "\n".join(
        [
            f"chr1\torographer\tgene\t100\t200\t.\t+\t.\t{gene_attrs}",
            f"chr1\torographer\tmRNA\t100\t200\t.\t+\t.\tID={tid};Parent={gene_id};",
            f"chr1\torographer\texon\t120\t180\t.\t+\t.\tParent={tid};exon_number=1;",
        ]
    )
    tmp_plain.write_text(gff3_text + "\n", encoding="utf-8")
    pysam.tabix_compress(str(tmp_plain), str(path_gff3_gz), force=True)
    pysam.tabix_index(str(path_gff3_gz), preset="gff", force=True)
    with contextlib.suppress(FileNotFoundError):
        tmp_plain.unlink()
