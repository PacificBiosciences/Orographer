from pathlib import Path
from types import SimpleNamespace

from orographer.plot_bokeh.data import generate_multi_region_filename, generate_output_filename
from orographer.utils import Region


def test_generate_output_filename_with_prefix(tmp_path: Path):
    seg = SimpleNamespace(chrom="chr1")
    segments_by_read = {"readA": [seg]}

    out = generate_output_filename(
        segments_by_read=segments_by_read,
        coordinate_start=100,
        coordinate_end=200,
        output_dir=str(tmp_path),
        prefix="sample_1",
    )

    assert out == str(tmp_path / "sample_1_chr1_100_200_bokeh.html")
    assert tmp_path.is_dir()


def test_generate_output_filename_without_prefix(tmp_path: Path):
    seg = SimpleNamespace(chrom="chrX")
    segments_by_read = {"readA": [seg]}

    out = generate_output_filename(
        segments_by_read=segments_by_read,
        coordinate_start=10,
        coordinate_end=20,
        output_dir=str(tmp_path),
        prefix=None,
    )

    assert out == str(tmp_path / "chrX_10_20_bokeh.html")


def test_generate_multi_region_filename_with_prefix(tmp_path: Path):
    # coordinate_str sanitization: ":" and "-" become "_"
    region1 = Region("chr1", 100, 200, "chr1:100-200")
    region2 = Region("chrX", 10, 20, "chrX:10-20")
    region_data_list = [{"region": region1}, {"region": region2}]

    out = generate_multi_region_filename(
        region_data_list=region_data_list,
        output_dir=str(tmp_path),
        prefix="pref",
    )

    # sanitize:
    #   chr1:100-200 -> chr1_100_200
    #   chrX:10-20    -> chrX_10_20
    expected_file = "pref_chr1_100_200_chrX_10_20_bokeh.html"
    assert out == str(tmp_path / expected_file)
    assert tmp_path.is_dir()


def test_generate_multi_region_filename_without_prefix(tmp_path: Path):
    region1 = Region("chr2", 1, 2, "2:1-2")
    region2 = Region("chr3", 5, 6, "chr3:5-6")
    region_data_list = [{"region": region1}, {"region": region2}]

    out = generate_multi_region_filename(
        region_data_list=region_data_list,
        output_dir=str(tmp_path),
        prefix=None,
    )

    expected_file = "2_1_2_chr3_5_6_bokeh.html"
    assert out == str(tmp_path / expected_file)
