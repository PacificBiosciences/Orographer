import pytest

from orographer.utils import parse_coordinate


def test_parse_coordinate_valid():
    chrom, start, end = parse_coordinate("chr1:1-2")
    assert chrom == "chr1"
    assert start == 1
    assert end == 2


@pytest.mark.parametrize(
    "coord",
    [
        "chr1-1-2",  # missing ':'
        "chr1:1_2",  # wrong delimiter
        "chr1:1:2",  # wrong delimiter count
        "chr1:1-1",  # start == end
        "chr1:2-1",  # start > end
        "chr1:abc-2",  # non-numeric
    ],
)
def test_parse_coordinate_invalid(coord):
    with pytest.raises(ValueError):
        parse_coordinate(coord)


def test_validate_output_dir_creates_dir(tmp_path):
    # validate_output_dir is used indirectly by CLI, but we can smoke it directly by importing
    from orographer.utils import validate_output_dir

    out = tmp_path / "new_out"
    assert not out.exists()
    validate_output_dir(str(out))
    assert out.exists()
    assert out.is_dir()


def test_validate_output_dir_rejects_non_dir(tmp_path):
    from orographer.utils import validate_output_dir

    p = tmp_path / "file"
    p.write_text("x")

    with pytest.raises(ValueError):
        validate_output_dir(str(p))
