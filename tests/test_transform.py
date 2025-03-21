from os.path import dirname, join
from unittest.mock import patch

import pytest

from agp import AgpRow, open_agp
from agp.agptools import main
from agp.bed import BedRange
from agp.transform import (
    BadOrientationError,
    CoordinateNotFoundError,
    NoSuchContigError,
    UnsupportedOperationError,
    create_contig_dict,
    find_agp_row,
    transform_bed_row,
    transform_single_position,
)


@pytest.fixture
def contig_dict():
    return create_contig_dict(open_agp(join(dirname(__file__), "data", "test.agp")))


def test_no_such_contig_error(contig_dict):
    with pytest.raises(NoSuchContigError) as err:
        find_agp_row("cotnig1", 100, contig_dict)
    assert 'contig called "cotnig1" in AGP' in str(err.value)


def test_bad_orientation_error():
    with pytest.raises(BadOrientationError) as err:
        transform_single_position(
            100, AgpRow("scaffold_16\t1\t167\t1\tW\tptg65l\t1\t167\t?")
        )
    assert "+/-:" in str(err.value)


def test_coordinate_not_found_error(contig_dict):
    with pytest.raises(CoordinateNotFoundError) as err:
        find_agp_row("tig00001012|arrow|arrow", 876238, contig_dict)
    assert "tig00001012|arrow|arrow:876238" in str(err.value)


def test_unsupported_operation_error_coordinateless(contig_dict):
    with pytest.raises(UnsupportedOperationError) as err:
        transform_bed_row(BedRange("tig00001012|arrow|arrow"), contig_dict)
    assert "coordinateless" in str(err.value)


def test_unsupported_operation_multiple_scaffolds():
    split_contig_dict = create_contig_dict(
        open_agp(join(dirname(__file__), "data", "test_split1.agp"))
    )
    with pytest.raises(UnsupportedOperationError) as err:
        transform_bed_row(
            BedRange("tig00005119|arrow|arrow", 500000, 600000), split_contig_dict
        )
    assert "spanning multiple" in str(err.value)


@pytest.mark.parametrize(
    "in_position, agp_row, out_position",
    [
        (1, AgpRow("scaffold_16\t1\t167\t1\tW\tptg65l\t1\t167\t+"), 1),
        (1, AgpRow("scaffold_16\t1\t167\t1\tW\tptg65l\t1\t167\t-"), 167),
        (150, AgpRow("scaffold_16\t1\t167\t1\tW\tptg65l\t1\t167\t+"), 150),
        (150, AgpRow("scaffold_16\t1\t167\t1\tW\tptg65l\t1\t167\t-"), 18),
        (150, AgpRow("scaffold_16\t1\t167\t1\tW\tptg65l\t1\t167\t-"), 18),
        (1, AgpRow("scaffold_16\t168\t200\t1\tW\tptg66l\t1\t33\t+"), 168),
        (1, AgpRow("scaffold_16\t168\t200\t1\tW\tptg66l\t1\t33\t-"), 200),
        (33, AgpRow("scaffold_16\t168\t200\t1\tW\tptg66l\t1\t33\t-"), 168),
        (33, AgpRow("scaffold_16\t168\t200\t1\tW\tptg66l\t1\t33\t+"), 200),
        (44, AgpRow("scaffold_16\t131\t150\t1\tW\tptg66l\t31\t50\t+"), 144),
        (44, AgpRow("scaffold_16\t131\t150\t1\tW\tptg66l\t31\t50\t-"), 137),
    ],
)
def test_transform_single_position(in_position, agp_row, out_position):
    assert transform_single_position(in_position, agp_row) == out_position


def test_transform_main(tmp_path):
    with patch(
        "sys.argv",
        [
            "agptools",
            "transform",
            join(dirname(__file__), "data", "test_transform.bed"),
            join(dirname(__file__), "data", "test_transform.agp"),
            "-o",
            join(tmp_path, "test_out.bed"),
        ],
    ):
        main()

    with (
        open(tmp_path / "test_out.bed", "r") as bed_out,
        open(join(dirname(__file__), "data", "transformed.bed"), "r") as correct_bed,
    ):
        for line1, line2 in zip(bed_out, correct_bed):
            assert line1 == line2
