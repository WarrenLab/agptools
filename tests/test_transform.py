from os.path import dirname, join

import pytest

from agp import AgpRow, open_agp
from agp.bed import open_bed
from agp.transform import run, transform_single_position


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


def test_run_transform(tmp_path):
    bed_in_path = join(dirname(__file__), "data", "test_transform.bed")
    agp_in_path = join(dirname(__file__), "data", "test_transform.agp")
    with open(tmp_path / "test_out.bed", "w") as outfile:
        run(open_bed(bed_in_path), open_agp(agp_in_path), outfile)

    with open(tmp_path / "test_out.bed", "r") as bed_out, open(
        join(dirname(__file__), "data", "transformed.bed"), "r"
    ) as correct_bed:
        for line1, line2 in zip(bed_out, correct_bed):
            assert line1 == line2
