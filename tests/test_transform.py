from os.path import dirname, join

from agp import open_agp
from agp.bed import open_bed
from agp.transform import run


def test_run_join(tmp_path):
    bed_in_path = join(dirname(__file__), "data", "test.bed")
    agp_in_path = join(dirname(__file__), "data", "test.agp")
    with open(tmp_path / "test_out.bed", "w") as outfile:
        run(open_bed(bed_in_path), open_agp(agp_in_path), outfile)

    with open(tmp_path / "test_out.bed", "r") as bed_out, open(
        join(dirname(__file__), "data", "transformed.bed"), "r"
    ) as correct_bed:
        for line1, line2 in zip(bed_out, correct_bed):
            assert line1 == line2
