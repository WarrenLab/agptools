import os

import pytest
from agp import open_agp
from agp.split import ParsingError, breakpoints_type, run


def test_parse_breakpoints(tmp_path):
    """Test parsing a good breakpoints file"""
    with open(tmp_path / "test_breaks.tsv", "w") as breakpoints_file:
        print(
            "scaffold_16\t4258995,21066364\nscaffold_17\t1234,5679",
            file=breakpoints_file,
        )

    breakpoints = breakpoints_type(tmp_path / "test_breaks.tsv")
    assert len(breakpoints["scaffold_16"]) == 2
    assert 1234 in breakpoints["scaffold_17"]


@pytest.mark.parametrize(
    "agp_string",
    [
        ("scaffold_16\t4258995,210r66364\nscaffold_17\t1234,5679"),
        ("scaffold_16\t4258995,21066364\nscaffold_17"),
        ("scaffold_16\t4258995,21066364\nscaffold_17    1234,5679"),
        ("scaffold_16\t4258995,21066364\nscaffold_16\t1234,5679"),
    ],
)
def test_parse_breakpoints_bad_format(agp_string, tmp_path):
    """Test parsing a few different types of misformatted breakpoints files"""
    with open(tmp_path / "test_breaks.tsv", "w") as breakpoints_file:
        print(agp_string, file=breakpoints_file)

    with pytest.raises(ParsingError):
        breakpoints_type(tmp_path / "test_breaks.tsv")


# TODO these three tests would be nice as they would allow pinpointing
# where a problem is, but they are not necessary from a code coverage
# perspective as test_split_run covers all of these functions
def test_unoffset_rows():
    pass


def test_convert_rows():
    pass


def test_split_scaffold():
    pass


def test_split_run(tmp_path):
    """Test a full run of the split module"""
    with open(tmp_path / "test_breaks.tsv", "w") as breakpoints_file:
        print(
            "scaffold_16\t4258995,21066364\nscaffold_17\t812345", file=breakpoints_file
        )

    agp_in_path = os.path.join(os.path.dirname(__file__), "data", "test.agp")
    with open(tmp_path / "test_out.agp", "w") as test_out_agp:
        run(
            breakpoints_type(tmp_path / "test_breaks.tsv"),
            test_out_agp,
            open_agp(agp_in_path),
        )

    with open(tmp_path / "test_out.agp", "r") as test_out_agp, open(
        os.path.join(os.path.dirname(__file__), "data", "test_split.agp")
    ) as correct_agp:
        for line1, line2 in zip(test_out_agp, correct_agp):
            assert line1 == line2
