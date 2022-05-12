import os

import pytest

from agp import AgpRow, open_agp
from agp.split import ParsingError, breakpoints_type, run, split_contig, split_scaffold


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


# TODO these two tests would be nice as they would allow pinpointing
# where a problem is, but they are not necessary from a code coverage
# perspective as test_split_run covers all of these functions
def test_unoffset_rows():
    pass


def test_convert_rows():
    pass


@pytest.mark.parametrize(
    "in_agp_strings,break_location,out_agp_strings",
    [
        (
            ["scfd1\t1\t1000\t1\tW\tctg1\t1\t1000\t+"],
            300,
            [
                "scfd1.1\t1\t300\t1\tW\tctg1\t1\t300\t+",
                "scfd1.2\t1\t700\t1\tW\tctg1\t301\t1000\t+",
            ],
        ),
        (
            ["scfd2\t1\t1000\t1\tW\tctg2\t1\t1000\t-"],
            300,
            [
                "scfd2.1\t1\t300\t1\tW\tctg2\t701\t1000\t-",
                "scfd2.2\t1\t700\t1\tW\tctg2\t1\t700\t-",
            ],
        ),
        (
            [
                "scfd3\t1\t1000\t1\tW\tctg3\t1\t1000\t+",
                "scfd3\t1001\t2000\t2\tW\tctg4\t1\t1000\t-",
                "scfd3\t2001\t4000\t3\tW\tctg5\t1\t2000\t+",
            ],
            1300,
            [
                "scfd3.1\t1\t1000\t1\tW\tctg3\t1\t1000\t+",
                "scfd3.1\t1001\t1300\t2\tW\tctg4\t701\t1000\t-",
                "scfd3.2\t1\t700\t1\tW\tctg4\t1\t700\t-",
                "scfd3.2\t701\t2700\t2\tW\tctg5\t1\t2000\t+",
            ],
        ),
    ],
)
def test_split_scaffold(in_agp_strings, break_location, out_agp_strings):
    """Test splitting a scaffold within a contig

    This pair of tests was provided by @songtaogui in issue #2. The
    first works as expected in commit 4813a52, but not the second.
    """
    agp_rows = [AgpRow(s) for s in in_agp_strings]
    broken_rows = split_scaffold(agp_rows, [break_location])

    for split_row, expected_split_row_string in zip(broken_rows, out_agp_strings):
        assert split_row == AgpRow(expected_split_row_string)


@pytest.mark.parametrize(
    "in_agp_string,out_agp_strings",
    [
        (
            "scfd1\t1\t1000\t1\tW\tctg1\t1\t1000\t+",
            [
                "scfd1\t1\t300\t1\tW\tctg1\t1\t300\t+",
                "scfd1\t301\t1000\t2\tW\tctg1\t301\t1000\t+",
            ],
        ),
        (
            "scfd2\t1\t1000\t1\tW\tctg2\t1\t1000\t-",
            [
                "scfd2\t1\t300\t1\tW\tctg2\t701\t1000\t-",
                "scfd2\t301\t1000\t2\tW\tctg2\t1\t700\t-",
            ],
        ),
    ],
)
def test_split_contig(in_agp_string, out_agp_strings):
    """Test splitting a contig

    This is similar to the pair of tests provided by @songtaogui
    (see test_split_scaffold earlier in this file), but goes a little
    bit deeper down into the mechanics of the program.
    """
    agp_row = AgpRow(in_agp_string)
    broken_rows = split_contig(agp_row, [300])

    for split_row, expected_split_row_string in zip(broken_rows, out_agp_strings):
        assert split_row == AgpRow(expected_split_row_string)


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
