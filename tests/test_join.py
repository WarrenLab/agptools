from os.path import dirname, join

import pytest

from agp import open_agp
from agp.join import (
    BadSequenceNameError,
    JoinGroup,
    ScaffoldNotFoundError,
    ScaffoldUsedTwiceError,
    joins_type,
    make_superscaffold_name,
    run,
)


def test_parse_joins(tmp_path):
    with open(tmp_path / "joins.txt", "w") as joins_file:
        print(
            "scaffold_1,-scaffold_2,+scaffold_3\nscaffold_4,-scaffold_5",
            file=joins_file,
        )

    joins = joins_type(tmp_path / "joins.txt")
    assert joins == [
        JoinGroup(["scaffold_1", "-scaffold_2", "+scaffold_3"]),
        JoinGroup(["scaffold_4", "-scaffold_5"]),
    ]


def test_parse_named_joins(tmp_path):
    with open(tmp_path / "joins.txt", "w") as joins_file:
        print(
            "scaffold_1,-scaffold_2\tchr1\nscaffold_4,-scaffold_5",
            file=joins_file,
        )
    joins = joins_type(tmp_path / "joins.txt")
    assert joins == [
        JoinGroup(["scaffold_1", "-scaffold_2"], name="chr1"),
        JoinGroup(["scaffold_4", "-scaffold_5"]),
    ]


@pytest.mark.parametrize(
    "joins_txt",
    [
        ("scaffold_1,-scaffold_2,+scaffold_3\nscaffold_3,+scaffold_5"),
        ("scaffold_1,-scaffold_2\n-scaffold_3,scaffold_2,scaffold_3"),
        ("scaffold_1,-scaffold_2,+scaffold_3\n+scaffold_3,+scaffold_5"),
        ("scaffold_1,-scaffold_2,+scaffold_3\n-scaffold_3,+scaffold_5"),
    ],
)
def test_scaffold_used_twice_error(joins_txt, tmp_path):
    with open(tmp_path / "joins.txt", "w") as joins_file:
        print(joins_txt, file=joins_file)
    with pytest.raises(ScaffoldUsedTwiceError) as error:
        joins_type(tmp_path / "joins.txt")
    assert "scaffold_3" in str(error.value)


def test_bad_sequence_name_error(tmp_path):
    with open(tmp_path / "joins.txt", "w") as joins_file:
        print("scaffold_1,scaffold_2\tchr 1", file=joins_file)
    with pytest.raises(BadSequenceNameError) as error:
        joins_type(tmp_path / "joins.txt")
    assert "chr 1" in str(error.value)


@pytest.mark.parametrize(
    "old_names,new_name",
    [
        (["scaffold_3", "scaffold_4", "scaffold_5"], "scaffold_3p4p5"),
        (["contig2", "scaffold_4", "chr1"], "contig2pscaffold_4pchr1"),
        (["scaffold_1", "contig_2", "scaffold_3"], "scaffold_1pcontig_2pscaffold_3"),
    ],
)
def test_make_superscaffold_name(old_names, new_name):
    assert make_superscaffold_name(old_names) == new_name


# TODO would be nice to implement, but not all that important
def test_join_scaffolds():
    pass


@pytest.mark.parametrize(
    "join_txt, correct_agp_filename",
    [
        ("scaffold_17,-scaffold_16", "test_join_1.agp"),
        ("+scaffold_18,-scaffold_16\n-scaffold_17\tchr2", "test_join_2.agp"),
        ("-scaffold_16,+scaffold_17\tchr2\n-scaffold_18", "test_join_3.agp"),
    ],
)
def test_run_join(join_txt, correct_agp_filename, tmp_path):
    with open(tmp_path / "joins.txt", "w") as joins_file:
        print(join_txt, file=joins_file)

    agp_in_path = join(dirname(__file__), "data", "test.agp")
    with open(tmp_path / "test_out.agp", "w") as outfile:
        run(joins_type(tmp_path / "joins.txt"), outfile, open_agp(agp_in_path))

    with open(tmp_path / "test_out.agp", "r") as test_agp, open(
        join(dirname(__file__), "data", correct_agp_filename)
    ) as correct_agp:
        for line1, line2 in zip(test_agp, correct_agp):
            assert line1 == line2


def test_scaffold_not_found(tmp_path):
    with open(tmp_path / "joins.txt", "w") as joins_file:
        print("scaffold_10,scaffold_17", file=joins_file)

    agp_in_path = join(dirname(__file__), "data", "test.agp")
    with pytest.raises(ScaffoldNotFoundError) as error:
        with open(tmp_path / "test_out.agp", "w") as outfile:
            run(joins_type(tmp_path / "joins.txt"), outfile, open_agp(agp_in_path))
    assert "scaffold_10" in str(error.value)
