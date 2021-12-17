import pytest

from agp.join import (
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
def test_parse_bad_joins(joins_txt, tmp_path):
    with open(tmp_path / "joins.txt", "w") as joins_file:
        print(joins_txt, file=joins_file)
    with pytest.raises(ScaffoldUsedTwiceError) as error:
        joins_type(tmp_path / "joins.txt")
    assert "scaffold_3" in str(error.value)


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


def test_join_scaffolds():
    pass


def test_run_join():
    pass
