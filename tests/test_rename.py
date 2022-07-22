from os.path import dirname, join
from unittest.mock import patch

import pytest

from agp import open_agp
from agp.agptools import main
from agp.rename import FormatError, ScaffoldNotFoundError, renaming_file_type, run


def test_parse_renaming_file(tmp_path):
    with open(tmp_path / "renames.tsv", "w") as rename_file:
        print(
            "scaffold_1\tchr1\t+\n\nscaffold_2\tchr2\t-\nscaffold_3\tchr3",
            file=rename_file,
        )

    renames = renaming_file_type(tmp_path / "renames.tsv")
    assert renames["scaffold_1"] == ("chr1", "+")
    assert renames["scaffold_2"] == ("chr2", "-")
    assert renames["scaffold_3"] == ("chr3", "+")


def test_format_error_not_enough_columns(tmp_path):
    with open(tmp_path / "renames.tsv", "w") as rename_file:
        print(
            "scaffold_1\tchr1\t+\nscaffold_2\nscaffold_3\tchr3",
            file=rename_file,
        )

    with pytest.raises(FormatError) as error:
        renaming_file_type(tmp_path / "renames.tsv")
    assert "Line 2" in str(error.value)
    assert "columns" in str(error.value)


def test_format_error_bad_orientation(tmp_path):
    with open(tmp_path / "renames.tsv", "w") as rename_file:
        print(
            "scaffold_1\tchr1\t+\nscaffold_2\tchr2\t?\nscaffold_3\tchr3",
            file=rename_file,
        )

    with pytest.raises(FormatError) as error:
        renaming_file_type(tmp_path / "renames.tsv")
    assert "Line 2" in str(error.value)
    assert "orientation" in str(error.value)


@pytest.mark.parametrize(
    "rename_map, testfile_name",
    [
        ({"scaffold_17": ("chr1", "+"), "scaffold_18": ("chr2", "-")}, "renamed_1.agp"),
        ({"scaffold_16": ("chr1", "+"), "scaffold_17": ("chr2", "-")}, "renamed_2.agp"),
        ({"scaffold_16": ("chr1", "+"), "scaffold_18": ("chr2", "-")}, "renamed_3.agp"),
    ],
)
def test_run_rename(rename_map, testfile_name, tmp_path):
    run(
        rename_map,
        open(tmp_path / "test_out.agp", "w"),
        open_agp(join(dirname(__file__), "data", "test.agp")),
    )

    correct_agp = open_agp(join(dirname(__file__), "data", testfile_name))
    test_agp = open_agp(tmp_path / "test_out.agp")
    for correct_line, test_line in zip(correct_agp, test_agp):
        assert correct_line == test_line


@pytest.mark.parametrize(
    "rename_map, testfile_name",
    [
        ({"scaffold_17": ("chr1", "+"), "scaffold_18": ("chr2", "-")}, "renamed_1.agp"),
        ({"scaffold_16": ("chr1", "+"), "scaffold_17": ("chr2", "-")}, "renamed_2.agp"),
        ({"scaffold_16": ("chr1", "+"), "scaffold_18": ("chr2", "-")}, "renamed_3.agp"),
    ],
)
def test_rename_main(tmp_path, rename_map, testfile_name):
    with open(tmp_path / "scaffolds_to_rename.tsv", "w") as rename_file:
        for old_name, new_name_and_orientation in rename_map.items():
            print(
                "\t".join(
                    [old_name, new_name_and_orientation[0], new_name_and_orientation[1]]
                ),
                file=rename_file,
            )

    with patch(
        "sys.argv",
        [
            "agptools",
            "rename",
            join(tmp_path, "scaffolds_to_rename.tsv"),
            join(dirname(__file__), "data", "test.agp"),
            "-o",
            join(tmp_path, "test_out.agp"),
        ],
    ):
        main()

    correct_agp = open_agp(join(dirname(__file__), "data", testfile_name))
    test_agp = open_agp(tmp_path / "test_out.agp")
    for correct_line, test_line in zip(correct_agp, test_agp):
        assert correct_line == test_line


def test_rename_scaffold_not_found_error(tmp_path):
    with pytest.raises(ScaffoldNotFoundError) as error:
        run(
            {"scaffold_16": ("chr5", "-"), "scaffold_20": ("chr1", "-")},
            open(tmp_path / "test_out.agp", "w"),
            open_agp(join(dirname(__file__), "data", "test.agp")),
        )
    assert "scaffold_20" in str(error.value)
