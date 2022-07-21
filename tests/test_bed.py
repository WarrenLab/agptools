from io import StringIO

import pytest

from agp.bed import BedRange, ParsingError, read


def test_read_bed():
    bed_text = (
        "scaffold_11\t1\t68766634\n"
        "scaffold_11\t68767135\t88214406\t+\n"
        "scaffold_13\n"
        "scaffold_12\t70344537\t87137069\n"
        "scaffold_14\t1\t12309\t-\tblahblah\tblah\n"
    )
    for i, bed_range in enumerate(read(StringIO(bed_text))):
        assert str(bed_range.chrom) in bed_text
        if bed_range.start and bed_range.end:
            assert str(bed_range.start) in bed_text
            assert str(bed_range.end) in bed_text

        if i == 1:
            assert bed_range.chrom == "scaffold_11"
            assert bed_range.start == 68767135
            assert bed_range.end == 88214406
            assert bed_range.strand == "+"

        if i == 4:
            assert bed_range.extra_fields[0] == "blahblah"
            assert bed_range.extra_fields[1] == "blah"


def test_read_bad_bed():
    bed_text = (
        "scaffold_11\t1\t68766634\n"
        "scaffold_11\t68767135   88214406\n"
        "scaffold_12\t70344537\t87137069"
    )
    with pytest.raises(ParsingError) as error:
        for bed_range in read(StringIO(bed_text)):
            pass

    assert "Line 2" in str(error.value)


@pytest.mark.parametrize(
    "bed_line",
    [
        "scaffold_11\t68766634",
        "scaffold_11\t68767135   88214406",
        "scaffold_12\t703q4537\t87137069",
    ],
)
def test_read_bad_bed_line(bed_line):
    with pytest.raises(ParsingError):
        for bed_range in read(StringIO(bed_line)):
            pass


@pytest.mark.parametrize(
    "bed_line, bed_line_str",
    [
        (BedRange("scaffold_11", 1, 68766634), "scaffold_11\t1\t68766634"),
        (BedRange("scaffold_11"), "scaffold_11"),
        (
            BedRange("scaffold_11", 1, 2, "+", ["blah", "blee", "bleh"]),
            "scaffold_11\t1\t2\t+\tblah\tblee\tbleh",
        ),
    ],
)
def test_bed_to_string(bed_line, bed_line_str):
    assert str(bed_line) == bed_line_str
