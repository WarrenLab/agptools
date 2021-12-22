from io import StringIO

import pytest

from agp.bed import ParsingError, read


def test_read_bed():
    bed_text = (
        "scaffold_11\t1\t68766634\n"
        "scaffold_11\t68767135\t88214406\n"
        "scaffold_13\n"
        "scaffold_12\t70344537\t87137069"
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
