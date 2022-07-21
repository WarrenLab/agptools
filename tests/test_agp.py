import os.path

import pytest

import agp

gap_row_string = "scaffold_18\t8679231\t8679730\t28\tN\t500\tscaffold\tyes\tna"
contig_row_string = (
    "scaffold_18\t7309250\t8679230\t27\tW\ttig00005023|arrow|arrow\t1\t1369981\t+"
)


def test_parse_gap_row():
    """Test that an example gap row is parsed correctly."""
    row = agp.AgpRow(gap_row_string)
    assert row.object == "scaffold_18"
    assert row.object_beg == 8679231
    assert row.object_end == 8679730
    assert row.part_number == 28
    assert row.component_type == "N"
    assert row.is_gap
    assert row.gap_length == 500
    assert row.gap_type == "scaffold"
    assert row.linkage == "yes"
    assert row.linkage_evidence == "na"


def test_parse_contig_row():
    """Test that an example contig row is parsed correctly."""
    row = agp.AgpRow(contig_row_string)
    assert row.object == "scaffold_18"
    assert row.object_beg == 7309250
    assert row.object_end == 8679230
    assert row.part_number == 27
    assert row.component_type == "W"
    assert not row.is_gap
    assert row.component_id == "tig00005023|arrow|arrow"
    assert row.component_beg == 1
    assert row.component_end == 1369981
    assert row.orientation == "+"


def test_parse_misformatted_row():
    """Test that a misformatted row raises the correct error."""
    with pytest.raises(agp.AgpFormatError) as excinfo:
        agp.AgpRow("scaffold_18\t8679q231\t8679730\t28\tN\t500\tscaffold\tyes\tna")
    assert "Invalid AGP:" in str(excinfo.value)


@pytest.mark.parametrize(
    "row1_string, row2_string, equality",
    [
        (gap_row_string, gap_row_string, True),
        (contig_row_string, contig_row_string, True),
        (gap_row_string, contig_row_string, False),
    ],
)
def test_row_equality(row1_string, row2_string, equality):
    """Test that row equality comparison function works as expected"""
    row1, row2 = agp.AgpRow(row1_string), agp.AgpRow(row2_string)
    assert (row1 == row2) == equality


def test_gap_row_to_str():
    """Test that a gap row can print back to itself"""
    assert str(agp.AgpRow(gap_row_string)) == gap_row_string


def test_contig_row_to_str():
    """Test that a contig row string can print back to itself"""
    assert str(agp.AgpRow(contig_row_string)) == contig_row_string


def test_contains():
    """Test that AgpRow.contains() works as expected"""
    row = agp.AgpRow(contig_row_string)
    assert row.contains(8000000)
    assert not row.contains(9000000)


def test_init_gaprow():
    row = agp.GapRow("scaffold_10", 1000000, 2000000, 2)
    assert (
        str(row)
        == "scaffold_10\t1000000\t2000000\t2\tN\t500\tscaffold\tyes\tpaired-end"
    )


def test_read():
    """Test reading a whole agp file with a comment on the first line"""
    agp_path = os.path.join(os.path.dirname(__file__), "data", "test.agp")
    for row_number, row in enumerate(agp.open_agp(agp_path)):
        if row_number == 0:
            assert agp.is_string(row)
        else:
            assert isinstance(row, agp.AgpRow)
