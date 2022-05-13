"""
Functions for transforming genomic coordinates from contig-based to
scaffold-based.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Iterator, TextIO, Union

from agp import AgpRow
from agp.bed import BedRange


class NoSuchContigError(Exception):
    def __init__(self, contig_name):
        self.contig_name = contig_name

    def __str__(self):
        return 'Cannot find contig called "{}" in AGP'.format(self.contig_name)


class UnsupportedOperationError(Exception):
    """Error for trying to transform an unsupported style of bed row"""

    pass


class CoordinateNotFoundError(Exception):
    """Error for bad coordinates specified"""

    pass


class BadOrientationError(Exception):
    """Error for bad orientation (not +/-)"""

    pass


ContigDict = dict[str, list[AgpRow]]


def create_contig_dict(agp_in: Iterator[Union[str, AgpRow]]) -> ContigDict:
    """
    Load the agp file into a dictionary mapping contig name to a list
    of rows in which that contig is the component_id.

    Args:
        agp_in: an agp file to read into a dictionary

    Returns:
        a dict mapping contig name to a list of rows in which that
        contig is the component_id
    """
    contig_dict: ContigDict = defaultdict(list)
    for row in (r for r in agp_in if isinstance(r, AgpRow) and not r.is_gap):
        contig_dict[row.component_id].append(row)
    return contig_dict


def transform_single_position(position: int, agp_row: AgpRow) -> int:
    """Transform a position from component to object coordinates

    Given a position in component coordinates (e.g., a position on a
    contig), and the agp row containing that position, transform the
    position to object coordinates (e.g., the same position, but on
    a scaffold of which that contig is a part).

    Args:
        position: a position in component/contig coordinates
        agp_row: the row containing that position

    Returns:
        the same position, but in object/scaffold coordinates
    """

    # first, offset the position by where it is on this _row_. This
    # is necessary because not all rows contain an entire component/
    # contig. E.g., if this row contains positions 101-200 on the
    # contig, then position 150 on the contig is the 50th base in
    # the piece of the contig enclosed by this row
    position = position - agp_row.component_beg + 1
    if agp_row.orientation == "+":
        # then, offset the position by where the row is on the object/
        # scaffold
        position += agp_row.object_beg - 1
    elif agp_row.orientation == "-":
        position = agp_row.object_end - position + 1
    else:
        raise BadOrientationError(f"Orientation must be +/-: {agp_row}")
    return position


def find_agp_row(
    contig_name: str, coordinate_on_contig: int, contig_dict: ContigDict
) -> AgpRow:
    """Find the agp row containing a coordinate.

    Given a dictionary mapping contig name to a list of agp rows
    containing that contig, and a single position in contig coordinates,
    find and return the AGP row containing that position.

    Args:
        contig_name:
            the contig containing the coordinate
        coordinate_on_contig:
            the position to look for, in contig coordinates
        contig_dict:
            a dictionary mapping contig name to a list of agp rows
            containing that contig

    Returns:
        a single AGP row containing the requested position
    """
    contig_list = contig_dict[contig_name]
    if not contig_list:
        raise NoSuchContigError(contig_name)
    for agp_row in contig_list:
        if (
            agp_row.component_beg <= coordinate_on_contig
            and agp_row.component_end >= coordinate_on_contig
        ):
            return agp_row

    raise CoordinateNotFoundError(f"{contig_name}:{coordinate_on_contig}")


def transform_bed_row(bed_row: BedRange, contig_dict: ContigDict) -> BedRange:
    """Transform a bed row to scaffold coordinates

    Transform a single row of a bed file from contig coordinates to
    scaffold coordinates.

    Args:
        bed_row: the row to transform
        contig_dict:
            a dictionary mapping contig name to the agp row(s)
            containing that contig

    Returns:
        the same bed row, but in scaffold coordinates now
    """
    if bed_row.start is None or bed_row.end is None:
        raise UnsupportedOperationError(
            f"Transforming coordinateless bed row: {bed_row}"
        )

    # retrieve the AGP rows containing the bed range's start and end points
    agp_row_start = find_agp_row(bed_row.chrom, bed_row.start, contig_dict)
    agp_row_end = find_agp_row(bed_row.chrom, bed_row.end, contig_dict)

    # only transform bed ranges that are entirely within one scaffold
    if agp_row_start.object == agp_row_end.object:
        bed_row.chrom = agp_row_start.object
        # transform the start and end positions of the bed range to the
        # new coordinate system
        bed_row.start, bed_row.end = sorted(
            [
                transform_single_position(bed_row.start, agp_row_start),
                transform_single_position(bed_row.end, agp_row_end),
            ]
        )
        # flip the orientation of the bed row, if necessary
        if (
            agp_row_start == agp_row_end
            and bed_row.strand
            and agp_row_start.orientation == "-"
        ):
            if bed_row.strand == "+":
                bed_row.strand = "-"
            else:
                bed_row.strand = "+"

    return bed_row


def run(
    bed_in: Iterator[BedRange], agp_in: Iterator[Union[str, AgpRow]], bed_out: TextIO
):
    contig_dict = create_contig_dict(agp_in)
    for bed_row in bed_in:
        print(transform_bed_row(bed_row, contig_dict), file=bed_out)
