"""
Functions for splitting a scaffold into subscaffolds at gaps.
"""

import re
from copy import deepcopy
from typing import Dict, Iterable, Iterator, List, TextIO, Union

from . import AgpRow


class ParsingError(Exception):
    """Raised when breakpoints file is misformatted."""

    pass


empty_line_regex = re.compile(r"^\s*$")


def breakpoints_type(filename: str) -> Dict[str, List[int]]:
    """
    Argparse type function for breakpoints file: first column is the
    scaffold name; second column is a comma-separated list of
    locations within gaps where scaffold should be broken.

    Args:
        filename: path to the breakpoints file

    Returns:
        breakpoints: a dict mapping scaffold name to a list
            of breakpoints (int) on that scaffold

    Raises:
        FileNotFoundError: if `filename` does not point to readable file
        ParsingError: if there is an error in the input file's format or
            content
    """
    breakpoints = {}
    with open(filename) as breakpoints_file:
        for i, line in enumerate(breakpoints_file):
            if not empty_line_regex.match(line):
                splits = line.strip().split("\t")
                try:
                    if splits[0] in breakpoints:
                        raise ParsingError(
                            f"{splits[0]} specified multiple times in breakpoints file"
                        )
                    breakpoints[splits[0]] = list(map(int, splits[1].split(",")))
                except (ValueError, IndexError):
                    raise ParsingError(f"Cannot parse line {i} of breakpoints: {line}")
    return breakpoints


def unoffset_rows(new_scaffold_name: str, rows: List[AgpRow]) -> List[AgpRow]:
    """
    Modifies some AGP rows so that they can be their own standalone
    scaffold. This requires changing their object names to a new
    scaffold name, and changing the part numbers and coordinates such
    that the first row starts with 1 and the rest follow.

    Args:
        new_scaffold_name: name for the new scaffold which will
            replace all 'object' fields
        rows: rows to modify so that they function as a
            standalone scaffold together. The first row will be used to
            calculate offsets.

    Returns:
        out_rows: input rows, but with all 'object'
            fields replaced with new_scaffold_name, and all positions
            and part numbers modified so that the first row is the
            beginning of a new scaffold.
    """
    position_offset = rows[0].object_beg - 1
    out_rows = []
    for i, row in enumerate(rows):
        row.object = new_scaffold_name
        row.object_beg -= position_offset
        row.object_end -= position_offset
        row.part_number = i + 1
        out_rows.append(row)
    return out_rows


def split_contig(contig_row: AgpRow, breakpoints: List[int]) -> List[AgpRow]:
    """
    Splits a single row containing a contig into multiple rows,
    each containing a piece of the contig.

    Args:
        contig_row: a single row to be split
        breakpoints: positions where contig should be split,
            in object coordinates, *not* component coordinates. The left
            part of the split includes the breakpoint: e.g., splitting a
            contig of length 100 at 43 will make two new contigs: one
            from 1-43 and the other from 44-100.

    Returns:
        a list of rows of length len(breakpoints) + 1, where each row is
        a piece of the contig after splitting. Note that while the new
        coordinates are correct in the coordinates of the original
        scaffold they came from, the object coordinates still need to be
        offset to match the new scaffold they will end up in.

    Examples:
        >>> import agp
        >>> r = agp.AgpRow('\\t'.join(map(str, ['scaf', 501, 1000, 5, 'W',
        ...                                    'ctg', 1, 500, '+'])))
        >>> [str(x).split('\\t') for x in split_contig(r, [750, 867])]
        [['scaf', '501', '750', '5', 'W', 'ctg', '1', '250', '+'],
         ['scaf', '751', '867', '6', 'W', 'ctg', '251', '367', '+'],
         ['scaf', '868', '1000', '7', 'W', 'ctg', '368', '500', '+']]
    """
    rows = [contig_row]
    for this_breakpoint in sorted(breakpoints):
        left_part = deepcopy(rows.pop())
        right_part = deepcopy(left_part)

        left_part.object_end = this_breakpoint
        right_part.object_beg = this_breakpoint + 1
        right_part.part_number += 1

        if contig_row.orientation == "+":
            left_part.component_end = left_part.component_beg + (
                this_breakpoint - left_part.object_beg
            )
            right_part.component_beg = left_part.component_end + 1
        else:
            left_part.component_beg = (
                left_part.object_beg + left_part.component_end - this_breakpoint
            )
            right_part.component_end = left_part.component_beg - 1

        rows += [left_part, right_part]
    return rows


def convert_rows(rows: List[AgpRow], subscaffold_counter: int) -> List[AgpRow]:
    """
    Converts rows that are part of a scaffold into their own standalone
    scaffold. Changes the positions and part numbers so that the new
    scaffold starts at 1, and names the new scaffold after the old
    scaffold, but with '.{subscaffold_counter}' at the end.

    Args:
        rows: rows to turn into their own scaffold.
        subscaffold_counter: suffix to add to the old scaffold
            name in order to turn it into the new scaffold name.

    Returns:
        the input rows, but with object names, positions, and part
        numbers changed so that they now function as a standalone
        scaffold
    """
    new_scaffold_name = "{}.{}".format(rows[0].object, subscaffold_counter)
    return unoffset_rows(new_scaffold_name, rows)


def split_scaffold(
    scaffold_rows: Iterable[AgpRow], breakpoints: List[int]
) -> List[AgpRow]:
    """
    Splits a scaffold at specified breakpoints.

    Args:
        scaffold_rows: all the rows for a given scaffold in an AGP file
        breakpoints: a list of locations where scaffold should be
            broken. All locations are specified in genomic coordinates
            and must be within the boundaries of a gap.

    Returns:
        broken_rows: rows of the input scaffold broken
            into len(breakpoints)+1 sub-scaffolds, with the gaps
            containing the breakpoints removed
    """
    out_rows = []
    rows_this_subscaffold: List[AgpRow] = []
    subscaffold_counter = 1
    for row in scaffold_rows:
        if any(map(row.contains, breakpoints)):
            # if the breakpoint is within a gap, our job is simple:
            # just forget about the gap row, output the previous
            # subscaffold, and start a new subscaffold
            if row.is_gap:
                out_rows += convert_rows(rows_this_subscaffold, subscaffold_counter)

                rows_this_subscaffold = []
                subscaffold_counter += 1
            # if the breakpoint is not within a gap, we need to actually
            # break a contig into pieces
            else:
                # split the contig into two or more rows
                contig_rows = split_contig(row, list(filter(row.contains, breakpoints)))

                # the first row goes at the end of the current scaffold
                rows_this_subscaffold.append(contig_rows[0])
                del contig_rows[0]
                out_rows += convert_rows(rows_this_subscaffold, subscaffold_counter)
                subscaffold_counter += 1

                # the last row goes at the beginning of the next
                # scaffold
                rows_this_subscaffold = [contig_rows.pop()]

                # if there are any rows in between, they each get their
                # own subscaffold
                for contig_part in contig_rows:
                    out_rows += convert_rows([contig_part], subscaffold_counter)
                    subscaffold_counter += 1

        else:  # only add this row if there are no breakpoints in it
            rows_this_subscaffold.append(row)

    out_rows += convert_rows(rows_this_subscaffold, subscaffold_counter)

    return out_rows


def run(
    breakpoints: Dict[str, List[int]],
    outfile: TextIO,
    agp_infile: Iterator[Union[str, AgpRow]],
):
    rows_this_scaffold: List[AgpRow] = []  # list of all agp rows in current scaffold
    for row in agp_infile:
        if isinstance(row, str):  # print out comment rows as-is
            print(row, file=outfile)
            continue

        # if we're on a new scaffold, do any necessary modification to
        # the previous scaffold, print it out, and clear the buffer
        if rows_this_scaffold and rows_this_scaffold[0].object != row.object:
            if rows_this_scaffold[0].object in breakpoints:
                rows_this_scaffold = split_scaffold(
                    rows_this_scaffold,
                    breakpoints[rows_this_scaffold[0].object],
                )
            for r in rows_this_scaffold:
                print(r, file=outfile)
            rows_this_scaffold = []

        rows_this_scaffold.append(row)

    if rows_this_scaffold[0].object in breakpoints:
        rows_this_scaffold = split_scaffold(
            rows_this_scaffold,
            breakpoints[rows_this_scaffold[0].object],
        )
    for r in rows_this_scaffold:
        print(r, file=outfile)
