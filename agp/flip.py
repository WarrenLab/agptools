#!/usr/bin/env python3
"""
Given an AGP file describing an assembly and a list of chromosomes to
flip, output a new AGP file with those chromosomes reoriented in
reverse-complement.
"""

from agp import bed


def reverse_rows(rows):
    """
    Given a list of AGP rows, return a new set of rows assembling
    the reverse complement of the input rows.

    Args:
        rows (list): an ordered list of AgpRow objects

    Returns:
        reverse_chromosome_rows (list): the input, except in reverse
            order, with the rows modified to be reverse-complements
    """
    # first, get the first and last part numbers and the start and end
    # position of the entire segment being flipped so we can count
    # backwards from there.
    first_part_number = rows[0].part_number
    last_part_number = rows[-1].part_number
    first_part_start = rows[0].object_beg
    last_part_end = rows[-1].object_end

    reversed_rows = []
    for row in reversed(rows):
        new_start = first_part_start + last_part_end - row.object_end
        new_end = first_part_start + last_part_end - row.object_beg
        row.object_beg, row.object_end = new_start, new_end

        row.part_number = first_part_number + last_part_number - row.part_number

        if not row.is_gap:
            if row.orientation == "+":
                row.orientation = "-"
            elif row.orientation == "-":
                row.orientation = "+"
            # if the orientation is something else, leave it be

        reversed_rows.append(row)

    return reversed_rows


def flip(agp_rows, ranges_to_flip):
    """
    Reverse-complements all rows of an AGP file that fall within a list
    of ranges.

    Args:
        agp_rows (Iterable(AgpRow)): all rows in an AGP file
        ranges_to_flip (list(BedRange)): a list of ranges in the AGP
            file to reverse-complement. Takes genomic (bed)
            coordinates, not array coordinates.

    Returns:
        agp_rows (list(AgpRow)): the same rows that were input, except
            the rows in the specified range are now reverse-
            complemented
    """
    # agp_rows is a generator, but we're going to have to store the
    # whole agp in memory
    agp_rows = list(agp_rows)

    for bed_range in ranges_to_flip:
        # list of tuples (i, agp_row) of rows in bed_range
        rows_to_reverse = []
        for i, row in enumerate(agp_rows):
            if not isinstance(row, str) and row.object == bed_range.chrom:
                # if no range in seq specified, flip the whole sequence
                if bed_range.start is None:
                    rows_to_reverse.append((i, row))
                # otherwise, flip only the part within range
                elif (
                    row.object_beg >= bed_range.start
                    and row.object_end <= bed_range.end
                ):
                    rows_to_reverse.append((i, row))
                # check for bad ranges (i.e., ones that only partially
                # contain a component)
                elif row.contains(bed_range.start) or row.contains(bed_range.end):
                    raise bed.BadRangeError(bed_range)

        # unzip and flip those rows we've collected
        if len(rows_to_reverse) == 0:
            raise bed.EmptyRangeError(bed_range)
        indices, rows = map(list, zip(*rows_to_reverse))
        reversed_rows = reverse_rows(rows)

        # replace the original rows with the reversed rows
        for i, row in zip(indices, reversed_rows):
            agp_rows[i] = row

    return agp_rows


def run(segments_to_flip, outfile, agp_rows):
    for row in flip(agp_rows, segments_to_flip):
        print(row, file=outfile)
