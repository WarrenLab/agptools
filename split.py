"""
Functions for splitting a scaffold into subscaffolds at gaps.
"""

from copy import deepcopy

class BreakpointError(Exception):
    def __init__(self, scaffold_name, bad_breakpoints):
        self.scaffold_name = scaffold_name
        self.bad_breakpoints = bad_breakpoints

    def __str__(self):
        return "Breakpoint(s) {} in {} are not in a gap!".format(
                self.scaffold_name,
                self.bad_breakpoints,
                )


def breakpoints_type(filename):
    """
    Argparse type function for breakpoints file: first column is the
    scaffold name; second column is a comma-separated list of
    locations within gaps where scaffold should be broken.

    Args:
        filename (str): path to the breakpoints file

    Returns:
        breakpoints (dict): a dict mapping scaffold name to a list
            of breakpoints (int) on that scaffold
    """
    breakpoints = {}
    with open(filename) as breakpoints_file:
        for line in breakpoints_file:
            splits = line.strip().split()
            breakpoints[splits[0]] = list(map(int, splits[1].split(',')))
    return breakpoints


def unoffset_rows(new_scaffold_name, rows):
    """
    Modifies some AGP rows so that they can be their own standalone
    scaffold. This requires changing their object names to a new
    scaffold name, and changing the part numbers and coordinates such
    that the first row starts with 1 and the rest follow.

    Args:
        new_scaffold_name (str): name for the new scaffold which will
            replace all 'object' fields
        rows (list(AgpRow)): rows to modify so that they function as a
            standalone scaffold together. The first row will be used to
            calculate offsets.

    Returns:
        out_rows (list(AgpRow)): input rows, but with all 'object'
            fields replaced with new_scaffold_name, and all positions
            and part numbers modified so that the first row is the
            beginning of a new scaffold.
    """
    position_offset = rows[0].object_beg - 1
    part_number_offset = rows[0].part_number - 1
    out_rows = []
    for row in rows:
        row.object = new_scaffold_name
        row.object_beg -= position_offset
        row.object_end -= position_offset
        row.part_number -= part_number_offset
        out_rows.append(row)
    return out_rows


def split_contig(contig_row, breakpoints):
    """
    Splits a single row containing a contig into multiple rows,
    each containing a piece of the contig.

    >>> import agp
    >>> r = agp.AgpRow('\\t'.join(map(str, ['scaf', 501, 1000, 5, 'W',
    ...                                    'ctg', 1, 500, '+'])))
    >>> list(map(lambda x: str(x).split('\\t'),
    ...          split_contig(r, [750, 867])))
    [['scaf', '501', '750', '5', 'W', 'ctg', '1', '250', '+'],
     ['scaf', '751', '867', '6', 'W', 'ctg', '251', '367', '+'],
     ['scaf', '868', '1000', '7', 'W', 'ctg', '368', '500', '+']]

    Args:
        contig_row (AgpRow): a single row to be split
        breakpoints (list(int)): positions where contig should be split,
            in object coordinates, *not* component coordinates. The left
            part of the split includes the breakpoint: e.g., splitting a
            contig of length 100 at 43 will make two new contigs: one
            from 1-43 and the other from 44-100.
    """
    rows = [contig_row]
    for breakpoint in sorted(breakpoints):
        left_part = deepcopy(rows.pop())
        right_part = deepcopy(left_part)

        left_part.object_end = breakpoint
        right_part.object_beg = breakpoint + 1
        right_part.part_number += 1
        left_part.component_end = (left_part.component_beg
                                   + (breakpoint - left_part.object_beg))
        right_part.component_beg = left_part.component_end + 1

        rows += [left_part, right_part]
    return rows


def split_scaffold(scaffold_rows, breakpoints):
    """
    Splits a scaffold at specified breakpoints.

    Args:
        scaffold_rows (list(AgpRow)): all the rows for a given scaffold
            in an AGP file
        breakpoints (list(int)): a list of locations where scaffold
            should be broken. All locations are specified in genomic
            coordinates and must be within the boundaries of a gap.

    Returns:
        broken_rows (list(AgpRow)): rows of the input scaffold broken
            into len(breakpoints)+1 sub-scaffolds, with the gaps
            containing the breakpoints removed
    """
    out_rows = []
    rows_this_subscaffold = []
    subscaffold_counter = 1
    for row in scaffold_rows:
        if any(map(row.contains, breakpoints)):
            if row.is_gap: # this is a good breakpoint
                new_scaffold_name = '{}.{}'.format(
                        rows_this_subscaffold[0].object,
                        subscaffold_counter,
                        )
                out_rows += unoffset_rows(
                        new_scaffold_name,
                        rows_this_subscaffold,
                        )

                rows_this_subscaffold = []
                subscaffold_counter += 1
            else: # this is a bad breakpoint
                raise BreakpointError(row.object[0],
                                      filter(row.contains, breakpoints))
        else: # only add this row if there are no breakpoints in it
            rows_this_subscaffold.append(row)

    new_scaffold_name = '{}.{}'.format(
            rows_this_subscaffold[0].object,
            subscaffold_counter,
            )
    out_rows += unoffset_rows(
            new_scaffold_name,
            rows_this_subscaffold,
            )

    return out_rows


def run(breakpoints, outfile, agp):
    rows_this_scaffold = [] # list of all agp rows in current scaffold
    for row in agp:
        if isinstance(row, str): # print out comment rows as-is
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
            for r in rows_this_scaffold: print(r, file=outfile)
            rows_this_scaffold = []

        rows_this_scaffold.append(row)

    if rows_this_scaffold[0].object in breakpoints:
        rows_this_scaffold = split_scaffold(
            rows_this_scaffold,
            breakpoints[rows_this_scaffold[0].object],
        )
    for r in rows_this_scaffold: print(r, file=outfile)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
