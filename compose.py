"""
Functions for composing two AGPs together.
"""

import collections
from functools import partial
from itertools import filterfalse

class BrokenComponentError(Exception):
    def __init__(self, outer_row, inner_row):
        self.outer_row, self.inner_row = outer_row, inner_row

    def __str__(self):
        return """{}:{}-{} does not fully contain {}:{}-{}, and breaking
                  components in this module is not supported at this
                  time.""".format(self.outer_row.object,
                                  self.outer_row.object_beg,
                                  self.outer_row.object_end,
                                  self.inner_row.object,
                                  self.inner_row.object_beg,
                                  self.inner_row.object_end)


def contains(outer_row, inner_row):
    """
    Determines whether outer_row fully contains inner_row. For example,
    if outer_row specifies that the current object contains the range
    500-1000 of a component1, and inner_row specifies how to assemble
    the range 600-800 of this component, then outer_row contains
    inner_row. Because breaking components is not supported, right now,
    this function raises a BrokenComponentError if the outer row
    partially contains the inner row.
    """
    if (outer_row.component_beg <= inner_row.object_beg
            and outer_row.component_end >= inner_row.object_end):
        return True
    elif (outer_row.component_beg <= inner_row.object_beg
          or outer_row.component_end >= inner_row.object_end):
        raise BrokenComponentError(outer_row, inner_row)
    else:
        return False


def run(outer_agp, inner_agp, outfile):
    """runs the compose module of agptools"""
    # read the inner agp into a dictionary so that we can look up the
    # definition of an outer component later
    inner_dict = collections.defaultdict(list)
    for inner_row in filter(lambda r: not isinstance(r, str), inner_agp):
        inner_dict[inner_row.object].append(row)

    previous_scaffold = None
    component_counter = 1
    for outer_row in outer_agp:
        # print the row as is if it's a comment
        if isinstance(outer_row, str):
            print(outer_row, file=outfile)
            continue

        # reset the component counter if we're on a new scaffold
        if previous_scaffold and previous_scaffold != outer_row.object:
            component_counter = 1

        previous_scaffold = outer_row.object

        # for gaps, we only need to change the part number
        if outer_row.is_gap:
            outer_row.part_number = component_counter
            print(outer_row, file=outfile)
            component_counter += 1
        # for non-gaps, we need to replace this single row with all of
        # the component rows, modified as necessary
        else:
            # check to make sure this component is actually defined
            if outer_row.component_id not in inner_dict:
                raise Exception # TODO specify!

            # make a list of all rows in this component, within range
            inner_rows = inner_dict[outer_row.component_id]
            component_rows = filter(partial(contains, outer_row), inner_rows)

            # remove anything we're using now from inner_dict, as we
            # will later be printing all components that did not get
            # used
            inner_dict[outer_row.component_id][:] = filterfalse(
                    partial(contains, outer_row), inner_rows)

            # reverse complement everything if necessary
            if outer_row.orientation == '-':
                component_rows = flip.reverse_rows(component_rows)

            # TODO fix offsets


