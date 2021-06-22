"""
Functions for composing two AGPs together.
"""

import collections
from itertools import chain, filterfalse

import agp
from agp import flip


class BrokenComponentError(Exception):
    def __init__(self, outer_row, inner_row):
        self.outer_row, self.inner_row = outer_row, inner_row

    def __str__(self):
        return """{}:{}-{} does not fully contain {}:{}-{}, and breaking
                  components in this module is not supported at this
                  time.""".format(
            self.outer_row.object,
            self.outer_row.object_beg,
            self.outer_row.object_end,
            self.inner_row.object,
            self.inner_row.object_beg,
            self.inner_row.object_end,
        )


class ComponentNotDefinedError(Exception):
    def __init__(self, component_name):
        self.component_name = component_name

    def __str__(self):
        return 'Component "{}" not defined.'.format(self.component_name)


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
    if (
        outer_row.component_beg <= inner_row.object_beg
        and outer_row.component_end >= inner_row.object_end
    ):
        return True
    elif (
        outer_row.component_beg <= inner_row.object_beg
        or outer_row.component_end >= inner_row.object_end
    ):
        raise BrokenComponentError(outer_row, inner_row)
    else:
        return False


def run(outer_agp, inner_agp, outfile, print_unused=False):
    """runs the compose module of agptools"""
    # read the inner agp into a dictionary so that we can look up the
    # definition of an outer component later. inner_dict maps object
    # names to a list of all AgpRow instances with that object name
    inner_dict = collections.defaultdict(list)
    for inner_row in filterfalse(agp.is_string, inner_agp):
        inner_dict[inner_row.object].append(inner_row)

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
                raise ComponentNotDefinedError(outer_row.component_id)

            # make a list of all rows in this component and make sure
            # the object is using the whole component
            inner_rows = inner_dict[outer_row.component_id]
            del inner_dict[outer_row.component_id]
            if (
                outer_row.component_beg != 1
                or outer_row.component_end != inner_rows[-1].object_end
            ):
                raise BrokenComponentError(outer_row, inner_rows[-1])

            # reverse complement everything if necessary
            if outer_row.orientation == "-":
                inner_rows = flip.reverse_rows(inner_rows)

            # fix offsets and print rows
            for inner_row in inner_rows:
                inner_row.object = outer_row.object
                inner_row.object_beg += outer_row.object_beg - 1
                inner_row.object_end += outer_row.object_beg - 1
                inner_row.part_number = component_counter
                component_counter += 1
                print(inner_row, file=outfile)

    # at the end, print all the objects in the inner AGP that were not
    # assembled into anything in the outer AGP
    if print_unused:
        for row in chain.from_iterable(inner_dict.values()):
            print(row, file=outfile)
