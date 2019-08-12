"""
Functions for joining scaffolds
"""

from itertools import chain
import re
import sys

import agp
import flip

scaffold_regex = re.compile(r'([a-zA-Z0-9]+)_([a-zA-Z0-9]+)')

def joins_type(filename):
    """
    argparse type function for file listing scaffold joins

    Args:
        filename (str): filename giving path to a file listing
            scaffolds to join. Each line should contain a comma-
            separated list of scaffolds to be joined into a single
            scaffold.

    Returns:
        joins (list(list(str))): list of join groups, where each join
            group is a list of the scaffolds in that group
    """
    joins = []
    with open(filename) as joins_file:
        for line in joins_file:
            joins.append(line.strip().split(','))
    return joins

def make_superscaffold_name(subscaffold_names):
    """
    Comes up with a nice superscaffold name based on the given
    subscaffold names. If the subscaffold names are all in the format
    '[prefix]_[suffix]', where prefix is the same for all subscaffolds,
    then the superscaffold name is '[prefix]_[suffix1]p[suffix2]p[etc.]'
    Otherwise, it the names of all scaffolds concatenated with 'p'.

    Args:
        subscaffold_names (list(str)): list of names of subscaffolds
            being combined into a superscaffold

    Returns:
        superscaffold_name (str): a name for the new scaffolds based
            on the subscaffold names
    """
    matches = list(map(scaffold_regex.search, subscaffold_names))
    if all(m.group(1) == matches[0].group(1) for m in matches):
        prefix = matches[0].group(1)
        suffix = 'p'.join([m.group(2) for m in matches])
        return '{}_{}'.format(prefix, suffix)
    else:
        return 'p'.join(subscaffold_names)


def join_scaffolds(superscaffold_rows, gap_size=500, gap_type='scaffold',
                   gap_evidence='na'):
    """
    Transforms agp rows from multiple scaffolds into agp rows for a
    single superscaffold containing all the given scaffolds in the
    specified order. Scaffolds should be oriented properly before
    using this function.

    Args:
        superscaffold_rows (list(list(AgpRow)): a list of lists of agp
            rows. Each sub-list is a list of all rows in a single
            scaffold. These scaffolds will be joined in the order given
        gap_size (int): length of the new gaps created by joining
            scaffolds together
        gap_type (str): what to call the new gaps in the agp file
        gap_evidence (str): evidence type for linkage across gap

    Returns:
        agp_rows (list(AgpRow)): a list of AgpRow instances containing
            the new scaffold specification.
    """
    # make a nice name for the new superscaffold we are creating
    subscaffold_names = map(lambda s: s[0].object, superscaffold_rows)
    superscaffold_name = make_superscaffold_name(subscaffold_names)

    # keeps track of what component number of the superscaffold we are
    # currently on
    component_number_counter = 1
    # keeps track of the ending position of the previous subscaffold,
    # so that we can add it as an offset to components of the current
    # subscaffold
    end_of_previous_scaffold = 0
    # list(AgpRow) containing the rows to be returned
    outrows = []
    # loop over the subscaffolds
    for i, this_scaffold_rows in enumerate(superscaffold_rows):
        # loop over the agp rows in this subscaffold
        for row in this_scaffold_rows:
            # update the current row and put it in the output list
            row.object = superscaffold_name
            row.part_number = component_number_counter
            component_number_counter += 1
            row.object_beg += end_of_previous_scaffold
            row.object_end += end_of_previous_scaffold
            outrows.append(row)
        end_of_previous_scaffold = outrows[-1].object_end

        # add a gap if we're not on the last subscaffold
        if i < len(superscaffold_rows) - 1:
            outrows.append(agp.GapRow(superscaffold_name,
                                      end_of_previous_scaffold + 1,
                                      end_of_previous_scaffold + gap_size,
                                      component_number_counter,
                                      length=gap_size,
                                      gap_type=gap_type,
                                      evidence=gap_evidence))
            component_number_counter += 1
            end_of_previous_scaffold = outrows[-1].object_end

    return outrows


def run(joins_list, outfile, agp, gap_size, gap_type, gap_evidence):
    """
    Runs the 'join' module of agptools.
    """
    # first, make a dict mapping the names of all scaffolds that will
    # be modified to an empty list which will later contain all agp rows
    # of that scaffold.
    scaffolds_to_join = {}
    for name in map(lambda s: s.lstrip('+-'), chain.from_iterable(joins_list)):
        scaffolds_to_join[name] = []

    # print all the rows to be output as-is and put the rows that need
    # to be modified first into the correct slot of the
    # scaffolds_to_join dict
    for row in agp:
        if row.object in scaffolds_to_join:
            scaffolds_to_join[row.object].append(row)
        else:
            print(row, file=outfile)

    # loop through each join group
    for join_group in joins_list:
        scaffold_rows = []
        # loop through the scaffolds in this join group
        for scaffold_name in join_group:
            if scaffold_name.startswith('-'):
                # reverse-complement if necessary
                scaffold_rows.append(flip.reverse_rows(
                    scaffolds_to_join[scaffold_name.lstrip('-')]))
            else:
                scaffold_rows.append(
                        scaffolds_to_join[scaffold_name.lstrip('+')])

        # print out all the rows for this join group
        for row in join_scaffolds(scaffold_rows, gap_size, gap_type):
            print(row, file=outfile)


