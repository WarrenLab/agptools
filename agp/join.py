"""
Functions for joining scaffolds
"""

import re
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import chain
from typing import Dict, Iterable, Iterator, List, Match, Optional, TextIO, Union, cast

from agp import AgpRow, GapRow
from agp.flip import reverse_rows

scaffold_regex = re.compile(r"([a-zA-Z0-9]+)_([a-zA-Z0-9.]+)")
sequence_name_regex = re.compile(r"^[a-zA-Z0-9.]+$")


class ScaffoldNotFoundError(Exception):
    """Scaffold in joins list not found in agp"""

    pass


class ScaffoldUsedTwiceError(Exception):
    """Scaffold(s) used more than once in joins file"""

    pass


class BadSequenceNameError(Exception):
    """Sequence name is not allowed"""

    def __init__(self, name: str):
        self.name = name

    def __str__(self) -> str:
        return f"Bad scaffold name: '{self.name}'. Only [a-zA-Z0-9._] allowed."


@dataclass
class JoinGroup(Sequence):
    """A group of scaffolds to be joined together

    This class represents a group of scaffolds that should be joined
    together as specified in a joins file. Optionally, a group can be given
    a name to use for the new superscaffold instead of using
    `make_superscaffold_name()` to make one up based on the names of the
    components. If accessed as a `Sequence`, it provides direct access to
    the underlying list of scaffolds, e.g.,

    >>> from agp.join import JoinGroup
    >>> join_group = JoinGroup(["scaffold_1", "-scaffold_2", "+scaffold_3"])
    >>> len(join_group)
    3
    >>> join_group[1]
    '-scaffold_2'

    Attributes:
        scaffolds: a list of scaffolds to join together, in the correct
            order. A scaffold name can be preceded by "+" or "-" to
            indicate orientation; if none is specified, "+" is assumed.
        name: a name to use for the new superscaffold in the output agp
    """

    scaffolds: List[str]
    name: Optional[str] = None

    def __getitem__(self, index) -> str:
        return self.scaffolds[index]

    def __len__(self) -> int:
        return len(self.scaffolds)


def joins_type(filename: str) -> List[JoinGroup]:
    """
    argparse type function for file listing scaffold joins

    Args:
        filename: filename giving path to a file listing
            scaffolds to join. Each line should contain a comma-
            separated list of scaffolds to be joined into a single
            scaffold. Optionally, the file can contain a second column
            after a tab that gives each group a name.

    Returns:
        list of join groups. Scaffold names may have "+" or "-"
            before their name; that is the job of downstream code to handle.

    Raises:
        ScaffoldUsedTwiceError: if a scaffold is used multiple times
    """
    joins: List[JoinGroup] = []
    with open(filename) as joins_file:
        for line in joins_file:
            columns = line.strip().split("\t")
            join_group = JoinGroup(columns[0].split(","))
            if len(columns) > 1:
                if not sequence_name_regex.match(columns[1]):
                    raise BadSequenceNameError(columns[1])
                join_group.name = columns[1]
            joins.append(join_group)

    # look for scaffolds used more than once
    scaffold_counts = Counter(s.lstrip("+-") for s in chain.from_iterable(joins))
    reused_scaffolds = [scaf for scaf, count in scaffold_counts.items() if count > 1]
    if reused_scaffolds:
        raise ScaffoldUsedTwiceError(f"Scaffolds used 2+ times: {reused_scaffolds}")

    return joins


def make_superscaffold_name(subscaffold_names: Iterable[str]) -> str:
    """
    Comes up with a nice superscaffold name based on the given
    subscaffold names. If the subscaffold names are all in the format
    '[prefix]_[suffix]', where prefix is the same for all subscaffolds,
    then the superscaffold name is '[prefix]_[suffix1]p[suffix2]p[etc.]'
    Otherwise, it's the names of all scaffolds concatenated with 'p'.

    Args:
        subscaffold_names (list(str)): list of names of subscaffolds
            being combined into a superscaffold

    Returns:
        superscaffold_name (str): a name for the new scaffolds based
            on the subscaffold names
    """
    # TODO this casting is necessary but ugly as sin. Find a nicer way
    matches = list(map(scaffold_regex.search, subscaffold_names))
    if all(matches):
        if all(
            cast(Match, m).group(1) == cast(Match, matches[0]).group(1) for m in matches
        ):
            prefix = cast(Match, matches[0]).group(1)
            suffix = "p".join([cast(Match, m).group(2) for m in matches])
            return "{}_{}".format(prefix, suffix)
    return "p".join(subscaffold_names)


def join_scaffolds(
    superscaffold_rows: List[List[AgpRow]],
    gap_size: int = 500,
    gap_type: str = "scaffold",
    gap_evidence: str = "na",
    name: str = None,
) -> List[AgpRow]:
    """
    Transforms agp rows from multiple scaffolds into agp rows for a
    single superscaffold containing all the given scaffolds in the
    specified order. Scaffolds should be oriented properly before
    using this function.

    Args:
        superscaffold_rows: a list of lists of agp rows. Each sub-list is a
            list of all rows in a single scaffold. These scaffolds will be
            joined in the order given
        gap_size: length of the new gaps created by joining scaffolds
            together
        gap_type: what to call the new gaps in the agp file
        gap_evidence: evidence type for linkage across gap
        name: name of new superscaffold. If None, this function will come
            up with a name based on the names of subscaffolds contained by
            this superscaffold

    Returns:
        a list of AgpRow instances containing the new scaffold specification
    """
    # make a nice name for the new superscaffold we are creating
    subscaffold_names = (s[0].object for s in superscaffold_rows)
    if name is None:
        superscaffold_name = make_superscaffold_name(subscaffold_names)
    else:
        superscaffold_name = name

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
            outrows.append(
                GapRow(
                    superscaffold_name,
                    end_of_previous_scaffold + 1,
                    end_of_previous_scaffold + gap_size,
                    component_number_counter,
                    length=gap_size,
                    gap_type=gap_type,
                    evidence=gap_evidence,
                )
            )
            component_number_counter += 1
            end_of_previous_scaffold = outrows[-1].object_end

    return outrows


def run(
    joins_list: List[JoinGroup],
    outfile: TextIO,
    agp_infile: Iterator[Union[str, AgpRow]],
    gap_size: int = 500,
    gap_type: str = "scaffold",
    gap_evidence: str = "na",
):
    """
    Runs the 'join' module of agptools.
    """
    # first, make a dict mapping the names of all scaffolds that will
    # be modified to an empty list which will later contain all agp rows
    # of that scaffold.
    scaffolds_to_join: Dict[str, List[AgpRow]] = {}
    for name in (s.lstrip("+-") for s in chain.from_iterable(joins_list)):
        scaffolds_to_join[name] = []

    # print all the rows to be output as-is and put the rows that need
    # to be modified first into the correct slot of the
    # scaffolds_to_join dict
    for row in agp_infile:
        if isinstance(row, AgpRow) and row.object in scaffolds_to_join:
            scaffolds_to_join[row.object].append(row)
        else:
            print(row, file=outfile)

    # make sure we found agp entries for all the scaffolds
    for scaffold_name, rows in scaffolds_to_join.items():
        if len(rows) == 0:
            raise ScaffoldNotFoundError(f"Scaffold {scaffold_name} not found in agp.")

    # loop through each join group
    for i, join_group in enumerate(joins_list):
        scaffold_rows = []
        # loop through the scaffolds in this join group
        for scaffold_name in join_group:
            if scaffold_name.startswith("-"):
                # reverse-complement if necessary
                scaffold_rows.append(
                    reverse_rows(scaffolds_to_join[scaffold_name.lstrip("-")])
                )
            else:
                scaffold_rows.append(scaffolds_to_join[scaffold_name.lstrip("+")])

        # print out all the rows for this join group
        for row in join_scaffolds(
            scaffold_rows, gap_size, gap_type, name=join_group.name
        ):
            print(row, file=outfile)
