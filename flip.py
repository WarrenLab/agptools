#!/usr/bin/env python3
"""
Given an AGP file describing an assembly and a list of chromosomes to
flip, output a new AGP file with those chromosomes reoriented in
reverse-complement.
"""

import argparse
from functools import partial
from itertools import filterfalse
import sys

import agp
import bed

def reverse_rows(rows):
    pass # TODO implement


def flip(agp_rows, ranges_to_flip):
    # agp_rows is a generator, but we're going to have to store the
    # whole agp in memory
    agp_rows = list(agp_rows)

    for bed_range in ranges_to_flip:
        # list of tuples (i, agp_row) of rows in bed_range
        rows_to_reverse = []
        for i, row in enumerate(agp_rows):
            if not isinstance(row, str) and row.object == bed_range.chrom:
                # if no range in seq specified, flip the whole sequence
                if row.object_beg is None:
                    rows_to_reverse.append((i, row))
                # otherwise, flip only the part within range
                elif (row.object_beg >= bed_range.start
                      and row.object_end <= bed_range.end):
                    rows_to_reverse.append((i, row))

        # unzip and flip those rows we've collected
        indices, rows = map(list, zip(*rows_to_reverse))
        reversed_rows = reverse_rows(rows)

        # replace the original rows with the reversed rows
        for i, row in zip(indices, reversed_rows):
            agp_rows[i] = row

    return agp_rows


def run(segments_to_flip, outfile, agp_rows):
    for row in flip(agp_rows, segments_to_flip):
        print(row, file=outfile)


