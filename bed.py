"""
Functions for parsing bed files
"""

import collections

BedRange = collections.namedtuple('BedRange', ['chrom', 'start', 'end'])

def open_bed(filename):
    return read(open(filename))

def read(bedfile):
    """
    Read a bed file line by line. If a line contains start and end
    positions, yield the sequence name and start and end positions as
    a BedRange; if the line contains only a sequence name, yield a
    BedRange with {start,end}=None.
    """
    for line in bedfile:
        splits = line.strip().split('\t')
        if len(splits) == 1:
            yield BedRange(splits[0], None, None)
        else:
            yield BedRange(splits[0], int(splits[1]), int(splits[2]))


