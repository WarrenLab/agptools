#!/usr/bin/env python3
"""
Given a fasta of contigs and a list of breakpoints, break the contigs
at these breakpoints.
"""
import argparse
import sys
from functools import partial

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def break_contig(contig, breakpoints):
    """
    Given a contig and a list of breakpoints, break the contig at these
    breakpoints and yield the pieces as SeqRecords with ids
    'contigID_1', 'contigID_2', etc.

    Args:
        contig (SeqRecord): a contig to be broken into pieces
        breakpoints (list(int)): a list of breakpoints where contig
            should be broken. Position of breakpoint is included in
            the contig before the breakpoint, not the contig after
            the breakpoint. E.g., a contig that is 100 bases long
            broken at position 53 will result in one contig from
            1-53 and another from 54-100 (genomic coordinates, not
            array coordinates)

    Returns:
        a generator of len(breakpoints)+1 SeqRecord instances, one for
        each broken piece of the contig, with id as the original contig
        id with '_{i}' suffixed to it, where i is the index of that
        piece, starting with 1.
    """
    previous_breakpoint = 0
    for i, breakpoint in enumerate(sorted(breakpoints)):
        contig_part_seq = contig.seq[previous_breakpoint:breakpoint]
        contig_part_id = "{}_{}".format(contig.id, i + 1)
        yield SeqRecord(contig_part_seq, id=contig_part_id, description="")
        previous_breakpoint = breakpoint
    # get last part of contig
    contig_part_seq = contig.seq[previous_breakpoint:]
    contig_part_id = "{}_{}".format(contig.id, len(breakpoints) + 1)
    yield SeqRecord(contig_part_seq, id=contig_part_id, description="")


def breakpoints_type(filename):
    """
    argparse type function for breakpoints file. Returns breakpoints
    as a dictionary with key as contig name and value as a list of
    breakpoint positions on that contig
    """
    breakpoints = {}
    with open(filename) as breakpoints_file:
        for line in breakpoints_file:
            contig_name, positions = line.strip().split()
            breakpoints[contig_name] = list(map(int, positions.split(",")))
    return breakpoints


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        default=sys.stdout,
        help="where to write broken contigs [STDOUT]",
    )
    parser.add_argument(
        "contigs",
        type=partial(SeqIO.parse, format="fasta"),
        help="fasta file containing contig sequences",
    )
    parser.add_argument(
        "breakpoints",
        type=breakpoints_type,
        help="list of breakpoints where each line is in the "
        'format "[contig_name] [comma-separated list of '
        'positions]"',
    )
    return parser.parse_args()


def main():
    args = parse_args()
    for contig in args.contigs:
        if contig.id in args.breakpoints:
            breakpoints = args.breakpoints[contig.id]
            for contig_piece in break_contig(contig, breakpoints):
                print(contig_piece.format("fasta"), end="", file=args.outfile)
        else:
            print(contig.format("fasta"), end="", file=args.outfile)


if __name__ == "__main__":
    main()
