"""
Functions for assembling scaffolds from contigs based on an agp file
"""

from itertools import filterfalse
import textwrap
from typing import Iterable, IO

import screed

import agp


class NoSuchContigError(Exception):
    def __init__(self, contig_name):
        self.contig_name = contig_name

    def __str__(self):
        return 'FATAL: No contig named "{}" found.'.format(self.contig_name)


def complement(base: str):
    complement_dict = {
        "A": "T",
        "a": "t",
        "C": "G",
        "c": "g",
        "G": "C",
        "g": "c",
        "T": "A",
        "t": "a",
    }
    if base in complement_dict:
        return complement_dict[base]
    else:
        return base


def reverse_complement(sequence: str):
    return "".join(complement(s) for s in sequence[::-1])


def sequence_to_fasta(name: str, sequence: str, wrap: int = 60) -> str:
    """
    Given a sequence name and the sequence itself, make it into a fasta-
    formatted string.

    Arguments:
        name: the sequence id
        sequence: the sequence itself
        wrap: the number of bases per line of sequence

    Returns: the sequence in fasta format

    >>> sequence_to_fasta('chr1', 'ACTGGAGACATATAGTCCCACG', wrap=10)
    '>chr1\nACTGGAGACA\nTATAGTCCCA\nCG'
    """
    fasta_lines = []
    fasta_lines.append(f">{name}")
    fasta_lines += textwrap.wrap(sequence, width=wrap)
    return "\n".join(fasta_lines)


def run(
    contigs_fasta: screed.openscreed.Open, outfile: IO, agp_rows: Iterable[agp.AgpRow]
):
    """
    Given contigs in fasta format and their order and orientation into
    scaffolds in AGP format, outputs the assembled scaffolds in fasta
    format.

    Args:
        contigs_fasta: fasta iterator in screed format containing contigs
        outfile (file): file where scaffolds fasta should be written
        agp (iterable): iterable returning agp.AgpRow objects, each
            containing a single row of the agp file
    """
    # unfortunately, the contigs fasta file I'm writing this for has
    # variable line-length and is thus not faidx-able, so we have to
    # load it into memory :(
    contigs = {record.name: record.sequence for record in contigs_fasta}

    current_sequence = None
    current_chrom = None
    # loop through AGP, skipping comment lines, which my agp library
    # yields as strings
    for row in filterfalse(agp.is_string, agp_rows):
        # check if starting a new chromosome
        if row.object != current_chrom:
            # if this is not the first chromosome, output the previous
            # chromosome
            if current_chrom is not None:
                print(sequence_to_fasta(current_chrom, current_sequence), file=outfile)
            # start the new chromosome as an empty sequence
            current_chrom = row.object
            current_sequence = ""

        if row.is_gap:
            current_sequence += "N" * row.gap_length
        else:
            start, end = row.component_beg - 1, row.component_end
            if row.component_id not in contigs.keys():
                raise NoSuchContigError(row.component_id)
            component = contigs[row.component_id][start:end]
            if row.orientation == "-":
                component = reverse_complement(component)
            current_sequence += component

    print(sequence_to_fasta(current_chrom, current_sequence), file=outfile)
