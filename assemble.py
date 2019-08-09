"""
Functions for assembling scaffolds from contigs based on an agp file
"""

import argparse
from itertools import filterfalse
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

import agp

def run(contigs_fasta, outfile, agp_rows):
    """
    Given contigs in fasta format and their order and orientation into
    scaffolds in AGP format, outputs the assembled scaffolds in fasta
    format.

    Args:
        contigs_fasta (iterable): iterable returning Bio.SeqRecord
            objects, each containing a single contig
        outfile (file): file where scaffolds fasta should be written
        agp (iterable): iterable returning agp.AgpRow objects, each
            containing a single row of the agp file
    """
    # unfortunately, the contigs fasta file I'm writing this for has
    # variable line-length and is thus not faidx-able, so we have to
    # load it into memory :(
    contigs = {record.id: record.seq for record in contigs_fasta}

    current_sequence = None
    current_chrom = None
    # loop through AGP, skipping comment lines, which my agp library
    # yields as strings
    for row in filter(lambda r: not isinstance(r, str), agp_rows):
        # check if starting a new chromosome
        if row.object != current_chrom:
            # if this is not the first chromosome, output the previous
            # chromosome
            if current_chrom is not None:
                record = SeqRecord(current_sequence, id=current_chrom,
                                   description='')
                print(record.format('fasta'), end='', file=outfile)
            # start the new chromosome as an empty sequence
            current_chrom = row.object
            current_sequence = Seq('', generic_dna)

        if row.is_gap:
            current_sequence += 'N' * row.gap_length
        else:
            start, end = row.component_beg - 1, row.component_end
            component = contigs[row.component_id][start:end]
            if row.orientation == '-':
                component = component.reverse_complement()
            current_sequence += component

    record = SeqRecord(current_sequence, id=current_chrom, description='')
    print(record.format('fasta'), end='', file=outfile)


# this is not meant to be run as a stand-alone executable but rather as
# part of the 'agptools' executable; nonetheless, who am I to stop the
# user from doing this?
if __name__ == '__main__':
    args = AssembleParser().parse_args()
    args.func(args)

