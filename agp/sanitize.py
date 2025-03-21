"""
Functions for cleaning up an AGP to conform to NCBI rules.
"""

from typing import IO, Iterator, List, Sequence

from pyfaidx import Fasta

from agp import AgpRow
from agp.fasta import print_fasta


def divide_into_scaffolds(agp_in: Sequence[AgpRow]) -> Iterator[List[AgpRow]]:
    """
    Take an iterator over the rows of an AGP and turn it into an
    iterator of the scaffolds of an AGP, where each item is a list of
    rows with the same object ID.
    """
    current_scaffold = None
    current_scaffold_rows: List[AgpRow] = []
    for row in agp_in:
        if row.object != current_scaffold:
            if current_scaffold_rows:
                yield current_scaffold_rows
            current_scaffold = row.object
            current_scaffold_rows = []
        current_scaffold_rows.append(row)
    yield current_scaffold_rows


def run(agp_in: Sequence[AgpRow], agp_out: IO, contigs_in: Fasta, contigs_out: IO):
    contig_counter = 1
    for scaffold_rows in divide_into_scaffolds(agp_in):
        # NCBI does not like single-component scaffolds with "-"
        # orientation, so force this not to be the case
        if len(scaffold_rows) == 1:
            scaffold_rows[0].orientation = "+"
        for row in scaffold_rows:
            if not row.is_gap:
                # make a new contig containing exactly the range of this
                # row, and adjust the AGP accordingly
                component_part_sequence = contigs_in[row.component_id][
                    row.component_beg - 1 : row.component_end
                ].seq
                new_contig_name = f"contig_{contig_counter}"
                print_fasta(new_contig_name, component_part_sequence, contigs_out)
                contig_counter += 1
                row.component_id = new_contig_name
                row.component_beg = 1
                row.component_end = len(component_part_sequence)
            print(row, file=agp_out)
