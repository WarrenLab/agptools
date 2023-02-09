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


def run(output_prefix: str, contigs_in: Fasta, agp_in: Sequence[AgpRow]):
    # open up the necessary output files
    chromosomes_agp_out = open(f"{output_prefix}.chromosomes.agp", "w")
    unplaced_agp_out = open(f"{output_prefix}.unplaced.agp", "w")
    contigs_out = open(f"{output_prefix}.contigs.fa", "w")

    contig_counter = 1
    for scaffold_rows in divide_into_scaffolds(agp_in):
        # check whether scaffold is placed or not so that we can output it to
        # the correct file, and also make sure it's oriented correctly if it is
        # an unplaced single-component scaffold
        if scaffold_rows[0].object.startswith("chr"):  # placed scaffold
            agp_out = chromosomes_agp_out
        else:  # unplaced scaffold
            agp_out = unplaced_agp_out
            # NCBI does not like unplaced single-component scaffolds with "-"
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
