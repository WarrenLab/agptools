"""
Functions for working with fasta format.
"""

from typing import IO


def print_fasta(name: str, sequence: str, outfile: IO, wrap: int = 60):
    """
    Given a sequence name and the sequence itself, print a fasta-
    formatted string of it.

    Arguments:
        name: the sequence id
        sequence: the sequence itself
        outfile: where to write the sequence to
        wrap: the number of bases per line of sequence

    Examples:
        >>> import sys
        >>> print_fasta(
        ...    'chr1', 'ATCGACTGATCGACTGACTGACTACTG', outfile=sys.stdout, wrap=10
        ... )
        >chr1
        ATCGACTGAT
        CGACTGACTG
        ACTACTG
    """
    print(f">{name}", file=outfile)
    for start_pos in range(0, len(sequence), wrap):
        print(sequence[start_pos : start_pos + wrap], file=outfile)
