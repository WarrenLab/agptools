"""
Functions for parsing bed files
"""
from dataclasses import dataclass
from typing import Iterator, List, Optional, TextIO, Union


class EmptyRangeError(Exception):
    """Error raised when bed range contains no contigs"""

    def __init__(self, bed_range: "BedRange"):
        self.bed_range = bed_range

    def __str__(self):
        return "The range {}:{}-{} does not contain any contigs.".format(
            self.bed_range.chrom, self.bed_range.start, self.bed_range.end
        )


class BadRangeError(Exception):
    """Error raised when range starts in middle of component"""

    def __init__(self, bed_range: "BedRange"):
        self.bed_range = bed_range

    def __str__(self):
        return (
            "The bed range {}:{}-{} starts and/or ends in the middle of a"
            "component. This is not supported in this module."
        ).format(
            self.bed_range.chrom,
            self.bed_range.start,
            self.bed_range.end,
        )


class ParsingError(Exception):
    """Error parsing bed file"""

    pass


@dataclass
class BedRange:
    """Single entry (i.e., line) of a bed file."""

    chrom: str
    """Name of the chromosome or sequence of the bed range"""
    start: Optional[int] = None
    """Start position of the bed range"""
    end: Optional[int] = None
    """End position of the bed range"""
    strand: Optional[str] = None
    """Strand of the bed range (either '+' or '-')"""

    def __str__(self) -> str:
        fields: List[Union[str, int]] = [self.chrom]
        if self.start and self.end:
            fields += [self.start, self.end]
            if self.strand:
                fields += [self.strand]
        return "\t".join(map(str, fields))


def open_bed(filename: str) -> Iterator[BedRange]:
    """Open and parse a bed file

    Args:
        filename: path to a bed file
    """
    return read(open(filename))


def read(bedfile: TextIO) -> Iterator[BedRange]:
    """Read and parse a bed file

    Read a bed file line by line. If a line contains start and end
    positions, yield the sequence name and start and end positions as
    a BedRange; if the line contains only a sequence name, yield a
    BedRange with {start,end}=None.

    Args:
        bedfile: a bed-formatted file

    Returns:
        an iterator that yields single lines of the file at a time
    """
    for i, line in enumerate(bedfile):
        splits = line.strip().split("\t")
        try:
            if len(splits) == 1:
                yield BedRange(splits[0])
            elif len(splits) == 2:
                raise ParsingError(f"Line {i+1} of bed misformatted.")
            elif len(splits) == 3:
                yield BedRange(
                    splits[0],
                    start=int(splits[1]),
                    end=int(splits[2]),
                )
            else:
                yield BedRange(
                    splits[0],
                    start=int(splits[1]),
                    end=int(splits[2]),
                    strand=splits[3],
                )
        except (ValueError, IndexError):
            raise ParsingError(f"Line {i+1} of bed misformatted.")
