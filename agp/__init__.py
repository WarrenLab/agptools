"""Functions for reading and writing AGP files."""
from typing import Iterator, TextIO, Union


class AgpFormatError(Exception):
    """Invalidly formatted AGP"""

    def __init__(self, line):
        self.line = line

    def __str__(self) -> str:
        return 'Invalid AGP: "{}"'.format(self.line)


class AgpRow:
    """A single non-comment row of an AGP

    A single non-comment row of an AGP file. Because AGP is a weird
    kind of TSV where fields can have different meanings depending on
    the value of the component_type field, I'm representing this as a
    class with instance variables rather than an OrderedDict or
    something along those lines. See the NCBI documentation for more
    information:

    https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
    """

    def __init__(self, line: str):
        """Creates a new instance of AgpRow

        Creates a new instance of AgpRow by parsing the given line of
        text.

        Args:
            line: a line of text in AGP format
        """
        try:
            splits = line.strip().split("\t")
            self.object = splits[0]
            self.object_beg = int(splits[1])
            self.object_end = int(splits[2])
            self.part_number = int(splits[3])
            self.component_type = splits[4]

            if self.component_type in ["N", "U"]:
                self.is_gap = True
                self.gap_length = int(splits[5])
                self.gap_type = splits[6]
                self.linkage = splits[7]
                self.linkage_evidence = splits[8]
            else:
                self.is_gap = False
                self.component_id = splits[5]
                self.component_beg = int(splits[6])
                self.component_end = int(splits[7])
                self.orientation = splits[8]
        except (ValueError, IndexError) as e:
            raise AgpFormatError(line) from e

    def __str__(self) -> str:
        """
        Returns the string representation of the AGP row as a line of
        text containing all the fields separated by tabs.
        """
        if self.is_gap:
            return "\t".join(
                map(
                    str,
                    [
                        self.object,
                        self.object_beg,
                        self.object_end,
                        self.part_number,
                        self.component_type,
                        self.gap_length,
                        self.gap_type,
                        self.linkage,
                        self.linkage_evidence,
                    ],
                )
            )
        else:
            return "\t".join(
                map(
                    str,
                    [
                        self.object,
                        self.object_beg,
                        self.object_end,
                        self.part_number,
                        self.component_type,
                        self.component_id,
                        self.component_beg,
                        self.component_end,
                        self.orientation,
                    ],
                )
            )

    def __eq__(self, other):
        if self.is_gap != other.is_gap:
            return False

        common_field_equality = (
            self.object == other.object
            and self.object_beg == other.object_beg
            and self.object_end == other.object_end
            and self.part_number == other.part_number
            and self.component_type == other.component_type
        )

        if self.is_gap:
            return (
                common_field_equality
                and self.gap_length == other.gap_length
                and self.gap_type == other.gap_type
                and self.linkage == other.linkage
                and self.linkage_evidence == other.linkage_evidence
            )
        else:
            return (
                common_field_equality
                and self.component_id == other.component_id
                and self.component_beg == other.component_beg
                and self.component_end == other.component_end
                and self.orientation == other.orientation
            )

    def contains(self, position: int) -> bool:
        """Check whether a row contains a position

        Determine whether this row contains a given coordinate.

        Args:
            position: a genomic position in base pairs

        Returns:
            true if position is within the bounds of this entry, false
            otherwise
        """
        return self.object_beg <= position and self.object_end >= position


class GapRow(AgpRow):
    def __init__(
        self,
        name: str,
        beginning: int,
        end: int,
        part_number: int,
        length: int = 500,
        gap_type: str = "scaffold",
        linkage: str = "yes",
        evidence: str = "paired-end",
    ):
        self.object = name
        self.object_beg, self.object_end = beginning, end
        self.part_number = part_number
        self.component_type, self.is_gap = "N", True
        self.gap_length, self.gap_type = length, gap_type
        self.linkage, self.linkage_evidence = linkage, evidence


def is_string(var) -> bool:
    """Determine whether var is a string

    Determine whether `var` is a string.

    Args:
        var: variable to check for stringiness

    Returns:
        true if `var` is a string, false otherwise
    """
    return isinstance(var, str)


def read(infile: TextIO) -> Iterator[Union[str, AgpRow]]:
    """Read an AGP file

    Read an AGP file, yielding rows as AgpRow instances and comment
    lines as plain strings.

    Args:
        infile: input file in AGP format

    Yields:
        each row as either a string containing the entire line for
        comments, or an AgpRow for non-comment lines
    """
    for line in infile:
        if line.startswith("#"):
            yield line.strip()
        else:
            yield AgpRow(line)


def open_agp(filename: str) -> Iterator[Union[str, AgpRow]]:
    """Open and read an AGP file

    Open an AGP file from a path, and then parse it

    Args:
        filename: path to an AGP file to read

    Returns:
        an Iterator over the rows of the file
    """
    return read(open(filename))
