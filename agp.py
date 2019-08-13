"""
Functions for reading and writing AGP files.
"""

class AgpRow:
    """
    A single non-comment row of an AGP file. Because AGP is a weird
    kind of TSV where fields can have different meanings depending on
    the value of the component_type field, I'm representing this as a
    class with instance variables rather than an OrderedDict or
    something along those lines. See the NCBI documentation for more
    information:

    https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
    """
    def __init__(self, line):
        """
        Creates a new instance of AgpRow by parsing the given line of
        text.
        """
        splits = line.strip().split('\t')
        self.object = splits[0]
        self.object_beg = int(splits[1])
        self.object_end = int(splits[2])
        self.part_number = int(splits[3])
        self.component_type = splits[4]

        if self.component_type in ['N', 'U']:
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

    def __str__(self):
        """
        Returns the string representation of the AGP row as a line of
        text containing all the fields separated by tabs.
        """
        if self.is_gap:
            return '\t'.join(map(str, [self.object, self.object_beg,
                                       self.object_end, self.part_number,
                                       self.component_type, self.gap_length,
                                       self.gap_type, self.linkage,
                                       self.linkage_evidence]))
        else:
            return '\t'.join(map(str, [self.object, self.object_beg,
                                       self.object_end, self.part_number,
                                       self.component_type, self.component_id,
                                       self.component_beg, self.component_end,
                                       self.orientation]))

    def contains(self, position):
        """
        Returns true if position is within the bounds of this entry,
        false otherwise.

        Args:
            position (int): a genomic position in base pairs
        """
        return self.object_beg <= position and self.object_end >= position

class GapRow(AgpRow):
    def __init__(self, name, beginning, end, part_number, length=500,
                 gap_type='scaffold', linkage='yes', evidence='paired-end'):
        self.object = name
        self.object_beg, self.object_end = beginning, end
        self.part_number = part_number
        self.component_type, self.is_gap = 'N', True
        self.gap_length, self.gap_type = length, gap_type
        self.linkage, self.linkage_evidence = linkage, evidence


def read(infile):
    """
    Reads an AGP file, yielding rows as AgpRow instances and comment
    lines as plain strings.
    """
    for line in infile:
        if line.startswith('#'):
            yield line.strip()
        else:
            yield AgpRow(line)


def open_agp(filename):
    return read(open(filename))


