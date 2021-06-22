"""
Functions for parsing bed files
"""


class EmptyRangeError(Exception):
    def __init__(self, bed_range):
        self.bed_range = bed_range

    def __str__(self):
        return "The range {}:{}-{} does not contain any contigs.".format(
            self.bed_range.chrom, self.bed_range.start, self.bed_range.end
        )


class BadRangeError(Exception):
    def __init__(self, bed_range):
        self.bed_range = bed_range

    def __str__(self):
        return (
            "The bed range {}:{}-{} starts and/or ends in the middle of a"
            "component. This is not supported in this module."
        ).format(self.bed_range.chrom, self.bed_range.start, self.bed_range.end,)


class BedRange:
    def __init__(self, chrom, start=None, end=None, extra_fields=[]):
        self.chrom = chrom
        self.start, self.end = start, end
        if len(extra_fields) >= 3:
            self.name, self.score, self.strand = extra_fields[:3]
            self.extra_fields = extra_fields[3:]
        else:
            self.name, self.score, self.strand = [None] * 3
            self.extra_fields = extra_fields

    def __str__(self):
        fields = [self.chrom]
        if self.start is not None and self.end is not None:
            fields += [self.start, self.end]
            if self.name and self.score and self.strand:
                fields += [self.name, self.score, self.strand]
            fields += self.extra_fields
        return "\t".join(map(str, fields))


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
        splits = line.strip().split("\t")
        if len(splits) == 1:
            yield BedRange(splits[0])
        else:
            yield BedRange(
                splits[0],
                start=int(splits[1]),
                end=int(splits[2]),
                extra_fields=splits[3:],
            )
