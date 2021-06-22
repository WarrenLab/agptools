"""
Functions for transforming genomic coordinates from contig-based to
scaffold-based.
"""

import agp


class NoSuchContigError(Exception):
    def __init__(self, contig_name):
        self.contig_name = contig_name

    def __str__(self):
        return 'Cannot find contig called "{}" in AGP'.format(self.contig_name)


def run(bed_in, agp_in, bed_out):
    # this dict maps contig name to a tuple containing:
    # 0) scaffold name (str),
    # 1) start position of contig on scaffold (int),
    # 2) end position of contig on scaffold (int), and
    # 3) orientation of contig on scaffold (str)
    contig_to_scaffold = {}
    for row in (r for r in agp_in if not (agp.is_string(r) or r.is_gap)):
        contig_to_scaffold[row.component_id] = (
            row.object,
            row.object_beg,
            row.object_end,
            row.orientation,
        )

    for bed_row in bed_in:
        try:
            scaffold, start, end, orient = contig_to_scaffold[bed_row.chrom]
            bed_row.chrom = scaffold
            if orient == "-":
                tmp_end = end - bed_row.start + 1
                tmp_start = end - bed_row.end + 1
                bed_row.start, bed_row.end = tmp_start, tmp_end
                if bed_row.strand == "-":
                    bed_row.strand = "+"
                elif bed_row.strand == "+":
                    bed_row.strand = "-"
            else:
                bed_row.start += start - 1
                bed_row.end += start - 1
            print(bed_row, file=bed_out)
        except KeyError as ke:
            raise NoSuchContigError(bed_row.chrom) from ke
