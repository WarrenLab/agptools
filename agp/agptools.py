#!/usr/bin/env python3
"""
Various tools for modifying and using agp files
"""

import argparse
from functools import partial
import sys

import screed

import agp
from agp import bed, assemble, compose, flip, join, split, transform


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(
        required=True, dest="command", help="command to run"
    )

    # --- 'split' command options ---
    split_parser = subparsers.add_parser(
        "split",
        help="split scaffolds into subscaffolds",
        description="""Split scaffolds into subscaffolds by breaking
                           at gaps.""",
    )
    split_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output AGP [STDOUT]",
        default=sys.stdout,
    )
    split_parser.add_argument(
        "breakpoints",
        type=split.breakpoints_type,
        help="""File listing all places where scaffolds
                              should be broken. First column is name of
                              scaffold; second column is a comma-separated
                              list of locations where scaffolds should be
                              broken. These positions must be within a gap.""",
    )
    split_parser.add_argument(
        "agp",
        nargs="?",
        type=agp.open_agp,
        help="AGP file to modify [STDIN]",
        default=agp.read(sys.stdin),
    )
    split_parser.set_defaults(func=lambda a: split.run(a.breakpoints, a.outfile, a.agp))

    # --- 'join' command options ---
    join_parser = subparsers.add_parser(
        "join",
        help="join multiple scaffolds into a superscaffold",
        description="Join groups of scaffolds in an AGP into single " "scaffolds",
    )
    join_parser.add_argument(
        "-n",
        "--gap-size",
        type=int,
        default=500,
        help="size of new gaps to be created [500]",
    )
    join_parser.add_argument(
        "-t",
        "--gap-type",
        default="scaffold",
        help="what to call new gaps in AGP [scaffold]",
    )
    join_parser.add_argument(
        "-e", "--gap-evidence", default="na", help="evidence type for gaps [na]"
    )
    join_parser.add_argument(
        "-r",
        "--renaming-file",
        type=join.renaming_type,
        help="file containing names of new objects to "
        "use, instead of making up names, one per line",
    )
    join_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output AGP [STDOUT]",
        default=sys.stdout,
    )
    join_parser.add_argument(
        "joins_list",
        type=join.joins_type,
        help="file listing scaffolds "
        "to join. Each line should contain a comma-separated list of "
        "scaffolds to be joined into a single scaffold. Scaffold names "
        "may be prefixed by + or - to indicate orientation; otherwise + "
        "is assumed.",
    )
    join_parser.add_argument(
        "agp",
        nargs="?",
        type=agp.open_agp,
        help="AGP file to modify [STDIN]",
        default=agp.read(sys.stdin),
    )
    join_parser.set_defaults(
        func=lambda a: join.run(
            a.joins_list,
            a.outfile,
            a.agp,
            a.gap_size,
            a.gap_type,
            a.gap_evidence,
            names=a.renaming_file,
        )
    )

    # --- 'flip' command options ---
    flip_parser = subparsers.add_parser(
        "flip",
        help="reverse complement scaffolds or components thereof",
        description="""
                Reverse complement scaffolds or components thereof in an agp
                file. Specify ranges to flip in scaffold coordinates.
                Breakpoints inside of contigs not yet supported.""",
    )
    flip_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output AGP [STDOUT]",
        default=sys.stdout,
    )
    flip_parser.add_argument(
        "segments_to_flip",
        type=bed.open_bed,
        help="""list of segments to flip, in bed format.
                             To flip a whole sequence, just put the name of
                             the sequence on a line by itself.""",
    )
    flip_parser.add_argument(
        "agp",
        nargs="?",
        type=agp.open_agp,
        help="AGP file to modify [STDIN]",
        default=agp.read(sys.stdin),
    )
    flip_parser.set_defaults(
        func=lambda a: flip.run(a.segments_to_flip, a.outfile, a.agp)
    )

    # --- 'assemble' command options ---
    assemble_parser = subparsers.add_parser(
        "assemble",
        help="output scaffolds in fasta " "format based on agp and contigs fasta",
    )
    assemble_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        help="where to write fasta of scaffolds [STDOUT]",
        default=sys.stdout,
    )
    assemble_parser.add_argument(
        "contigs_fasta", type=screed.open, help="Assembly to flip scaffolds in"
    )
    assemble_parser.add_argument(
        "agp",
        type=agp.open_agp,
        default=agp.read(sys.stdin),
        help="AGP file assembling contigs into scaffolds [STDIN]",
    )
    assemble_parser.set_defaults(
        func=lambda a: assemble.run(a.contigs_fasta, a.outfile, a.agp)
    )

    # --- 'transform' command options ---
    transform_parser = subparsers.add_parser(
        "transform",
        help="transform coordinates in bed to scaffold-based",
        description="""Given a bed file with coordinates on contigs and an
            agp file assembling contigs into scaffolds, transform the bed file
            into scaffold coordinates.""",
    )
    transform_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output bed [STDOUT]",
        default=sys.stdout,
    )
    transform_parser.add_argument(
        "bed", type=bed.open_bed, help="bed file to transform"
    )
    transform_parser.add_argument(
        "agp",
        nargs="?",
        type=agp.open_agp,
        help="AGP file to modify [STDIN]",
        default=agp.read(sys.stdin),
    )
    transform_parser.set_defaults(func=lambda a: transform.run(a.bed, a.agp, a.outfile))

    # --- 'compose' command options ---
    compose_parser = subparsers.add_parser(
        "compose",
        help="turn two agps assembling a genome to "
        "different levels into a single agp.",
        description="""Given two agp
            files, each representing a different layer of an assembly, compose
            them into a single agp file assembling the lowest level into the
            highest level. For example, if you have one AGP assembling contigs
            into scaffolds, and another AGP assembling scaffolds into
            chromosomes, you could use this to make a new AGP assembling
            contigs directly into chromosomes.""",
    )
    compose_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output bed [STDOUT]",
        default=sys.stdout,
    )
    compose_parser.add_argument(
        "-p",
        "--print-unused",
        action="store_true",
        help="""print all objects from inner AGP that
                                were not used in outer AGP""",
    )
    compose_parser.add_argument(
        "outer_agp",
        type=agp.open_agp,
        help="""AGP file of the outer assembly step
                                (e.g., scaffolds into chromosomes)""",
    )
    compose_parser.add_argument(
        "inner_agp",
        type=agp.open_agp,
        help="""AGP file of the inner assembly step
                                (e.g., contigs into scaffolds)""",
    )
    compose_parser.set_defaults(
        func=lambda a: compose.run(a.outer_agp, a.inner_agp, a.outfile, a.print_unused)
    )

    return parser.parse_args()


def main():
    """main method of script"""
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()