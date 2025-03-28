#!/usr/bin/env python3
"""
Various tools for modifying and using agp files
"""
from __future__ import annotations

import argparse
import sys

import pyfaidx
import screed

import agp
from agp import assemble, bed, flip, join, remove, rename, sanitize, split, transform


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

    remove_parser = subparsers.add_parser(
        "remove", help="remove scaffolds from the assembly"
    )
    remove_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output AGP [STDOUT]",
        default=sys.stdout,
    )
    remove_parser.add_argument(
        "scaffolds_to_remove",
        type=remove.scaffolds_list_type,
        help="list of scaffolds to remove, one per line",
    )
    remove_parser.add_argument(
        "agp",
        nargs="?",
        type=agp.open_agp,
        help="AGP file to modify [STDIN]",
        default=agp.read(sys.stdin),
    )
    remove_parser.set_defaults(
        func=lambda a: remove.run(a.scaffolds_to_remove, a.outfile, a.agp)
    )

    rename_parser = subparsers.add_parser(
        "rename",
        help="rename scaffolds",
        description="Given a list of name changes, rename scaffolds accordingly",
    )
    rename_parser.add_argument(
        "-o",
        "--outfile",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output AGP [STDOUT]",
        default=sys.stdout,
    )
    rename_parser.add_argument(
        "renaming_file",
        type=rename.renaming_file_type,
        help="tsv containing columns current name, new name, and optional orientation",
    )
    rename_parser.add_argument(
        "agp",
        nargs="?",
        type=agp.open_agp,
        help="AGP file to modify [STDIN]",
        default=agp.read(sys.stdin),
    )
    rename_parser.set_defaults(
        func=lambda a: rename.run(a.renaming_file, a.outfile, a.agp)
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

    # --- 'transform' command options ---
    sanitize_parser = subparsers.add_parser(
        "sanitize",
        help="clean up agp and contigs so NCBI will accept them",
        description="""AGP format allows using only part of a contig in a scaffold,
            but NCBI does not allow this. This module splits contigs into parts
            according to the AGP, and makes a new agp and fasta that NCBI will
            not complain about.""",
    )
    sanitize_parser.add_argument(
        "contigs_fasta_in",
        type=pyfaidx.Fasta,
        help="input contigs corresponding to input AGP",
    )
    sanitize_parser.add_argument(
        "contigs_fasta_out",
        type=argparse.FileType("w"),
        help="path where sanitized contigs should be written",
    )
    sanitize_parser.add_argument(
        "agp_in",
        nargs="?",
        type=agp.open_agp,
        help="AGP file to modify [STDIN]",
        default=agp.read(sys.stdin),
    )
    sanitize_parser.add_argument(
        "-o",
        "--agp-out",
        type=argparse.FileType("w"),
        nargs="?",
        help="where to write output AGP [STDOUT]",
        default=sys.stdout,
    )
    sanitize_parser.set_defaults(
        func=lambda a: sanitize.run(
            a.agp_in, a.agp_out, a.contigs_fasta_in, a.contigs_fasta_out
        )
    )

    return parser.parse_args()


def main():
    """main method of script"""
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()  # pragma: no cover
