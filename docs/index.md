# agptools

Tools for working with agp files

## Introduction
The [AGP][agp] format is a tab-separated table format describing how components
of a genome assembly fit together. NCBI accepts assemblies for submission in
the format of a fasta file giving the sequences of components (usually contigs)
along with an AGP file showing how these components are assembled into larger
pieces like scaffolds or chromosomes. For this reason, scaffolders such as
[SALSA][salsa] output this format alongside the fasta file of the scaffolded
assembly containing gaps and all.

Unfortunately, scaffolding never really works perfectly, so you invariably have
to correct mistakes or add in other sources of data such as synteny with a
related species or a phyiscal map to get an assembly to chromosome level. While
you could perform this manual curation process by editing the fasta file of
scaffolds, I think it is a lot easier to leave the fasta file of contigs intact
and move things around in the AGP file. For example, let's say the scaffolder
misorients a contig. Fixing this in the fasta would require taking a chunk from
the middle of a sequence, reverse orienting it, and pasting it back together.
Fixing it in the AGP file is as simple as finding the line corresponding to
that contig and changing the character in the orientation column from '+' to
'-'.

agptools is a suite of scripts for performing edits to an AGP file during this
manual curation stage of genome assembly.

## Installation
I'll get this on PyPI soon so that you can install this with one command
instead of three, but for now, you have to clone the repository first:
```bash
git clone https://github.com/esrice/agptools.git
cd agptools
pip install .
```

## Commands
agptools contains modules for the following operations:

* `split`: given one or more breakpoint coordinates, split a scaffold into two or
  more scaffolds
* `join`: given an ordered list of scaffolds and their desired orientations, join
  them into a single super-scaffold (or even a chromosome)
* `flip`: given two coordinates on a scaffold, reverse complement everything
  between these two coordinates
* `assemble`: given an AGP file and a fasta file containing the sequences of its
  components (e.g., contigs), output a fasta file containing the sequences of
  the components assembled into objects (e.g., scaffolds or chromosomes) as
  specified by the AGP file
* `rename`: rename scaffolds
* `remove`: remove scaffolds that you don't want in the final assembly
* `transform`: transform a bed file from contig coordinates to scaffold
  coordinates

Full documentation for each of them is available on the [commands](commands.md)
page.

[agp]: https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
[salsa]: https://github.com/marbl/SALSA

