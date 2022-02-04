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
manual curation stage of genome assembly. It contains modules for the following
operations:
* split: given one or more breakpoint coordinates, split a scaffold into two or
  more scaffolds
* join: given an ordered list of scaffolds and their desired orientations, join
  them into a single super-scaffold (or even a chromosome)
* flip: given two coordinates on a scaffold, reverse complement everything
  between these two coordinates
* assemble: given an AGP file and a fasta file containing the sequences of its
  components (e.g., contigs), output a fasta file containing the sequences of
  the components assembled into objects (e.g., scaffolds or chromosomes) as
  specified by the AGP file.

## Installation
I'll get this on PyPI soon so that you can install this with one command
instead of three, but for now, you have to clone the repository first:
```bash
git clone https://github.com/esrice/agptools.git
cd agptools
pip install .
```

## Example AGP
Here is an example input AGP file that will be used to demonstrate the usage
of the modules in the following sections:
```
scaffold_16  1         1096465   1   W  tig00005080 1    1096465 -
scaffold_16  1096466   1096965   2   N  500     scaffold        yes  na
scaffold_16  1096966   1973201   3   W  tig00001012 1    876236  +
scaffold_16  1973202   1973701   4   N  500     scaffold        yes  na
scaffold_16  1973702   4258994   5   W  tig00182876 1    2285293 -
scaffold_16  4258995   4259494   6   N  500     scaffold        yes  na
scaffold_16  4259495   11764263  7   W  tig00000113 1    7504769 +
scaffold_16  11764264  11764763  8   N  500     scaffold        yes  na
scaffold_16  11764764  13768005  9   W  tig00004962 1    2003242 -
scaffold_16  13768006  13768505  10  N  500     scaffold        yes  na
scaffold_16  13768506  17994060  11  W  tig00004933 1    4225555 -
scaffold_16  17994061  17994560  12  N  500     scaffold        yes  na
scaffold_16  17994561  21066363  13  W  tig00000080 1    3071803 -
scaffold_16  21066364  21066863  14  N  500     scaffold        yes  na
scaffold_16  21066864  21807834  15  W  tig00183148 1    740971  +
```
It contains eight contigs (and seven gaps, naturally).

## flip
The flip command reverse complements a segment of a scaffold. The use case for
this is the scaffolder putting piece of a scaffold in the wrong orientation.

The required arguments for the flip command are a list of flips to make and the
input agp you want to modify. The list of flips has three columns:
1. The name of the scaffold you want to flip a piece of
2. The beginning position in base pairs, in scaffold coordinates, of the piece
   you want to flip.
3. The end position in base pairs, in scaffold coordinates, of the piece you
   want to flip.

Both coordinates are in the usual DNA range format (i.e., "1-100" takes
everything from the first base of the sequence to the 100th base of the
sequence, including the first and 100th base). The begin coordinate must be
at the beginning of a component and the end coordinate must be at the end of
a component. To do something more complicated, like reverse complement only
part of a contig, use the split and join modules instead.

Here is an example flips file:
```
scaffold_16   1        1973201
scaffold_16   11764764 21066363
```

Here is the command you would use to perform these flips on the example file
shown in the introduction:
```bash
agptools flip flips.tsv test_in.agp > test_flip.agp
```

And here is the output:
```
scaffold_16  1         876236    1    W   tig00001012 1    876236      -
scaffold_16  876237    876736    2    N   500     scaffold        yes     na
scaffold_16  876737    1973201   3    W   tig00005080 1    1096465     +
scaffold_16  1973202   1973701   4    N   500     scaffold        yes     na
scaffold_16  1973702   4258994   5    W   tig00182876 1    2285293     -
scaffold_16  4258995   4259494   6    N   500     scaffold        yes     na
scaffold_16  4259495   11764263  7    W   tig00000113 1    7504769     +
scaffold_16  11764264  11764763  8    N   500     scaffold        yes     na
scaffold_16  11764764  14836566  9    W   tig00000080 1    3071803     +
scaffold_16  14836567  14837066  10   N   500     scaffold        yes     na
scaffold_16  14837067  19062621  11   W   tig00004933 1    4225555     +
scaffold_16  19062622  19063121  12   N   500     scaffold        yes     na
scaffold_16  19063122  21066363  13   W   tig00004962 1    2003242     +
scaffold_16  21066364  21066863  14   N   500     scaffold        yes     na
scaffold_16  21066864  21807834  15   W   tig00183148 1    740971      +
```

## split
The split module breaks a scaffold into multiple scaffolds. There are two
common cases where this operation may be necessary:
1. The scaffolder (or even contig assembler) joined two sequences together that
   aren't actually on the same chromosome.
2. The scaffolder joined two sequences together that are in fact on the same
   chromosome, but it did it in the wrong order, so you want to split the
   scaffold up into pieces and then put the pieces back together in a different
   order later with the `join` module.

The two required inputs to the split module are the AGP you want to operate on
and a list of splits you want to make. The output is the edited AGP file. The
format of the list of splits is a tab-separated file with the following
columns:
1. Name of the scaffold you want to split into parts
2. A comma-separated list of breakpoint coordinates where you want to break the
   scaffold

Let's say you want to split the example AGP into three pieces: one containing
the first three contigs, one containing the next four, and one containing the
last contig. Here is what the relevant line of your splits file would look
like:
```
scaffold_16	 4258995,21066364
```
Note that I used the BEGIN coordinate of the gaps, but you can use any
coordinate inside the gap and it will have the same result.

Here is the command:
```
agptools split splits.txt test.agp > split_out.agp
```
Here is the output:
```
scaffold_16.1	1	     1096465	1	W	tig00005080	 1	1096465	-
scaffold_16.1	1096466	 1096965	2	N	500	scaffold	yes	na
scaffold_16.1	1096966	 1973201	3	W	tig00001012	 1	876236	+
scaffold_16.1	1973202	 1973701	4	N	500	scaffold	yes	na
scaffold_16.1	1973702	 4258994	5	W	tig00182876	 1	2285293	-
scaffold_16.2	1	     7504769	1	W	tig00000113	 1	7504769	+
scaffold_16.2	7504770	 7505269	2	N	500	scaffold	yes	na
scaffold_16.2	7505270	 9508511	3	W	tig00004962	 1	2003242	-
scaffold_16.2	9508512	 9509011	4	N	500	scaffold	yes	na
scaffold_16.2	9509012	 13734566	5	W	tig00004933	 1	4225555	-
scaffold_16.2	1373456	 13735066	6	N	500	scaffold	yes	na
scaffold_16.2	13735067 16806869	7	W	tig00000080	 1	3071803 -
scaffold_16.3	1	     740971     1	W	tig00183148	 1	740971	+
```

## join
The join module is for taking two different scaffolds and joining them into one
scaffold. Common use-cases for this include:
* The scaffolder failed to join two contigs that belong together
* You want to split up the pieces of a scaffold and put them back together
  again in a different order

The two required arguments for this module are a file specifying what joins
you want to make, and the AGP file you want to modify. The joins list contains
one join per line. Each line is a comma-separated list of scaffolds you want
to put together in the correct order. You can prefix the name of a scaffold with
'+' or '-' to specify its orientation; scaffolds with no orientation specified
are '+' by default.

You can also change the default size, type, and evidence of the newly created
gaps with command-line arguments. See help message for details.

Here is an example joins file:
```
scaffold_16.2,-scaffold_16.3,+scaffold_16.1
```

Here is an example command:
```bash
agptools join joins.txt split_out.agp > join_out.agp
```
And here is the output of that command:
```
scaffold_16.2p16.3p16.1 1        7504769    1   W   tig00000113 1   7504769 +
scaffold_16.2p16.3p16.1 7504770  7505269    2   N   500 scaffold    yes na
scaffold_16.2p16.3p16.1 7505270  9508511    3   W   tig00004962 1   2003242 -
scaffold_16.2p16.3p16.1 9508512  9509011    4   N   500 scaffold    yes na
scaffold_16.2p16.3p16.1 9509012  1373456    5   W   tig00004933 1   4225555 -
scaffold_16.2p16.3p16.1 13734567 1373506    6   N   500 scaffold    yes na
scaffold_16.2p16.3p16.1 13735067 1680686    7   W   tig00000080 1   3071803 -
scaffold_16.2p16.3p16.1 16806870 1680736    8   N   500 scaffold    yes na
scaffold_16.2p16.3p16.1 16807370 1754834    9   W   tig00183148 1   740971 -
scaffold_16.2p16.3p16.1 17548341 1754884    1   N   500 scaffold    yes na
scaffold_16.2p16.3p16.1 17548841 1864530    1   W   tig00005080 1   1096465 -
scaffold_16.2p16.3p16.1 18645306 1864580    1   N   500 scaffold    yes na
scaffold_16.2p16.3p16.1 18645806 1952204    1   W   tig00001012 1   876236 +
scaffold_16.2p16.3p16.1 19522042 1952254    1   N   500 scaffold    yes na
scaffold_16.2p16.3p16.1 19522542 2180783    1   W   tig00182876 1   2285293 -
```

You can also specify a new name to use instead of the "16.2p16.3p16.1" scheme.
Add a column to the joins file after a tab giving the new name. For example,
```
scaffold_16.2,-scaffold_16.3,+scaffold_16.1       chr1
```

will result in
```
chr1   1        7504769    1   W   tig00000113 1   7504769 +
chr1   7504770  7505269    2   N   500 scaffold    yes na
chr1   7505270  9508511    3   W   tig00004962 1   2003242 -
chr1   9508512  9509011    4   N   500 scaffold    yes na
chr1   9509012  1373456    5   W   tig00004933 1   4225555 -
chr1   13734567 1373506    6   N   500 scaffold    yes na
chr1   13735067 1680686    7   W   tig00000080 1   3071803 -
chr1   16806870 1680736    8   N   500 scaffold    yes na
chr1   16807370 1754834    9   W   tig00183148 1   740971 -
chr1   17548341 1754884    1   N   500 scaffold    yes na
chr1   17548841 1864530    1   W   tig00005080 1   1096465 -
chr1   18645306 1864580    1   N   500 scaffold    yes na
chr1   18645806 1952204    1   W   tig00001012 1   876236 +
chr1   19522042 1952254    1   N   500 scaffold    yes na
chr1   19522542 2180783    1   W   tig00182876 1   2285293 -
```

## assemble
Once you've got your _final_ final agp, you probably want to create a fasta of
these nice new corrected scaffolds. The `assemble` module takes a fasta file
containing the original contigs and an agp of how you want to assemble these
contigs into scaffolds, and outputs a fasta of the assembled scaffolds. For
example,
```bash
agptools assemble contigs.fa corrected_scaffolds.agp > corrected_scaffolds.fa
```
Please note that if using SALSA as your scaffolder, it makes some breaks to
input contigs based on the Hi-C data, and then gives the pieces different names
(e.g., "contig1" to "contig1\_1" and "contig1\_2"), so you should use the fasta
containing broken contigs as the contigs argument to `assemble` rather than the
actual original contigs you started out with. This file is called
`assembly.clean.fasta` and lives in the same directory as the rest of the SALSA
output.

[agp]: https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
[salsa]: https://github.com/marbl/SALSA
