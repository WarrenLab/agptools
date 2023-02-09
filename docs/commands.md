# agptools commands

agptools's command-line interface consists of several sub-commands, in the
style of, e.g., samtools or git. Here we document each of the commands
including an example of their usage.

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

## remove
You may want to get rid of some scaffolds because they contain contamination,
duplication, or something else. The remove module can help with that. Just
give it a list of scaffolds you don't want in the final assembly and it will
remove them. The list of scaffolds to remove should have one scaffold per line,
e.g.,
```
scaffold_5
scaffold_7
```
Running this command:
```bash
agptools remove scaffolds_to_remove.txt scaffolds.agp > corrected_scaffolds.agp
```
would output `scaffolds.agp` but with all lines corresponding to `scaffolds_5`
and `scaffolds_7` removed.

## rename
Often, you end up with scaffolds that correspond to whole chromosomes, and you
want to therefore give them names befitting chromosomes. The input file for
this command has two required columns and an optional third one:

1. Current name of scaffold
2. New name of scaffold
3. Orientation, either `+` or `-` (optional). If this column is `-`, the new
   scaffold will be reverse-oriented compared to the input. If this column is
   left blank, the new scaffold will be exactly the same as the old one, just
   with a different name.

## transform
You may have bed files of alignments or annotations where the coordinates are
given in reference to the original contigs, but you want those coordinates
transformed to your new scaffolds. This module makes those transformations.
For example, if you had this bed file:
```
tig00005080   10952    10960
tig00004962   1        2003242
```
and wanted to convert it to scaffold coordinates, you would run the command:
```bash
agptools transform contig_coordinates.bed example.agp > scaffold_coordinates.bed
```

and the output would be
```
scaffold_16   10952     10960
scaffold_16   11764764  13768005
```

## sanitize
NCBI has some pretty strict rules about how genomes in the contig fasta + agp
format combination should be submitted:
* There must be a separate agp for chromosomes and for unplaced sequences.
* Unplaced single-component scaffolds must have their single component in the
  '+' orientation.
* Each sequence in the contigs file must be used fully and contiguously in the
  AGP layout.
