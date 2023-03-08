# motif-mark

Kobe Ikegami
08MAR2023

Repository Description:

Inputs/Argparse arguments:
- FASTA file with sequences split over multiple lines and exon nucleotides are capitalized within the sequenced, generated in the UCSC Genome Browser
- Motif file: a list of motifs to be referenced in a search of the sequence lines of the FASTA file

Intermediate file:
- A One-line FASTA file is generated where all sequenced lines are concatenated into a single line for easier parsing.

Outputs:
- png
- pdf

Argparse Arguments:
- -f: input FASTA file
- -m: motif file
- -o: output png
- -i: intermediate FASTA file

Packages Used:
- argrparse
- pycairo
- re
- math

Unit Tests:
- test Fasta file
- test motif file

Code Function:
The motif-mark-oop.py uses object oreinted programming to take a multi-line FASTA file with capitalized exonic bases and a file containing a list of motifs and return a png annotating the gene sequences of the FASTA file. It creates a color-coded image that depicts Genes, their exons, and all of their motifs of each gene. The script is capable of handling 5 motifs and sequences < 1000 in length. 
