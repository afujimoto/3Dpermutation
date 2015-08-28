# 3Dpermute
This is a program to test distances between somatic mtautions in 3D protein structure.

Please reference our manuscripts whenever you utilize this package in your analysis.

Fujimoto A. et al. Systematic analysis of mutation distribution in three dimensional protein structures identifies cancer driver genes (submitted)

*Setting

This program needs two software and data files.

(1) samtools (software)

samtools; http://samtools.sourceforge.net


(2) MAFFT (software)

MAFFT; http://mafft.cbrc.jp/alignment/software/


(3) Gene annotaiton file

Downloadable from the same repository (Annotation.gencode.v19.txt)


(4) Indexed reference genome sequence file

Please download reference genome sequence file and indexed by "samtools index <reference genome sequence file>".

(5) Protein strucuture file (PDB format)

Please download protein structure file from "http://www.rcsb.org/pdb/home/home.do".


*Usage

perl 3Dpermute3.pl -D -G <Transcript ID> -I <PDB file>  -M <Mutation position on the transcript (Amino Acid position)> -C <Chain of PDB> -P <Number of permutation> -R <Reference genome sequence file> -A <Gene annotaiton file>

example
The following command is 3D permutation analysis for five mutations (11, 24, 27, 138 and 142 in amino acid position) of PTEN gene (transcript id ENST00000371953.3) for 3D structure in "pdb1d5r.ent" file.

perl 3Dpermute.pl -G ENST00000371953.3 -I pdb1d5r.ent  -M 11:24:27:138:142 -C A -N 3 -P 10000 -R All.fa -A Annotation.gencode.v19.txt > result.txt

This is a new soure repository, currently under construction!
29/7/2015

