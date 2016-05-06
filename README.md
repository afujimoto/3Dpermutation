# 3Dpermutation
This is a program to test distribution of somatic mutations in 3D protein structure.

Please reference our manuscript whenever you utilize this package in your analysis.

Fujimoto A, Okada Y, Boroevich KA, Tsunoda T and Taniguchi H, and Nakagawa H. (2016) Systematic analysis of mutation distribution in three dimensional protein structures identifies cancer driver genes. Sci Rep (accepted)


##Requirement

This program needs two software and data files.

(1) samtools (software)

samtools; http://samtools.sourceforge.net


(2) MAFFT (software)

MAFFT; http://mafft.cbrc.jp/alignment/software/


(3) Gene annotation file

Downloadable from the same repository (Annotation.gencode.v19.txt)

If you want to use a different annotation information, please make gene annotation file by yourself.

(4) Indexed reference genome sequence file

Please download reference genome sequence file and create index file by "samtools index <reference genome sequence file>".

(5) Protein structure file (PDB format)

Please download protein structure file from "http://www.rcsb.org/pdb/home/home.do".


##Usage

perl 3Dpermutation.pl -N \<Minimum number of mutations in the 3D structure\> -G \<Transcript ID\> -I \<PDB file\>  -M \<Mutation position on the transcript (Amino Acid position)\> -C \<Chain of the protein in the PDB file\> -P \<Number of permutations\> -R \<Reference genome sequence file\> -A \<Gene annotation file\>

*Output format

gene_symbol; Gene symbol of the transcript ID


chain; Chain of the protein in the PDB file


mRNA_id; Transcript ID


PDB_file; File name of the PDB file


Length_ratio; Ratio of aligned sequence length to total sequence length


Gap_ratio; Ratio of gaps to length of the alignment 


Number_of_mapped_mutation; Number of mapped mutations on the amino acid sequence in the 3D structure


Mutation_position; Position of mutations


Mutation_position_on_protein; Position of mutations on the PDB sequence


Average_distance_between_mutations; Average distance between mutations


Permutaiton_P_value; P-value obtained by the permutation test


Two fasta files of the amino acid sequences will be generated in the current directory.


##Example

The following command performs 3D permutation analysis for five mutations (11, 24, 27, 138 and 142 in amino acid position) in PTEN gene (transcript id ENST00000371953.3) for 3D structure in "pdb1d5r.ent" file. Permutation repeats 10000 times.

perl 3Dpermutation.pl -G ENST00000371953.3 -I pdb1d5r.ent  -M 11:24:27:138:142 -C A -N 3 -P 10000 -R All.fa -A Annotation.gencode.v19.txt > result.txt


##Contact

Akihiro Fujimoto - fujimoto@ddm.med.kyoto-u.ac.jp