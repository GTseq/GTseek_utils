GTseq Analysis scripts

Use these scripts to test new primer sets to identify primers that are producing primer artifacts and off-target reads

GTseq_PrimerCheck2.pl -- This script tests a set of designed primers and reports any primers that are predicted to produce artifacts in multiplex PCR
execute this script by passing a csv file with primer sequences.  The script only uses the 3' ends of the primer sequences for testing so no need to omit
the Illumina tags.

$ GTseq_PrimerCheck2.pl primers.csv

The script will then filter any primers that are predicted to interact and produce primer artifacts in mulitplex PCR.  
A file with the same name as the primer file with the extension "LocusDropList.txt" will be produced
This file contains locus names that should be omitted or redesigned to avoid primer dimers.

To analyze sequencing data for proper amplification in multiplex PCR, I prefer to use paired-end data.  This allows us to capture the full primer sequence for
both the forward and reverse primers.  For this reason, when I'm analyzing data using the scripts below, I run this script first to produce a file that has both 
R1 and R2 sequences.  

$ paste_fastq.pl R1.fastq R2.fastq > pasted.fq

GTseq_Primer-Interaction-Test_v3.pl  --  This script analyzes sequence data and identifies primers that are producing large numbers of off-target sequences
In order to use this script, your fastq files will first need to be converted to a hash file using the HashSeqs.pl script.  This script collapses all unique 
sequences into a single entry with the number of occurrances.

$ HashSeqs.pl pasted.fq FISH > file.hash  # the second command line argument is just a tag to identify the project 

Then execute the interactions script on the new hash file

$ GTseq_Primer-Interaction-Test_v3.pl file.hash > file_PI-test.txt

The output is a very large file with a TON of information, most of which is useless, so I suggest truncating the file like this.

$ grep -A20 '[Pp]rimer\|[Bb]lack' file_PI-test.txt > file_PI-test-trunc.txt

This will show you several categories of primer artifacts and their relative abundance in the sequencing data.  Keep in mind though that the script is using partial 
primer sequences to make these calculations and sometimes things get counted just because they happen to match something in the sequencing data.  So take the output 
with a grain of salt.  Also, ideally you'll have either paired-end data or the reads are long enough to read through the reverse primer.  This analysis doesn't work
unless the reverse primer sequence is present in the reads.  

The last script I put in this zip folder checks for internal binding sites within the other target amplicons that could produce nested artifact sequences in multiplex
PCR.  You feed this script the designed primer sequences and a fastq file with the sequences used to design the primer set.  The output is a fasta file that report the
expected amplicon sequence when using the designed primers.  It also reports the number of cycles necessary to reach the target SNP so that you could swap the Illumina
tags if it's closer to the target SNP from the other direction.

$ misprime_filter.pl <designed_primer_seqs.csv> <fasta_file_with_target_seqs.fa>  # note that the primer file must be in this format LocusID,fwd-seq,rev-seq.
										     # note that the fasta file names must match exactly the LocusID in the csv file.

Any offending primers will be flagged with a message that says something like "Locus123 is whack yo!".  Which just means that there is an internal binding site for at
least one of the primers for this locus among the other target amplicons.  It also reports which locus has the binding site.
