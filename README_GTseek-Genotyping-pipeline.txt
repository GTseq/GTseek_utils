GT-seq Genotyping Pipeline

This pipeline genotypes single SNP loci from GT-seq library fastq data.  Haplotyping or genotyping of multiple SNPs in an amplicon sequence require using different analyses.  All files should have linux line endings to avoid warnings and error messages.

1-  Sequencing data must be split into individual fastq files using the dual index sequences used for tagging.  This can be accomplished using either the Illumina demultiplexing software OR using the barcode split python script.  
– Illumina demultiplexing:  Populate the “SampleSheet.csv” file with each sampleID and barcode sequence.
– GTseq_BarcodeSplit_MP.py script:  Produce an input file in the format below and execute the script. Note: Dual index sequences must appear on header line of the raw fastq file or the script won’t run.
– GTseq_BarcodeSplit_MP.py script input file format:
SampleID,PlateID,i7_name,i7_seq,i5_name,i5_seq
S0001,P001,i5-001,ACGTGC,i012,TTGCCA
S0002,… etc

Execute the barcode split script as shown
$ GTseq_BarcodeSplit_MP.py
type the path to input file
Format= /home/user/…
<inputfile.csv>
type the path to the fastq file to split
Format= /home/user/…
<Undetermined_S0_R1_001.fastq>
… Script will run and output as many samples as you have lines in the input file (not including the header line).  The header line is ignored by the script, so include it in all input files.  Note that workflow B instruments will capture the reverse complement of the 5’-3’ index sequence in I2 (P5 – index seq).  Check that each individual file has been captured as expected using a simple line count command.
$ wc -l *fastq
– each sequence read has 4 lines in fastq format.  Samples should return approximately 100,000 reads per individual but there will be a distribution around this number.
– If everything looks good, move all the individual fastq files into their own directory for genotyping.

2-  Genotype the individual fastq files using GTseq genotyper script
$ GTseq_Genotyper_v3.pl <input_file.csv> <sample.fastq> > sample.genos
– This must be run for each individual sample.  To make it easier to execute, set up a shell script to process all the samples using a “one-liner” command.
$ ls *fastq | sed ‘s/.*/GTseq_Genotyper_v3.pl  <path to input file> & > &/’ | sed ‘s/fastq$/genos/’ > genotyper.sh
$ sh genotyper.sh
– once finished, you should have one .genos file for every .fastq file in the directory.

3-  At this point, all the samples have been genotyped.  You may wish to compile them into a single file though.  To output a .csv with all the genotype results, you can use the GTseq_GenoCompile_v3.pl script.
$ GTseq_GenoCompile_v3.pl 
– This will output every captured genotype at each locus.
If you wish to suppress output or output in different formats, there are options for that using this script.
$ GTseq_GenoCompile_v3.pl C 90   # this command will output read count for each locus but will only report for samples with call rates at or above 90%.  
– Output types:  S = SNP; C = counts; N = numeric genotypes; A = allele counts
– Genotype capture percentage is a second argument you can use to filter data.  It will only report no-calls for any sample with call-rates below the set threshold.

3-  You may want to check the quality of the data or the overall performance of the run.  For this, I’ve written a python script to summarize the data for a given library.  This script uses the .genos files that were output by the GTseq_Genotyper_v3.pl script to produce summary figures for the run.  The output is a series of .png image files of the summary figures.  These may be compiled into a pdf using imagemagick’s convert command.
$ GTseq_SummaryFigures_v3.pl 
type the path to directory containing .genos files for library
./
type the library name
LIB001
## the script will take a little while to run ##
Once finished, use the convert command to produce a summary pdf for the run.
$ convert *png LIB001.pdf

