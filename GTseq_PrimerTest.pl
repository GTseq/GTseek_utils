#!/usr/bin/perl
# GTseq_PrimerTest.pl
# by Nate Campbell
# count the occurrences of primers in GTseq libraries and do some calculations

use strict; use warnings;

my %Fprimers = ();
my %Rprimers = ();
my %a1probe = ();
my %a2probe = ();
my %a1rc = ();
my %a2rc = ();
my %f_primer_ct = ();
my %r_primer_ct = ();
my %probe_ct = ();
my %both_ct = ();
my %all_together = ();

die "provide a list of fwd, rev, and probe seqs and R1.fastq and R2.fastq files\n" unless @ARGV == 3;

`paste $ARGV[1] $ARGV[2] > pasted_temp.fq`;
my $file = "pasted_temp.fq";  #create temporary file with R1 and R2 fastq data...

open (PRIMERS, "<$ARGV[0]") or die "Error opening primer file\n";

while (<PRIMERS>) {
	chomp;
	my @stuff = split "\t", $_;
	my $fpr = substr $stuff[1], 0, 18;
	my $rpr = substr $stuff[2], 0, 18;
	$Fprimers{$fpr} = $stuff[0];
	$Rprimers{$rpr} = $stuff[0];
	$a1probe{$stuff[0]} = $stuff[3];
	$a2probe{$stuff[0]} = $stuff[4];
	$a1rc{$stuff[0]} = reverse $stuff[3];
	$a1rc{$stuff[0]} =~ tr/ACGT/TGCA/;
	$a1rc{$stuff[0]} =~ tr/\]\[/\[\]/;
	$a2rc{$stuff[0]} = reverse $stuff[4];
	$a2rc{$stuff[0]} =~ tr/ACGT/TGCA/;
	$a2rc{$stuff[0]} =~ tr/\]\[/\[\]/;
	$f_primer_ct{$stuff[0]} = 0;
	$r_primer_ct{$stuff[0]} = 0;
	$probe_ct{$stuff[0]} = 0;
	$both_ct{$stuff[0]} = 0;
	$all_together{$stuff[0]} = 0;
	}
close PRIMERS;

open (FASTQ, "<$file") or die "Error opening fastq file\n";

while (<FASTQ>) {
	my $header = $_;
	my $seq = <FASTQ>;
	my $h2 = <FASTQ>;
	my $qual = <FASTQ>;
	chomp($seq);
	my $primer_hit = 0;
	my @r1n2 = split "\t", $seq;
	my $test_seq1 = substr $r1n2[0], 0, 18;
	if (exists $Fprimers{$test_seq1}) {
		$f_primer_ct{$Fprimers{$test_seq1}}++;
		if (($seq =~ m/$a1probe{$Fprimers{$test_seq1}}/) or ($seq =~ m/$a2probe{$Fprimers{$test_seq1}}/)) {
			$probe_ct{$Fprimers{$test_seq1}}++;
			$primer_hit++;
			}
		elsif (($seq =~ m/$a1rc{$Fprimers{$test_seq1}}/) or ($seq =~ m/$a2rc{$Fprimers{$test_seq1}}/)) {
			$probe_ct{$Fprimers{$test_seq1}}++;
			$primer_hit++;
			}
		}
	my $test_seq2 = substr $r1n2[1], 0, 18;
	if (exists $Rprimers{$test_seq2}) {
		$r_primer_ct{$Rprimers{$test_seq2}}++;
		$primer_hit++;
		}
	if ($primer_hit == 2) {
		$both_ct{$Rprimers{$test_seq2}}++;
		if (($seq =~ m/$a1probe{$Fprimers{$test_seq1}}/) or ($seq =~ m/$a2probe{$Fprimers{$test_seq1}}/)) {
			$all_together{$Fprimers{$test_seq1}}++;
			}
		elsif (($seq =~ m/$a1rc{$Fprimers{$test_seq1}}/) or ($seq =~ m/$a2rc{$Fprimers{$test_seq1}}/)) {
			$all_together{$Fprimers{$test_seq1}}++;
			}
		}
	}
close FASTQ;

`rm pasted_temp.fq`;  #remove temporary file...

print "LocusID\tFWD\tREV\tPAIRED\tPROBES\tPAIRED+PROBE\tPrimer Agreement\n";

foreach my $locus (sort keys %f_primer_ct) {
	if ($f_primer_ct{$locus} == 0) {
		$f_primer_ct{$locus} = 1;
		}
	my $prim_agree = 100 - (abs($f_primer_ct{$locus} - $r_primer_ct{$locus})/($f_primer_ct{$locus} + $r_primer_ct{$locus}))*100;
	$prim_agree = sprintf("%.2f", $prim_agree);
	print "$locus\t$f_primer_ct{$locus}\t$r_primer_ct{$locus}\t$both_ct{$locus}\t$probe_ct{$locus}\t$all_together{$locus}\t$prim_agree\%\n";
	}
	
