#!/usr/bin/perl
#Omy192_SeqTest.pl
#Test .hash file to see how many times loci fwd primer plus probe seqs occur

use strict; use warnings;

die "usage: provide <ProbeSeqs file> and <*.hash>\n" unless @ARGV == 2;

my %fwd_seq = ();
my %fwd_seqAssay = ();
my %probe1 = ();
my %probe2 = ();
my %probe1RC = ();
my %probe2RC = ();
my %fwd_count = ();
my %probe_count = ();
my %both_count = ();
my %OT_Ratio = ();

#read in assay and allele information and push to arrays...

open(SEQ, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<SEQ>) {
	chomp;
	my @info = split(/,/, $_);
	my $tring = substr $info[5], 0, 14;
	$fwd_seq{$info[0]} = $tring;
	$fwd_seqAssay{$tring} = $info[0];
	$probe1{$info[0]} = $info[3];
	$probe2{$info[0]} = $info[4];
	my $p1RC = reverse $info[3];
	my $p2RC = reverse $info[4];
	$p1RC =~ tr/ACGT][/TGCA[]/;
	$p2RC =~ tr/ACGT][/TGCA[]/;
	$probe1RC{$info[0]} = $p1RC;
	$probe2RC{$info[0]} = $p2RC;
		}
close SEQ;

#initialize stats at zero...
foreach my $assays (sort keys %fwd_seq){
$fwd_count{$assays} = 0;
$probe_count{$assays} = 0;
$both_count{$assays} = 0;
$OT_Ratio{$assays} = 0;
}
open(HASH, "<$ARGV[1]") or die "error reading $ARGV[1]\n";

	while (<HASH>) {
		my $hash = $_;
		my $R1_seq = <HASH>;
		my $seq_start = substr $R1_seq, 0, 14;
		chomp ($hash);
		chomp ($R1_seq);
		my $RC_seq = reverse $R1_seq;
		$RC_seq =~ tr/ACGT/TGCA/;
		#print "$R1_seq\n$RC_seq\n";   #testing...
		my $seq_start2 = substr $RC_seq, 0, 14;
		#print "$seq_start2\n";
		my @info = split(/;/, $hash);
		my $count = $info[2];
		my $count2 = $info[2];

			if (exists $fwd_seqAssay{$seq_start}) {
				my $assay_hit = $fwd_seqAssay{$seq_start};
				$count = $fwd_count{$assay_hit} + $count;
				$fwd_count{$assay_hit} = $count;
				$count = $info[2];

				if ($R1_seq =~ m/$probe1{$assay_hit}|$probe2{$assay_hit}|$probe1RC{$assay_hit}|$probe2RC{$assay_hit}/){
				$count = $probe_count{$assay_hit} + $count;
				$probe_count{$assay_hit} = $count;
				}
				$count = $info[2];
				if (($R1_seq =~ m/^$fwd_seq{$assay_hit}/) && ($R1_seq =~ m/$probe1{$assay_hit}|$probe2{$assay_hit}|$probe1RC{$assay_hit}|$probe2RC{$assay_hit}/)) {
				$count = $both_count{$assay_hit} + $count;
				$both_count{$assay_hit} = $count;
				}
			}
			if (exists $fwd_seqAssay{$seq_start2})  {
				my $assay_hit2 = $fwd_seqAssay{$seq_start2};
				$count2 = $fwd_count{$assay_hit2} + $count2;
				$fwd_count{$assay_hit2} = $count2;
				$count2 = $info[2];

				if ($RC_seq =~ m/$probe1{$assay_hit2}|$probe2{$assay_hit2}|$probe1RC{$assay_hit2}|$probe2RC{$assay_hit2}/){
				$count2 = $probe_count{$assay_hit2} + $count2;
				$probe_count{$assay_hit2} = $count2;
				}
				$count2 = $info[2];
				if (($RC_seq =~ m/^$fwd_seq{$assay_hit2}/) && ($RC_seq =~ m/$probe1{$assay_hit2}|$probe2{$assay_hit2}|$probe1RC{$assay_hit2}|$probe2RC{$assay_hit2}/)) {
				$count2 = $both_count{$assay_hit2} + $count2;
				$both_count{$assay_hit2} = $count2;
				}
			}
	}
close HASH;

print "LOCUS,FWD_PRIMER,PROBE,BOTH,OT_RATIO\n";

foreach my $assays2 (sort keys %fwd_seq) {
	if ($fwd_count{$assays2} == 0) {$fwd_count{$assays2} = 0.1}
	$OT_Ratio{$assays2} = $probe_count{$assays2}/$fwd_count{$assays2};
	print "$assays2,$fwd_count{$assays2},$probe_count{$assays2},$both_count{$assays2},$OT_Ratio{$assays2}\n";
}

