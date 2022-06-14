#!/usr/bin/perl
#GTseq_Primer-Interaction-Test.pl
#Test .hash file to count Primer Pair sequences in readthroughs
# Provide a tab delimited file of LocusName\tFWD-Primer sequence\tREV-Primer sequence\n and a .hash file

use strict; use warnings;

die "usage: provide <tab delimited txt file of locus\tfwd\trev seqs> and <*.hash>\n" unless @ARGV == 2;

my $raw_reads = 0;
my $ot = 0;
my $dc = 0;
my $fr = 0;
my $ff = 0;
my $rf = 0;
my $rr = 0;
my @assays = ();
my %fwd_seq = ();
my %fwd_RC = ();
my %fwd_start = ();
my %rev_seq = ();
my %rev_RC = ();
my %rev_start = ();
my %rev_RC_Full = ();
my %OT_PrimerComb = ();
my %DC_Artifact = ();
my %PrimerComb = ();
my %PrimerCombFW = ();
my %PrimerCombRV = ();
my %PrimerCombRV_RV = ();
my %LocusCount = ();
my %PercentRaw = ();
my %hash_check = ();  # holds a value if either a forward or reverse primer has been identified within the sequence
my $last_hash_id = 0;  # holds the value of the highest hashID

#read in assay and allele information and push to arrays...

open(SEQ, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<SEQ>) {
	chomp;
	my @info = split(/\t/, $_);
	push @assays, $info[0];
	$fwd_seq{$info[0]} = $info[1];
	$fwd_start{$info[0]} = substr $fwd_seq{$info[0]}, 0, 10;
	$fwd_RC{$info[0]} = reverse $fwd_seq{$info[0]};
	$fwd_RC{$info[0]} =~ tr/ACGT/TGCA/;
	$fwd_RC{$info[0]} = substr $fwd_RC{$info[0]}, -10;
	$rev_seq{$info[0]} = $info[2];
	$rev_start{$info[0]} = substr $rev_seq{$info[0]}, 0, 10;
	$rev_RC{$info[0]} = reverse $rev_seq{$info[0]};
	$rev_RC{$info[0]} =~ tr/ACGT/TGCA/;
	$rev_RC_Full{$info[0]} = $rev_RC{$info[0]};
	$rev_RC{$info[0]} = substr $rev_RC{$info[0]}, -10;
	#print "$info[0];$fwd_seq{$info[0]};$rev_seq{$info[0]};$fwd_RC{$info[0]};$rev_RC{$info[0]};$fwd_start{$info[0]};$rev_start{$info[0]};$rev_RC_Full{$info[0]}\n";  #Testing...
		}
close SEQ;

my $iteration = 0;

foreach my $Locus (@assays) {

$iteration++;

open (HASH, "<$ARGV[1]") or die "error reading $ARGV[1]\n";

	while (<HASH>) {
		my $hash = $_;
		my $R1seq = <HASH>;
		chomp ($hash);
		chomp ($R1seq);
		my @info = split(/;/, $hash);
		my $count = $info[2];
		my $hashID = $info[1];
		$last_hash_id = $hashID;
		if ($iteration == 1) {$raw_reads = $raw_reads + $count;}
		if ($R1seq =~ m/^$fwd_seq{$Locus}/) {
			$hash_check{$hashID} = "Yes";
			foreach my $Revs (@assays) {
				if (($R1seq =~ m/$rev_RC_Full{$Revs}/) && ($Locus !~ m/$Revs/)) {
					$dc = $dc + $count;
					if (exists $LocusCount{$Locus}) {$LocusCount{$Locus} = $LocusCount{$Locus} + $count;}
					else {$LocusCount{$Locus} = $count;}
					if (exists $LocusCount{$Revs}) {$LocusCount{$Revs} = $LocusCount{$Revs} + $count;}
					else {$LocusCount{$Revs} = $count;}
					my $pair = "$Locus $Revs";
					if (exists $DC_Artifact{$pair}) {$DC_Artifact{$pair} = $DC_Artifact{$pair} + $count;}
					else {$DC_Artifact{$pair} = $count;}
					}
				elsif (($R1seq =~ m/$rev_RC{$Revs}/) && ($Locus !~ m/$Revs/)) {
					$fr = $fr + $count;
					if (exists $LocusCount{$Locus}) {$LocusCount{$Locus} = $LocusCount{$Locus} + $count;}
					else {$LocusCount{$Locus} = $count;}
					if (exists $LocusCount{$Revs}) {$LocusCount{$Revs} = $LocusCount{$Revs} + $count;}
					else {$LocusCount{$Revs} = $count;}
					my $pair = "$Locus $Revs";
					if (exists $PrimerComb{$pair}) {$PrimerComb{$pair} = $PrimerComb{$pair} + $count;}
					else {$PrimerComb{$pair} = $count;}
					}
				elsif (($R1seq =~ m/$rev_RC{$Revs}/) && ($Locus =~ m/$Revs/)) {
					$ot = $ot + $count;
					if (exists $LocusCount{$Locus}) {$LocusCount{$Locus} = $LocusCount{$Locus} + $count;}
					else {$LocusCount{$Locus} = $count;}
					if (exists $LocusCount{$Revs}) {$LocusCount{$Revs} = $LocusCount{$Revs} + $count;}
					else {$LocusCount{$Revs} = $count;}
					my $pair = "$Locus $Revs";
					if (exists $OT_PrimerComb{$pair}) {$OT_PrimerComb{$pair} = $OT_PrimerComb{$pair} + $count;}
					else {$OT_PrimerComb{$pair} = $count;}
					}
				}
			foreach my $FWDRC (@assays) {
				if (($R1seq =~ m/$fwd_RC{$FWDRC}/) && ($Locus !~ m/$FWDRC/)) {
					$ff = $ff + $count;
					if (exists $LocusCount{$Locus}) {$LocusCount{$Locus} = $LocusCount{$Locus} + $count;}
					else {$LocusCount{$Locus} = $count;}
					if (exists $LocusCount{$FWDRC}) {$LocusCount{$FWDRC} = $LocusCount{$FWDRC} + $count;}
					else {$LocusCount{$FWDRC} = $count;}
					my $pair = "$Locus $FWDRC";
					if (exists $PrimerCombFW{$pair}) {$PrimerCombFW{$pair} = $PrimerCombFW{$pair} + $count;}
					else {$PrimerCombFW{$pair} = $count;}
					}
				}
			}
		elsif ($R1seq =~ m/$rev_RC_Full{$Locus}/) {
			$hash_check{$hashID} = "Yes";
			foreach my $FWD_starts (@assays) {
				if (($R1seq =~ m/$fwd_start{$FWD_starts}/) && ($Locus !~ m/$FWD_starts/)) {
					$rf = $rf + $count;
					if (exists $LocusCount{$Locus}) {$LocusCount{$Locus} = $LocusCount{$Locus} + $count;}
					else {$LocusCount{$Locus} = $count;}
					if (exists $LocusCount{$FWD_starts}) {$LocusCount{$FWD_starts} = $LocusCount{$FWD_starts} + $count;}
					else {$LocusCount{$FWD_starts} = $count;}
					my $pair = "$Locus $FWD_starts";
					if (exists $PrimerCombRV{$pair}) {$PrimerCombRV{$pair} = $PrimerCombRV{$pair} + $count;}
					else {$PrimerCombRV{$pair} = $count;}
					}
				}
			foreach my $revs2 (@assays) {
				if (($R1seq =~ m/$rev_start{$revs2}/) && ($Locus !~ m/$revs2/)) {
					$rr = $rr + $count;
					if (exists $LocusCount{$Locus}) {$LocusCount{$Locus} = $LocusCount{$Locus} + $count;}
					else {$LocusCount{$Locus} = $count;}
					if (exists $LocusCount{$revs2}) {$LocusCount{$revs2} = $LocusCount{$revs2} + $count;}
					else {$LocusCount{$revs2} = $count;}
					my $pair = "$Locus $revs2";
					if (exists $PrimerCombRV_RV{$pair}) {$PrimerCombRV_RV{$pair} = $PrimerCombRV_RV{$pair} + $count;}
					else {$PrimerCombRV_RV{$pair} = $count;}
					}
				}
			}
		
		}
		close HASH;
	}

my %BlackSeqs = ();  # holds information for sequences in which primer combinations can't be identified
my %BlackCounts = (); #holds read counts for each black seq...
my $BlackReads = 0; # holds read counts for all black reads...
my $PercentBlack = 0;  # holds the percentage of black sequences in the hash file...


# check all hashIDs to see if they have primers identified...
open (HASH2, "<$ARGV[1]") or die "Error opening hash file\n";
while (<HASH2>) {
	my $hash = $_;
	my $R1seq = <HASH2>;
	chomp ($hash);
	chomp ($R1seq);
	my @info = split(/;/, $hash);
	my $count = $info[2];
	my $hashID = $info[1];
	if (exists $hash_check{$hashID}) {}
	else {
		$BlackSeqs{$hashID} = $R1seq;
		$BlackCounts{$hashID} = $count;
		$BlackReads = $BlackReads + $count;
		#print "$hashID\tnot found\n";
		}
	}
close HASH2;	

print "Locus\tCounts with locus primers\tPercentage of raw reads attributible\n";
#sort all loci by number of attributible reads...
foreach my $targets (sort { $LocusCount{$b} <=> $LocusCount{$a} } keys %LocusCount) {
	$PercentRaw{$targets} = $LocusCount{$targets}/($raw_reads*2)*100;
	$PercentRaw{$targets} = sprintf("%.2f",$PercentRaw{$targets});
	print "$targets\t$LocusCount{$targets}\t$PercentRaw{$targets}\%\n";
	}

#print read counts from detected primer interactions...
my $Percentage = sprintf("%.2f",$ot/$raw_reads*100);
print "***Proper On-Target Primer combinations (%total reads: $Percentage)***\n";
foreach my $Combos1 (sort { $OT_PrimerComb{$b} <=> $OT_PrimerComb{$a} } keys %OT_PrimerComb) {
	$PercentRaw{$Combos1} = $OT_PrimerComb{$Combos1}/$raw_reads*100;
	$PercentRaw{$Combos1} = sprintf("%.2f",$PercentRaw{$Combos1});
	print "$Combos1\t$OT_PrimerComb{$Combos1}\t$PercentRaw{$Combos1}\%\n";
	}
$Percentage = sprintf("%.2f",$dc/$raw_reads*100);
print "***Double-Complement Primer-Artifacts (%total reads: $Percentage)***\n";
foreach my $CombosDC (sort { $DC_Artifact{$b} <=> $DC_Artifact{$a} } keys %DC_Artifact) {
	$PercentRaw{$CombosDC} = $DC_Artifact{$CombosDC}/$raw_reads*100;
	$PercentRaw{$CombosDC} = sprintf("%.2f",$PercentRaw{$CombosDC});
	print "$CombosDC\t$DC_Artifact{$CombosDC}\t$PercentRaw{$CombosDC}\%\n";
	}
$Percentage = sprintf("%.2f",$fr/$raw_reads*100);
print "****FWD-REV Primer mis-primes (%total reads: $Percentage)****\n";
foreach my $Combos (sort { $PrimerComb{$b} <=> $PrimerComb{$a} } keys %PrimerComb) {
	$PercentRaw{$Combos} = $PrimerComb{$Combos}/$raw_reads*100;
	$PercentRaw{$Combos} = sprintf("%.2f",$PercentRaw{$Combos});
	print "$Combos\t$PrimerComb{$Combos}\t$PercentRaw{$Combos}\%\n";
	}
$Percentage = sprintf("%.2f",$ff/$raw_reads*100);
print "***FWD-FWD Primer mis-primes (%total reads: $Percentage)***\n";
foreach my $Combos2 (sort { $PrimerCombFW{$b} <=> $PrimerCombFW{$a} } keys %PrimerCombFW) {
	$PercentRaw{$Combos2} = $PrimerCombFW{$Combos2}/$raw_reads*100;
	$PercentRaw{$Combos2} = sprintf("%.2f",$PercentRaw{$Combos2});
	print "$Combos2\t$PrimerCombFW{$Combos2}\t$PercentRaw{$Combos2}\%\n";
	}
$Percentage = sprintf("%.2f",$rf/$raw_reads*100);
print "***REV-FWD Primer mis-primes (%total reads: $Percentage)***\n";
foreach my $Combos3 (sort { $PrimerCombRV{$b} <=> $PrimerCombRV{$a} } keys %PrimerCombRV) {
	$PercentRaw{$Combos3} = $PrimerCombRV{$Combos3}/$raw_reads*100;
	$PercentRaw{$Combos3} = sprintf("%.2f",$PercentRaw{$Combos3});
	print "$Combos3\t$PrimerCombRV{$Combos3}\t$PercentRaw{$Combos3}\%\n";
	}
$Percentage = sprintf("%.2f",$rr/$raw_reads*100);
print "***REV-REV Primer mis-primes (total reads: $Percentage)***\n";
foreach my $Combos4 (sort { $PrimerCombRV_RV{$b} <=> $PrimerCombRV_RV{$a} } keys %PrimerCombRV_RV) {
	$PercentRaw{$Combos4} = $PrimerCombRV_RV{$Combos4}/$raw_reads*100;
	$PercentRaw{$Combos4} = sprintf("%.2f",$PercentRaw{$Combos4});
	print "$Combos4\t$PrimerCombRV_RV{$Combos4}\t$PercentRaw{$Combos4}\%\n";
	}

# print read counts from sequences for which primer combinations can't be identified...
$PercentBlack = sprintf("%.2f", $BlackReads/$raw_reads*100);
print "*** Black Reads (total reads: $PercentBlack) ***\n";
foreach my $blacks (sort { $a <=> $b } keys %BlackSeqs) {
	print "$blacks\t$BlackCounts{$blacks}\t$BlackSeqs{$blacks}\n";
	}

