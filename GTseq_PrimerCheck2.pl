#!/usr/bin/perl
# GTseq_PrimerCheck.pl
# Checks a list of forward and reverse primers for possible problems with primer-dimer artifacts in multiplex PCR.
# Usage: provide a list of primers in the format Name,FWD-Primer,REV-Primer.

use strict; use warnings;

die "Provide a list of Primers\n" unless @ARGV == 1;

my %FWD_Search = ();
my %REV_Search = ();
my %pair_counts = ();
my $tot_ints = 0;

open (FILE, "<$ARGV[0]") or die "Error opening $ARGV[0]\n";

while (<FILE>) {
	chomp;
	my @info = split ",", $_;
	$FWD_Search{$info[0]} = substr $info[1], -9;
	$FWD_Search{$info[0]} =~ tr/ACGT/TGCA/;
	$FWD_Search{$info[0]} = reverse $FWD_Search{$info[0]};
	$REV_Search{$info[0]} = substr $info[2], -9;
	$REV_Search{$info[0]} =~ tr/ACGT/TGCA/;
	$REV_Search{$info[0]} = reverse $REV_Search{$info[0]};
	$pair_counts{$info[0]} = 0;
	#print "$info[0],$FWD_Search{$info[0]},$REV_Search{$info[0]}\n";  #for testing...
	}
close FILE;

sub FindInts {
	
	my $int_tot = 0;  # defines the total number of interactions detected...
	
	foreach my $loci (sort keys %FWD_Search) {
		open (FILE2, "<$ARGV[0]") or die "Error opening $ARGV[0] for reading\n";

		$pair_counts{$loci} = 0;  # initialize number of interactions at a given locus at zero...

		while (<FILE2>) {
			chomp;
			my @info2 = split ",", $_;
			$info2[2] = substr $info2[2], 24;
			$info2[1] = substr $info2[1], 24;
	
			if ($info2[2] =~ m/$FWD_Search{$loci}/) {
				#print "FWD\t$loci\t$FWD_Search{$loci}\tREV\t$info2[0]\t$info2[2]\n";
				$pair_counts{$loci}++;
				$int_tot++;
				#else {print "No Match\n"}
				}
			if ($info2[1] =~ m/$REV_Search{$loci}/) {
				#print "REV\t$loci\t$REV_Search{$loci}\tFWD\t$info2[0]\t$info2[1]\n";
				$pair_counts{$loci}++;
				$int_tot++;
				#else {print "No Match\n"}
				}
			}
		}
	close FILE2;
	$tot_ints = $int_tot;
	}

FindInts();

my $drop = $ARGV[0];
$drop =~ s/.csv$/_LocusDropList.txt/;

open (DROP, ">$drop") or die "Error opening output drop list file $drop\n";

foreach my $names (sort {$pair_counts{$b} <=> $pair_counts{$a}} keys %pair_counts) {
	#print "$names\t$tot_ints\t$pair_counts{$names}\n";
	my $fwd = $FWD_Search{$names};
	my $rev = $REV_Search{$names};
	my $ints = $tot_ints;
	delete $FWD_Search{$names};
	delete $REV_Search{$names};
	my $counts = $pair_counts{$names};  # defines the number of primer interactions for this locus before it's deleted and recalculated
	FindInts();
	if (($tot_ints < $ints) && ($tot_ints > 0)) {
		print "Drop Locus: $names\t$tot_ints\t$pair_counts{$names}\n";
		print DROP "$names\n";
		}
	elsif (($tot_ints == $ints) && ($tot_ints > 0)) {
		$FWD_Search{$names} = $fwd;
		$REV_Search{$names} = $rev;
		}
	elsif ($tot_ints == 0) {
		last;
		}
	}
close DROP;
	
	
