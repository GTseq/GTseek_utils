#!/usr/bin/perl
# GTseq_PrimerAnalysis.pl
# by Nate Campbell
# use a list of forward and reverse primer sequences to evaluate loci in GT-seq panel
# Requires that libraries are run as paired end and pasted using paste_fastq.pl script

use strict; use warnings;
use List::MoreUtils qw(uniq);

my %fwd_primer = ();
my %rev_primer = ();
my %fprimer = ();
my %rprimer = ();

# hashes for overall performance by sample and locus...
my %samp_otp = ();  # sample on-target percentage
my %samp_raw = ();  # sample raw reads
my %samp_ot = ();   # sample on-target reads
my %loc_otp = ();   # locus percentage of on-target reads
my %loc_praw = ();  # locus percentage of raw reads
my %sig_seqs = ();  # stores unique sequences captured by the correct fwd-rev primers 
my %loc_psv = ();   # stores the number of times a locus (on-targets) reports more than 2 unique sequences
my %str_loc_otp = ();  # stores the locus on-target percentage for each sample
my %str_loc_praw = (); # stores the locus raw read percentage for each sample
my %str_loc_potp = ();  # stores a string of the primer on-target percentage per locus
my %str_loc_alct = ();  # stores a string of on-target counts per locus at each sample

die "Please provide a tab separated txt file with locusID fwd_primer rev_primer\n" unless @ARGV == 1;

open (LIST, "<$ARGV[0]") or die "Error opening $ARGV[0]\n";

while (<LIST>) {
	chomp;
	my @info = split "\t", $_;
	$fwd_primer{$info[0]} = substr $info[1], 0, 15;
	$rev_primer{$info[0]} = substr $info[2], 0, 15;
	$rev_primer{$info[0]} = reverse $rev_primer{$info[0]};
	$rev_primer{$info[0]} =~ tr/ACGT/TGCA/;
	$loc_psv{$info[0]} = 0;  # initialize variable
	$fprimer{$fwd_primer{$info[0]}} = $info[0];
	$rprimer{$rev_primer{$info[0]}} = $info[0];
	#print "$info[0]\t$fwd_primer{$info[0]}\t$rev_primer{$info[0]}\n";
	}
close LIST;

my @fastq = `ls *fastq`;
chomp(@fastq);

foreach my $files (@fastq) {
	my $raw_reads = 0;  # initialize counts for raw reads
	my $ot_reads = 0;  # initialize counts for on-target reads
	open (FASTQ, "<$files") or die "Error opening $files\n";
	my $name = $files;
	$name =~ s/.fastq//;
	my %fwd_counts = ();
	my %rev_counts = ();
	my %sequences = ();
	while (<FASTQ>) {
		my $header = $_;
		my $seq = <FASTQ>;
		my $h2 = <FASTQ>;
		my $qual = <FASTQ>;
		chomp($seq);
		my $fwd = substr $seq, 0, 15;
		$raw_reads++;
		#print "$fwd\n";
		if (exists $fprimer{$fwd}) {
			if (exists $fwd_counts{$fprimer{$fwd}}) {
				$fwd_counts{$fprimer{$fwd}}++;
				}
			else {
				$fwd_counts{$fprimer{$fwd}} = 1;
				}
			if ($seq =~ m/$rev_primer{$fprimer{$fwd}}/) {
				$ot_reads++;
				if (exists $rev_counts{$fprimer{$fwd}}) {
					$rev_counts{$fprimer{$fwd}}++;
					my $cap_seq = substr $seq, 0, 75;
					$sequences{$fprimer{$fwd}} = "$sequences{$fprimer{$fwd}},$cap_seq";
					}
				else {
					$rev_counts{$fprimer{$fwd}} = 1;
					my $cap_seq = substr $seq, 0, 75;
					$sequences{$fprimer{$fwd}} = $cap_seq;
					}
				}
			}
		}
	close FASTQ;
	foreach my $locs (sort keys %fwd_primer) {
		my $psv_seen = 0;  # initialize counts for when 3 or more unique sequences are captured by the primers & each represents at least 10% of the captured seqs
		if (exists $rev_counts{$locs}) {} #accounting for primers that weren't encountered...
		else {$rev_counts{$locs} = 0.001}
		if (exists $fwd_counts{$locs}) {}
		else {$fwd_counts{$locs} = 10}
		my $p_otp = $rev_counts{$locs}/$fwd_counts{$locs} * 100;
		$p_otp = sprintf("%.2f", $p_otp);
		if (exists $str_loc_potp{$locs}) {
			$str_loc_potp{$locs} = "$str_loc_potp{$locs},$p_otp";
			}
		else {
			$str_loc_potp{$locs} = $p_otp;
			}
		my $overall_otp = $ot_reads/$raw_reads*100;  #overall on-target percentage per sample
		$overall_otp = sprintf("%.2f", $overall_otp);
		my $p_loc_raw = $rev_counts{$locs}/$raw_reads*100;  #percentage of raw reads attributable to to locus
		$p_loc_raw = sprintf("%.2f", $p_loc_raw);
		if ($ot_reads == 0) {$ot_reads = 1}
		my $p_loc_ot = $rev_counts{$locs}/$ot_reads*100;
		$p_loc_ot = sprintf("%.2f", $p_loc_ot);
		$samp_otp{$name} = $overall_otp;
		$samp_raw{$name} = $raw_reads;
		$samp_ot{$name} = $ot_reads;
		$loc_otp{$locs} = $p_loc_ot;
		$loc_praw{$locs} = $p_loc_raw;
		if (exists $str_loc_otp{$locs}) {
			$str_loc_otp{$locs} = "$str_loc_otp{$locs},$p_loc_ot";
			$str_loc_praw{$locs} = "$str_loc_praw{$locs},$p_loc_raw";
			}
		else {
			$str_loc_otp{$locs} = $p_loc_ot;
			$str_loc_praw{$locs} = $p_loc_raw;
			}
		if (exists $sequences{$locs}) {
			my @seqs = split ",", $sequences{$locs};
			my @uniqs = uniq @seqs;
			my %useq_cts = ();
			my $total_seqs = @seqs;
			foreach my $useqs (@uniqs) {
				foreach my $all (@seqs) {
					if ($useqs eq $all) {
						if (exists $useq_cts{$useqs}) {
							$useq_cts{$useqs}++;
							}
						else {
							$useq_cts{$useqs} = 1;
							}
						}
					}
				}
			foreach my $sqs (sort keys %useq_cts) {
				my $p_uniq = $useq_cts{$sqs}/$total_seqs * 100;
				$p_uniq = sprintf("%.2f", $p_uniq);
				if ($p_uniq > 10) {
					$psv_seen++;
					if (exists $sig_seqs{$locs}) {
						$sig_seqs{$locs} = "$sig_seqs{$locs},$sqs";
						}
					else {
						$sig_seqs{$locs} = $sqs;
						}
					#print "$name\t$samp_raw{$name}\t$samp_ot{$name}\t$overall_otp\%\t$locs\t$p_loc_raw\%\t$p_loc_ot\%\t$p_otp\%\t$p_uniq\%\t$sqs\n";
					}
				}
			}
		else {
			#print "$name\t$samp_raw{$name}\t$samp_ot{$name}\t$overall_otp\%\t$locs\t$p_loc_raw\%\t$p_loc_ot\%\t$p_otp\%\t0\%\tNo Detected Sequences with these primers\n";
			}
		if ($psv_seen > 2) {
			$loc_psv{$locs}++;
			}
		}
	}

foreach my $targs (sort keys %fwd_primer) {
	my $unq_cts = 0;
	my $uqseqs = 0;
	my @praw = split ",", $str_loc_praw{$targs};
	my @potp = split ",", $str_loc_otp{$targs};
	my @pr_otp = split ",", $str_loc_potp{$targs};
	my $praw_total = 0;
	my $potp_total = 0;
	my $praw_ct = @praw;
	my $potp_ct = @potp;
	my $pr_otp_total = 0;
	my $pr_otp_ct = @pr_otp;
	foreach my $vals (@praw) {
		$praw_total = $praw_total + $vals;
		}
	foreach my $vals2 (@potp) {
		$potp_total = $potp_total + $vals2;
		}
	foreach my $vals3 (@pr_otp) {
		$pr_otp_total = $pr_otp_total + $vals3;
		}
	#print "$targs\t$potp_total\t$potp_ct\n";
	my $ave_praw = $praw_total/$praw_ct;
	my $ave_potp = $potp_total/$potp_ct;
	my $ave_pr_otp = $pr_otp_total/$pr_otp_ct;
	$ave_praw = sprintf("%.2f", $ave_praw);
	$ave_potp = sprintf("%.2f", $ave_potp);
	$ave_pr_otp = sprintf("%.2f", $ave_pr_otp);
	my @unq_sqs = ();
	if (exists $sig_seqs{$targs}) {
		my @the_seqs = split ",", $sig_seqs{$targs};
		@unq_sqs = uniq @the_seqs;
		$uqseqs = @unq_sqs;
		}
	#print "$targs\t$ave_praw\t$ave_potp\t$loc_psv{$targs}\t$ave_pr_otp\%\t$uqseqs\t@unq_sqs\n";  # use this print command to summarize loci and print unique seqs
	#print "$targs\t$ave_praw\t$ave_potp\t$loc_psv{$targs}\t$ave_pr_otp\%\t$uqseqs\n";  # use this print command to summarize loci
	}

foreach my $inds (sort keys %samp_otp) {
	print "$inds\t$samp_raw{$inds}\t$samp_ot{$inds}\t$samp_otp{$inds}\n";  # use this print command to summarize sample performance
	}


	

