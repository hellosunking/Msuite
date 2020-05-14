#!/usr/bin/perl

#
# Author: Kun Sun @ SZBL (hellosunking@foxmail.com)
#

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <in.genome> <out.prefix> <lambda.genome>\n";
	print STDERR "\nThis program is designed to treat CpG sites.\n\n";
	exit 2;
}

my @falist;

if( -d "$ARGV[0]" ) {	## parameter 1 is a directory, then load all FASTA files under this directory
	print "INFO: The given parameter is a directory! I will load all the fasta files in it!\n";
	opendir DIR, "$ARGV[0]" or die( "$!" );
	while( my $fa = readdir(DIR) ) {
		push @falist, "$ARGV[0]/$fa" if $fa =~ /\.fa(sta)?$/;
	}
	closedir DIR;
} else {	## parameter 1 is a file, which I assume should contain all sequences
	print "INFO: The given parameter is a file and I will load it directly!\n";
	push @falist, $ARGV[0];
}

## add lambda genome
push @falist, $ARGV[2];

my %g;
my $chr;
print "Loading genome files ...\n";
foreach my $fa ( @falist ) {
#	print STDERR "Loading $fa \n";
	if( $fa =~ /\.gz$/ ) {
		open IN, "gzip -cd $fa |" or die( "$!" );
	} elsif( $fa =~ /\.bz2$/ ) {
		open IN, "bzip2 -cd $fa |" or die( "$!" );
	} else {
		open IN, "$fa" or die( "$!" );
	}

	while( <IN> ) {
		chomp;
		if( s/^>// ) {
			$chr = $_;
			$chr =~ s/\s.*$//;
			$chr = "chr$chr" unless $chr=~/^chr/;
#			print STDERR "Adding $chr\n";
			$g{$chr} = '';
		} else {
			$g{$chr} .= uc $_;
		}
	}
	close IN;
}

print "Processing genome files ...\n";
my $bp_per_len = 50;

open CG2TG, ">$ARGV[1]/CG2TG.fa" or die( "$!" );
open CG2CA, ">$ARGV[1]/CG2CA.fa" or die( "$!" );

open C2T, ">$ARGV[1]/C2T.fa" or die( "$!" );
open G2A, ">$ARGV[1]/G2A.fa" or die( "$!" );

open SIZE, ">$ARGV[1]/chr.info" or die( "$!" );
open ORI,  ">$ARGV[1]/genome.fa" or die( "$!" );

my ($cnt1, $cnt2, $cnt3, $cnt4) = ( 0, 0, 0, 0 );
foreach $chr ( sort keys %g ) {
	my $raw = $g{$chr};
	print SIZE "$chr\t", length($raw), "\n";

	## record the raw sequence and chr size information
	my $i = 0;
	my $seq;
	print ORI ">$chr\n";
	while( 1 ) {
		$seq = substr( $raw, $i, $bp_per_len );
		my $len = length($seq);
		print ORI "$seq\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}

	## CG->TG
	$cnt1 += ($raw =~ s/CG/TG/g);
	print CG2TG ">$chr\n";
	$i = 0;
	while( 1 ) {
		$seq = substr( $raw, $i, $bp_per_len );
		my $len = length($seq);
		print CG2TG "$seq\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}

	## CG->CA
	$raw = $g{$chr};
	$cnt2 += ($raw =~ s/CG/CA/g);
	print CG2CA ">$chr\n";
	$i = 0;
	while( 1 ) {
		$seq = substr( $raw, $i, $bp_per_len );
		my $len = length($seq);
		print CG2CA "$seq\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}

	## C->T
	$raw = $g{$chr};
	$cnt3 += ($raw =~ s/C/T/g);
	print C2T ">$chr\n";
	$i = 0;
	while( 1 ) {
		$seq = substr( $raw, $i, $bp_per_len );
		my $len = length($seq);
		print C2T "$seq\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}

	## G->A
	$raw = $g{$chr};
	$cnt4 += ($raw =~ s/G/A/g);
	print G2A ">$chr\n";
	$i = 0;
	while( 1 ) {
		$seq = substr( $raw, $i, $bp_per_len );
		my $len = length($seq);
		print G2A "$seq\n" if $len != 0;

		last if $len != $bp_per_len;
		$i += $bp_per_len;
	}
}

print "Done. Coversions:\n",
	  "CG->TG: $cnt1\n",
	  "CG->CA: $cnt2\n",
	  "C -> T: $cnt3\n",
	  "G -> A: $cnt4\n";

close CG2TG;
close CG2CA;
close C2T;
close G2A;
close SIZE;
close ORI;

