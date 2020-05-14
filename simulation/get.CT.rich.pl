#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;
#use KSLIB::loadGenome qw/loadGenome/;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <genome.fa> <bin=1k> <CT/GA_threshold=0.7>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $bin = $ARGV[1] || 1000;
my $threshold = $ARGV[2] || 0.7;

## load genome
open IN, "$ARGV[0]" or die( "$!" );
my %g;
my $chr='';
while( <IN> ) {
	chomp;
	if( /^>(chr\S+)/ ) {
		$chr = $1;
		$g{$chr} = '';
	} else {
		$g{$chr} .= uc $_;
	}
}
close IN;

foreach $chr ( keys %g ) {
	my $len = length( $g{$chr} );
#	print STDERR "Loading $chr: $len\n";
	for( my $i=0; $i<$len; $i+=$bin ){
		my $seq = substr( $g{$chr}, $i, $bin );
		my $extracted = length( $seq );
		next if $extracted < $bin*0.9;	## usually the end of a chromosome
		my $Ns = ($seq=~s/N//g);
		if( $Ns > 0.1 * $extracted ) {	## too many Ns
#			print STDERR "Skip Bin $i due to $Ns Ns.\n";
			next;
		}
		$extracted = length( $seq );
		my $CT = ($seq=~s/[CT]//g);
		my $ratio = $CT/$extracted;
		if( $ratio>=$threshold || $ratio<=1-$threshold ) {
			print join("\t", $chr, $i, $i+$bin), "\n";
		}
#		print STDERR "Bin $i: CT=$CT\n";
	}
}

