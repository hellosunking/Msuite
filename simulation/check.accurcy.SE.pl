#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.SE.sam> [discard.Q0=0|1]\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

my $discardQ0 = $ARGV[1] || 0;

my $total = 0;
my $Q30 = 0;
my $correct = 0;

if( $ARGV[0] =~ /bam$/ ) {
	open IN, "samtools view $ARGV[0] |" or die( "$!" );
} else {
	open IN, "$ARGV[0]" or die( "$!" );
}
while( <IN> ) {
	next if /^@/;	## sam header
	my @l = split /\t/;
	##18_chr15:56313259-56313575_R1	83	chr15	56313476	20	100M	*	0	0	GCCTTGCATCCCAGGGCTGAAGCCCACTTGAGCATGGTGGATAAGCTTTTTGATGTGCTGCTGGATTCGTTTTGCCAGTATTTTATTGAGGATTTTTGCA	5567

	next if $discardQ0 && $l[4]==0;

	++ $total;
	++ $Q30 if $l[4] >= 30;

	my @preset = split /[_:-]/, $l[0];
	my ($chr, $start, $end) = ( $preset[1], $preset[2], $preset[3]);

	next if $l[2] ne $chr;
	if( $l[1] & 0x10 ) {	## mapped to crick chain
		++ $correct if $end == $l[3] + length($l[9]) -1;
	} else {
		++ $correct if $start == $l[3];
	}
}
close IN;

print join("\t", $ARGV[0], $total, $Q30, $correct, $correct/$total*100), "\n";


