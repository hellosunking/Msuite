#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <in.PE.sam> [discard.Q0=0|1]\n\n";
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
while( my $r1 = <IN> ) {
	next if $r1 =~ /^@/;	## sam header
	my $r2 = <IN>;
	my @a = split /\t/, $r1;
	my @b = split /\t/, $r2;
	##18_chr15:56313259-56313575_R1	83	chr15	56313476	20	100M	*	0	0	GCCTTGCATCCCAGGGCTGAAGCCCACTTGAGCATGGTGGATAAGCTTTTTGATGTGCTGCTGGATTCGTTTTGCCAGTATTTTATTGAGGATTTTTGCA	5567

	next if $discardQ0 && $a[4]==0;

	++ $total;
	++ $Q30 if $a[4] >= 30;

	my @preset = split /[_:-]/, $a[0];
	my ($chr, $start, $end) = ( $preset[1], $preset[2], $preset[3]);

	next if $a[2] ne $chr;
	my ( $left, $right );
	if( $a[3] <= $b[3] ) {
		$left = $a[3];
		$right = $b[3] + length($b[9]) -1;
	} else {
		$left = $b[3];
		$right = $a[3] + length($a[9]) -1;
	}
	next if $left != $start || $right != $end;

	++ $correct;
}
close IN;

print join("\t", $ARGV[0], $total, $Q30, $correct, $correct/$total*100), "\n";


