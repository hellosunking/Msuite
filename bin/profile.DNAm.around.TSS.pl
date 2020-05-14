#!/usr/bin/perl
#
# Author: Kun Sun (sunkun@szbl.ac.cn)
# This program is part of Msuite.
# Date: Dec 2019
#

use strict;
use warnings;

if( $#ARGV < 2 ) {
	print STDERR "\nUsage: $0 <Msuite.meth.call> <tss.ext.bed> <TAPS|BS>\n\n";
	exit 2;
}

my (%C, %T);
open IN, "$ARGV[0]" or die( "$!" );
<IN>;	## skip header
while( <IN> ) {
	chomp;
	my @l = split /\t/;	##chr	Locus	Total	wC	wT	wOther	Context	cC	cT	cOther
	$C{$l[0]}->{$l[1]} += $l[3] + $l[7];
	$T{$l[0]}->{$l[1]} += $l[4] + $l[8];
}
close IN;

my (%wC, %wT, %cC, %cT);
open IN, "$ARGV[1]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/; ##chr7	137613934	137614134	CREB3L2	3800	-
	next unless exists $C{$l[0]};
	my ($c, $t) = (0, 0);
	my $chrC = $C{$l[0]};
	my $chrT = $T{$l[0]};

	for(my $i=$l[1]+1; $i<=$l[2]; ++$i) {
		if( exists $chrC->{$i} ) {
			$c += $chrC->{$i};
			$t += $chrT->{$i};
		}
	}

	if( $l[5] eq '+' ) {
		$wC{$l[4]} += $c;
		$wT{$l[4]} += $t;
	} else {
		$cC{$l[4]} += $c;
		$cT{$l[4]} += $t;
	}
}
close IN;

print "Distance\tWatson\tCrick\n";
foreach my $i ( sort {$a<=>$b} keys %wC ) {
	my ($mW, $mC) = ('NA', 'NA');
	if( $ARGV[2] =~ /TAPS/i ) {
		$mW = $wT{$i}/($wC{$i}+$wT{$i})*100 if $wC{$i}+$wT{$i}!= 0;
		$mC = $cT{$i}/($cC{$i}+$cT{$i})*100 if $cC{$i}+$cT{$i}!= 0;
	} else {
		$mW = $wC{$i}/($wC{$i}+$wT{$i})*100 if $wC{$i}+$wT{$i}!= 0;
		$mC = $cC{$i}/($cC{$i}+$cT{$i})*100 if $cC{$i}+$cT{$i}!= 0;
	}
	print join("\t", $i, $mW, $mC), "\n";
}


