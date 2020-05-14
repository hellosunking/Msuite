#!/usr/bin/perl
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :

use strict;
use warnings;

if( $#ARGV < 0 ) {
	print STDERR "\nUsage: $0 <dir>\n\n";
#	print STDERR "\nThis program is designed to \n\n";
	exit 2;
}

## running time
my @start  = stat( "$ARGV[0]/analysis.start" );
my @finish = stat( "$ARGV[0]/Msuite.merge.log" );

my $runtime = $finish[9] - $start[9];

## mapping infor
open IN, "$ARGV[0]/Msuite.rmdup.log" or die( "$!" );
my $line = <IN>;
close IN;
$line =~ /(\d+)/;
my $aligned = $1;

## methy call
open LOG, "$ARGV[0]/Msuite.meth.log" or die( "$!" );
my ($wC, $wT, $cC, $cT) = ( 0, 0, 0, 0 );
while( <LOG> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;     #chr Reads Total.wC Total.wT Total.cC Total.cT
	if( $l[0] ne 'chrL' ) {
		$wC += $l[2];
		$wT += $l[3];
		$cC += $l[4];
		$cT += $l[5];
	}
}
close LOG;
my $meth = sprintf("%.4f", ($wT+$cT)/($wC+$wT+$cC+$cT)*100);

print join("\t", $ARGV[0], $runtime, $aligned, $meth), "\n";

