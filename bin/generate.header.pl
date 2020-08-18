#!/usr/bin/perl
#
# Author: Kun Sun (sunkun@szbl.ac.cn)
# This program is part of Msuite.
# Date: Aug 2020
#

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);
use lib $Bin;
use MsuiteVersion qw($version $ver);

if( $#ARGV < 4 ) {
	print STDERR "\nUsage: $0 <in.chr.info> <index> <protocol> <mode> <read1.fq> [read2.fq]\n\n";
	exit 2;
}

my $Msuite = dirname($0);
$Msuite =~ s/bin\/?$/Msuite/;

print "\@HD\tVN:1.0\tSO:coordinate\n";
my %info;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> ) {
	chomp;
	my @l = split /\t/;
	print "\@SQ\tSN:$l[0]\tLN:$l[1]\n";
}
close IN;

print "\@PG\tID:Msuite\tPN:$Msuite\tVN:$ver\tDS:genome=$ARGV[1];protocol=$ARGV[2];mode=$ARGV[3];";
if( $#ARGV > 4 ) {
	print "reads=$ARGV[4]:$ARGV[5]\n";
} else {
	print "reads=$ARGV[4]\n";
}

