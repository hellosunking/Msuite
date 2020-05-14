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
my ( @start, @finish );
@start  = stat( "$ARGV[0]/analysis.start" );
if( -e "$ARGV[0]/analysis.finished" ) {
	@finish = stat( "$ARGV[0]/analysis.finished" );
} else {
	@finish = stat( "$ARGV[0]/Msuite.merge.log" );
}

my $runtime = $finish[9] - $start[9];
print join("\t", $ARGV[0], $runtime), "\n";

