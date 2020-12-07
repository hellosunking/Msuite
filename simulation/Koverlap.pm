#!/usr/bin/perl

package Koverlap;
#by ahfyth
#2020/4/7

use strict;
use warnings;

#function: buildmap
#
#parameters: <file_name> [ idpos chrpos startpos endpos ]
#
#default: line in file : id chr start end
#or you may give me some information on the position of each field
#
#return    : \%map
#%map      : chr -> \@genelist
# @genelist: \@gene1, \@gene2 ... 
#	@gene  : geneid start end
#				they are sorted by position on chromsome
#

sub buildmap {
	my $filename = shift;
	my $idpos    = shift;

	my $chrpos   = shift;
	my $startpos = shift;
	my $endpos   = shift;

	open FILE, "$filename" or die("open file $filename failed!\n");
	print STDERR "Use file: $filename\n";
#	print STDERR "Use field postion : id=>$idpos chr=>$chrpos start=>$startpos end=>$endpos\n\n";

	my %map = ();
	my @line;
	my ($id, $chr, $start, $end);
	if( defined( $endpos ) ) {
		while( <FILE> ) {
			chomp;
			@line	= split /\t/;
			$id		= (defined $line[$idpos]) ? $line[$idpos] : "MAP_$.";
			$chr	= $line[ $chrpos];
			$start	= $line[ $startpos ];
			$end	= $line[ $endpos ];
#			print "load $_";
#			print "Add gene : id=>$id chr=>$chr start=>$start end=>$end\n";
			my @gene= ( $id, $start, $end );
			if( defined( $map{ $chr } ) ) {
#				print "find a defined for $chr\n";
				my $genelist = $map{ $chr };
				push( @$genelist, \@gene );
#				print "find a defined for $chr : ", $#$genelist+1, "\n";
			} else {
#				print "build genelist for $chr\n";
				my @genelist = ( \@gene );
				$map{ $chr } = \@genelist;
			}
		}
	} else {
#		view the file as a bed file
		while( <FILE> ) {
			chomp;
			@line	= split /\t/;
			$id		= (defined $line[$idpos]) ? $line[$idpos] : "MAP_$.";
			$chr	= $line[ 0 ];
			$start	= $line[ 1 ];
			$end	= $line[ 2 ];
			my @gene= ( $id, $start, $end );
			if( defined( $map{ $chr } ) ) {
				my $genelist = $map{ $chr };
				push( @$genelist, \@gene );
			} else {
				my @genelist = ( \@gene );
				$map{ $chr } = \@genelist;
			}
		}
	}
	print STDERR "Done: $. records loaded.\n";
	close FILE;

#	printmap( \%map );
	foreach $chr ( keys %map ) {
#		print STDERR "Sorting $chr\n";
		my $genelist = $map{ $chr };
		my @sorted = sort { $$a[1] <=> $$b[1] || $$a[2] <=> $$b[2] } @$genelist;
		$map{ $chr } = \@sorted;
	}
	return \%map;
}


#
#function : overlap
#
#parameters : \@query, \%map
#
#@query    : chr start end [OTHER]
#%map      : chr -> \@genelist
# @genelist: \@gene1, \@gene2 ... 
#	@gene  : geneid start end
#				they are sorted by position on chromsome
#
#return    : \@hit
#@hit      : [ \@match1 [ \@match2 ... ] ]
#  @match  : query_chr query_start query_end geneid genestart geneend \@otherqueryinfo
#
#
#this function is to find the genes that overlaps with the query bed style fragment
#
#use binarySearch to get the match genes
#

sub overlap {
	my $query    = shift;
	my $map      = shift;

#	print "Query : ", join("\t", @$query), "\n";
	my @hit = ();

	my $chr = $$query[0];
	return \@hit if( ! defined( $$map{ $chr } ) );

	my $genelist = $$map{ $chr };

	my $qstart = $$query[1];
	my $qend   = $$query[2];

	my $gstart = 0;
	my $gend   = $#$genelist;

	while( $gend - $gstart > 1 ) {
		my $gcenter = ( $gstart + $gend ) >> 1;
		my $detail = $$genelist[$gcenter];	# geneid start end
		if( $$detail[1] > $qend ){ $gend = $gcenter - 1; }	# to find the last gene that starts before fragment end
		#elsif( $$detail[1] < $qend ){ $gstart = $gcenter + 1; }
		else { $gstart = $gcenter; }
	}

	$gend ++;
	$gend = $#$genelist if $gend > $#$genelist;
	while( $gend >= 0 ) {
		my $detail = $$genelist[$gend];
		last if( $$detail[2] < $qstart );

		if( $$detail[1] <= $qend ) {
#			print STDERR "Find $gend : $$detail[1]-$$detail[2]\n";
			my @extra;
			for( my $i=3; $i<=$#$query; ++$i ) {
				push @extra, $$query[$i];
			}
			my @match = ( $chr, $qstart, $qend, $$detail[0], $$detail[1], $$detail[2], \@extra );
			push( @hit, \@match );
		}
		-- $gend;
	}

	return \@hit;
}


#
#function : printmap
#parameters : \%map
#
#

sub printmap {
	my $map = shift;

	foreach my $chr ( keys %$map ) {
		print "chr $chr:\n";
		my $genelist = $$map{ $chr };
		foreach my $gene ( @$genelist ) {
			print join("\t", @$gene), "\n";
		}

		print "\n";
	}
}


sub gethithead {
	return "chr\tfstart\tfend\tgid\tgstart\tgend\textra";
}

#
#function : printhit
#parameters : \@hit
#
#

sub printhit {
	my $hit = shift;

	for( my $i=0; $i<=$#$hit; ++$i ) {
		my $here = $$hit[$i];
		for( my $j=0; $j<$#$here; ++$j ) {
			print "$$here[$j]\t";
		}
		my $extra = $$here[-1];
		print join(",", @$extra ), "\n";
	}
	return $#$hit+1 ;
}

1;

