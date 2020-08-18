#!/usr/bin/perl
#
# Author: Kun Sun (sunkun@szbl.ac.cn)
# This program is part of Msuite.
# Date: Apr 2020
#

use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);
use lib $Bin;
use MsuiteVersion qw($version $ver);

if( $#ARGV < 0 ) {
	print STDERR "Usage: $0 <data.dir> [alignonly=0|1]\n";
	exit 2;
}

my $dir = $ARGV[0];
my @color = ( '#E8DAEF', '#D6EAF8' );

my $alignonly = $ARGV[1] || 0;

open OUT, ">$dir/Msuite.report/index.html" or die("$!");
select OUT;

## load configuration file
## this file should be prepared by Msuite program while NOT the user
my $absPath=`readlink -m $0`;
chomp( $absPath );
my $Msuite = dirname($absPath);
$Msuite =~ s/bin\/?$/Msuite/;

open CONF, "$dir/Msuite.conf" or die( "$!" );
my %conf;
my $parameter = "<tr bgcolor=\"$color[0]\"><td>Msuite path</td><td>$Msuite</td></tr>\n";
my $i = 1;
while( <CONF> ){
	chomp;
	next if /^#/;
	next unless /\S/;
	my @l = split /\t/;
	$conf{ $l[0] } = $l[1];
	$l[1] =~ s/:/<br \/>/;	## process R1:R2 files
	$parameter .= "<tr bgcolor=\"$color[$i]\"><td>$l[0]</td><td>$l[1]</td></tr>\n";
	$i = 1 - $i;
}
close CONF;

my $pe   = ( $conf{"Sequencing mode"}  =~ /^P/i   ) ? 1 : 0;
my $TAPS = ( $conf{"Library protocol"} =~ /TAPS/i ) ? 1 : 0;

print <<HTMLHEADER;
<html>
<head>
<title>Msuite Analysis Report</title>
<style type="text/css">
td {
	text-align: left;
	padding-left: 10px;
}
</style>
</head>
<body>
<h1>Msuite Analysis Report</h1>

HTMLHEADER

############################################
print "<h2>Alignment statistics</h2>\n";

## load trim log
open LOG, "$dir/Msuite.trim.log" or die( "$!" );
my $total = <LOG>;
chomp( $total );
$total --;
my $line = <LOG>;	##Dropped : 0
$line =~ /(\d+)/;
my $dropped = $1;
my $trim = $total - $dropped;
close LOG;

## load aligner log
open LOG, "$dir/Msuite.merge.log" or die( "$!" );
$line = <LOG>;	## WATSON XXX
$line =~ /(\d+)/;
my $watson = $1;
$line = <LOG>;	## CRICK XXX
$line =~ /(\d+)/;
my $crick = $1;
my $aligned = $watson+$crick;

## load rmdup log
open LOG, "$dir/Msuite.rmdup.log" or die( "$!" );
<LOG>;	##All XXX
$line = <LOG>;	## Unique XXX
$line =~ /(\d+)/;
my $unique = $1;
$line = <LOG>;	## Duplicate XXX
$line =~ /(\d+)/;
my $dup = $1;
$line = <LOG>;	## Discard XXX
$line =~ /(\d+)/;
my $discard = $1;

print "<table id=\"alignStat\" width=\"75%\">\n",
		"<tr bgcolor=\"$color[0]\"><td width=\"70%\"><b>Total input reads</b></td>",
		"<td width=\"30%\"><b>", digitalize($total), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td><b>After preprocessing</b></td>",
		"<td><b>", digitalize($trim), sprintf(" (%.2f %%)", $trim/$total*100), "</b></td></tr>\n",

		"<tr bgcolor=\"$color[0]\"><td><b>Total aligned reads</b></td>",
		"<td><b>", digitalize($aligned), sprintf(" (%.2f %%)", $aligned/$trim*100), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td>&nbsp;&nbsp;Forward chain</td>",
		"<td>&nbsp;&nbsp;", digitalize($watson), sprintf(" (%.2f %%)", $watson/$trim*100), "</td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td>&nbsp;&nbsp;Reverse chain</td>",
		"<td>&nbsp;&nbsp;", digitalize($crick), sprintf(" (%.2f %%)", $crick/$trim*100), "</td></tr>\n",

		"<tr bgcolor=\"$color[1]\"><td><b>Unique reads</b></td>",
		"<td><b>", digitalize($unique), sprintf(" (%.2f %%)", $unique/$aligned*100), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td><b>PCR duplicates</b></td>",
		"<td><b>", digitalize($dup), sprintf(" (%.2f %%)", $dup/$aligned*100), "</b></td></tr>\n";

if( $discard ) {
print "<tr bgcolor=\"$color[1]\"><td>Discarded reads (may due to errors in index)</td>",
		"<td>", digitalize($discard), sprintf(" (%.2f %%)", $discard/$aligned*100), "</td></tr>\n";
}
print "</table>\n\n";

############################################
unless( $alignonly ) {
print "<h2>Methylation statistics</h2>\n";
## load meth.caller log
open LOG, "$dir/Msuite.CpG.meth.log" or die( "$!" );
my ($wC, $wT, $cC, $cT) = ( 0, 0, 0, 0 );
my ( $cntL, $conversionL ) = ( 0, 'NA' );
while( <LOG> ) {
	next if /^#/;
	chomp;
	my @l = split /\t/;	#chr Reads Total.wC Total.wT Total.cC Total.cT
	if( $l[0] ne 'chrL' ) {
		$wC += $l[2];
		$wT += $l[3];
		$cC += $l[4];
		$cT += $l[5];
	} else {	## reads mapped to the lambda genome
		$cntL = $l[1];
		my $coverage = $l[2]+$l[3]+$l[4]+$l[5];
		if( $coverage ) {
			$conversionL = sprintf( "%.2f %%", ($l[3]+$l[5])/$coverage*100 );
		}
	}
}
my ($wm, $cm, $tm);
if( $TAPS ) {
	$wm = sprintf("%.2f", $wT/($wC+$wT)*100);
	$cm = sprintf("%.2f", $cT/($cC+$cT)*100);
	$tm = sprintf("%.2f", ($wT+$cT)/($wC+$wT+$cC+$cT)*100);
} else {
	$wm = sprintf("%.2f", $wC/($wC+$wT)*100);
	$cm = sprintf("%.2f", $cC/($cC+$cT)*100);
	$tm = sprintf("%.2f", ($wC+$cC)/($wC+$wT+$cC+$cT)*100);
}

print "<table id=\"methStat\" width=\"75%\">\n",
		"<tr bgcolor=\"$color[0]\"><td width=\"70%\"><b>Overall methylation density</b></td>" ,
			"<td width=\"30%\"><b>$tm %</b></td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td>&nbsp;&nbsp;Forward chain</td><td>&nbsp;&nbsp;$wm %</td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td>&nbsp;&nbsp;Reverse chain</td><td>&nbsp;&nbsp;$cm %</td></tr>\n",
		"<tr bgcolor=\"$color[1]\"><td><b>Reads mapped to Lambda genome</b></td>",
			"<td><b>", digitalize($cntL), "</b></td></tr>\n",
		"<tr bgcolor=\"$color[0]\"><td><b>C-&gt;T conversion rate for lambda genome</b></td>",
			"<td><b>$conversionL</b></td></tr>\n",
	  "</table>\n\n";
}

###################################################
print "<h2>Base composition in the sequenced reads</h2>\n<ul>\n",
		"<li>Read 1</li><br /><img src=\"R1.fqstat.png\" alt=\"Base composition in read 1\"><br />\n";
print "<li>Read 2</li><br /><img src=\"R2.fqstat.png\" alt=\"Base composition in read 2\"><br />\n" if $pe;
print "</ul>\n\n";

if( $pe ) {
	print "<h2>Fragment size distribution</h2>\n",
			"<img src=\"Msuite.rmdup.size.dist.png\" alt=\"fragment size distribution\"><br />\n\n";
}

unless( $alignonly ) {
print "<h2>Methylation level per chromosome</h2>\n",
		"<img src=\"DNAm.per.chr.png\" alt=\"DNAm.per.chr\"><br />\n\n";

print "<h2>Methylation level around TSS</h2>\n",
		"<img src=\"DNAm.around.TSS.png\" alt=\"Methylation level around TSS\"><br />\n\n";

print "<h2>M-bias plot</h2>\n<ul>\n",
		"<li>Read 1</li><br /><img src=\"Msuite.R1.mbias.png\" alt=\"M-bias in read 1\"><br />\n";
print "<li>Read 2</li><br /><img src=\"Msuite.R2.mbias.png\" alt=\"M-bias in read 2\"><br />\n" if $pe;
print "</ul>\n\n";	
}

print "<h2>Analysis Parameters</h2>\n",
		"<table id=\"para\" width=\"80%\">\n",
		"<tr bgcolor=\"#888888\"><td width=\"30%\"><b>Option</b></td><td width=\"70%\"><b>Value</b></td></tr>\n",
		"$parameter</table>\n";

## HTML tail
my ($sec, $min, $hour, $day, $mon, $year, $weekday, $yeardate, $savinglightday) = localtime();
$year+= 1900;
my @month = qw/Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec/;
$hour = "0$hour" if $hour < 10;
$min  = "0$min"  if $min  < 10;
$sec  = "0$sec"  if $sec  < 10;
$mon = $month[$mon];
my $time = "$hour:$min:$sec, $mon-$day-$year";

my $url = 'https://github.com/hellosunking/Msuite/';
print "<HR align=\"left\" width=\"80%\"/>\n",
		"<h4>Generated by <a href=\"$url\" target=\"_blank\">Msuite</a> (version $ver) on $time.</h4>\n",
		"</body>\n</html>\n";

close OUT;

###############################################################
sub digitalize {
	my $v = shift;

	while($v =~ s/(\d)(\d{3})((:?,\d\d\d)*)$/$1,$2$3/){};
	return $v;

	## using bash
	## $v=`printf "%'d\n" $v`
	## chomp $v;
}

