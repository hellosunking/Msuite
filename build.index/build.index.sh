#!/bin/bash

#
# Author: Kun Sun @ SZBL (hellosunking@foxmail.com)
#

set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

if [ $# -lt 3 ]
then
	echo -e "\nUsage: $0 <genome folder|fa> <refSeq.txt> <index.identifier> [thread=auto]\n" >/dev/stderr
	echo -e "This is a utility program of Msuite, designed to build genome references." >/dev/stderr
	echo -e "Please refer to README file for more information.\n" >/dev/stderr
	exit 2
fi

PRG=`dirname $0`
indexDIR=$PRG/../index
lambdaGenome=$PRG/lambda.genome.fa

id=$3

if [ $# -gt 3 ]
then
	thread=$4
else
	thread=`grep processor /proc/cpuinfo | wc -l`
fi
echo "=> INFO: Threads used: $thread"

bb=`which bowtie2-build 2>/dev/null`
if [ -z "$bb" ]; then
	echo -e "\e[31mFatal error: Could not find 'bowtie2-build' in your path!\e[39m" >/dev/stderr
	exit 1
fi

ver=`$bb --version | perl -ne 'print $1 if /bowtie2-build version ([\d\.]+)/'`
if [ -z "$ver" ]; then
	ver=unknown
fi

echo -e "\n\e[32mINFO: bowtie2-build (version $ver) found at '$bb'.\e[39m\n"
echo -e "\e[34mInput genome:      $1"
echo -e "RefSeq annotation: $2"
echo -e "Lambda genome:     $lambdaGenome"
echo -e "Index identifier:  $3"
echo -e "Index location:    $indexDIR/$3\e[39m\n"

echo "Preprocessing genome ..."
## prepare directories
mkdir -p $indexDIR/$id
mkdir -p $indexDIR/$id/Mode3 $indexDIR/$id/Mode4
perl $PRG/process.genome.pl $1 $indexDIR/$id/ $lambdaGenome

echo "Building 4-letter indices ..."
$bb --threads $thread $indexDIR/$id/CG2TG.fa $indexDIR/$id/Mode4/CG2TG >/dev/null
$bb --threads $thread $indexDIR/$id/CG2CA.fa $indexDIR/$id/Mode4/CG2CA >/dev/null

echo "Building 3-letter indices ..."
$bb --threads $thread $indexDIR/$id/C2T.fa   $indexDIR/$id/Mode3/C2T >/dev/null
$bb --threads $thread $indexDIR/$id/G2A.fa   $indexDIR/$id/Mode3/G2A >/dev/null

rm -f $indexDIR/$id/CG2TG.fa $indexDIR/$id/CG2CA.fa $indexDIR/$id/C2T.fa $indexDIR/$id/G2A.fa

echo "Processing refSeq gene annotation ..."
perl $PRG/process.refGene.pl $2 >$indexDIR/$id/tss.ext.bed

echo
echo -e "\e[35mDone. Now you can use \"-x $3\" to align your data to this genome in Msuite.\e[39m"
echo

