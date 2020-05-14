#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : Apr 24, 2020
#

set -o nounset
set -o errexit

## IMPORTANT NOTE: YOU NEED TO ADD Msuite, Bismark AND BWA-meth INTO YOUR PATH,
## BUILD THE INDICES FOR ALL OF THEM,
## AND UPDATE THE PATH TO INDICES FOR Bismark AND BWA-meth BEFORE RUNNING THIS PROGRAM
GenomeFolder=/path/to/human/genome/
BismarkIndex=/path/to/bismark/index
BWAmethIndex=/path/to/BWA-meth/index

## read length
size=${1:-36}

thread=8
## note that bismark always used two times the given thread parameter
bismark_thread=4

for CG in {0..100..10}
do
	for N in {0..9}
	do
		./Sherman -o simu.$CG.$N -l $size -I 200 -X 400 --genome_folder $GenomeFolder -pe -e 0.1 -CG $CG -CH 0.5 -n 1000000 &
	done
	wait

	for N in {0..9}
	do
		Msuite -x hg38 -m TAPS -4 \
			-1 simu.$CG.$N.read1.fq -2 simu.$CG.$N.read2.fq \
			-c $size -o Msuite.m4.$CG.$N -p $thread
		cd Msuite.m4.$CG.$N
		touch analysis.start
		make
		cd ..

		Msuite -x hg38 -m TAPS -3 \
			-1 simu.$CG.$N.read1.fq -2 simu.$CG.$N.read2.fq \
			-c $size -o Msuite.m3.$CG.$N -p $thread
		cd Msuite.m3.$CG.$N
		touch analysis.start
		make
		cd ..

		## note that bismark always used two times the given thread parameter
		mkdir -p Bismark.$CG.$N
		cd Bismark.$CG.$N
		touch analysis.start
		bismark $BismarkIndex --parallel $bismark_thread -o . \
			-1 ../simu.$CG.$N.read1.fq -2 ../simu.$CG.$N.read2.fq
		touch analysis.finished
		cd ../

		mkdir -p BWA-meth.$CG.$N
		cd BWA-meth.$CG.$N
		touch analysis.start
		bwameth.py --reference $BWAmethIndex -t $thread \
			../simu.$CG.$N.read1.fq ../simu.$CG.$N.read2.fq > simu.$CG.$N.sam 
		touch analysis.finished
		cd ../
	done
done

