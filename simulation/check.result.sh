#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : Apr 24, 2020
#

set -o nounset
set -o errexit

seqMode=${1:-PE}
size=${2:-100}

for CG in {0..100..10}
do
	echo $CG
	for N in {0..9}
	do
		## running time
		perl extract.Msuite.stat.pl Msuite.m3.$CG.$N >> Msuite.m3.$size.time.stat
		perl extract.Msuite.stat.pl Msuite.m4.$CG.$N >> Msuite.m4.$size.time.stat
		perl extract.running.time.pl Bismark.$CG.$N  >> Bismark.$size.time.stat
		perl extract.running.time.pl BWA-meth.$CG.$N >> BWA-meth.$size.time.stat

		## alignment accuracy
		perl check.accurcy.$seqMode.pl Msuite.m3.$CG.$N/Msuite.merged.sam >> Msuite.m3.$size.acc.stat &
		perl check.accurcy.$seqMode.pl Msuite.m4.$CG.$N/Msuite.merged.sam >> Msuite.m4.$size.acc.stat &
		perl check.accurcy.$seqMode.pl Bismark.$CG.$N/simu.*.bam >> Bismark.$size.acc.stat &
		perl check.accurcy.$seqMode.pl BWA-meth.$CG.$N/simu.*.sam 1 >> BWA-meth.$size.acc.stat &
		wait
	done
done

