#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  :
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

## IMPORTANT NOTE: YOU NEED TO PREPARE THE REFERENCE GENOME FILE TO RUN THIS SCRIPT
GENOME=/path/to/hg38.fa
for C in 0.7 0.8
do
	perl get.CT.rich.pl $GENOME 500 $C > hg38.bin500.CT.$C.bed
	bedtools merge -i hg38.bin500.CT.$C.bed > hg38.bin500.CT.$C.bed
done

