#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : Apr 24, 2020
#

set -o nounset
set -o errexit

PRGHERE=`readlink -m $0`
PRG=`dirname $PRGHERE`
indexDIR=$PRG/index

RED="\033[0;31m"
BLUE="\033[0;34m"
BG="\033[0m"

echo "${BLUE}This program is designed to run Msuite on the testing dataset (testing_dataset/simu*.fq)."
echo "This testing dataset is generated by SHERMAN and only those in the CT/GA-rich regions are kept."
echo "Fore more details of this dataset, please refer to README file in the package."
echo

## determine thread
totalThread=`grep processor /proc/cpuinfo | wc -l`
if [ $totalThread -gt 8 ]
then
	let thread=$totalThread-4
else
	thread=4
fi
echo -e "INFO: Running threads used: $thread.${BG}"

## check hg38 index
if [ -s $indexDIR/hg38/tss.ext.bed ]
then
	echo -e "${BLUE}INFO: It seems that the index for hg38 genome has been built before, I will skip this step.${BG}"
else
	echo -e "${BLUE}INFO: The index for hg38 genome is not available, I will build it first."
	echo -e "This step may need some time, please be patient ...${BG}"
	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -O $PRG/index/hg38.fa.gz
	$PRG/build.index/build.index.sh $PRG/index/hg38.fa.gz $PRG/build.index/hg38.refGene.txt.gz hg38 $thread

	if [ $? != 0 ]
	then
		echo -e "${RED}ERROR: Build index failed! Please check whether you have correctly installed 'bowtie2'!${BG}"
		exit 1
	fi
fi
echo

## Uncompress FASTQ files
echo Decompressing read files ...
bzip2 -cd $PRG/testing_dataset/simu.read1.fq.bz2 > $PRG/testing_dataset/simu.read1.fq &
bzip2 -cd $PRG/testing_dataset/simu.read2.fq.bz2 > $PRG/testing_dataset/simu.read2.fq &
wait

########### 4-letter mode ##########
## build makefile
echo -e "${BLUE}Generating makefile for 4-letter alignment ...."
CMD="$PRG/msuite -x hg38 -1 $PRG/testing_dataset/simu.read1.fq -2 $PRG/testing_dataset/simu.read2.fq -4 -m TAPS -c 36 -q 0 -p $thread -o $PRG/testing_dataset/Msuite.Mode4"
echo -e "Command line: $CMD${BG}"
$CMD
if [ $? != 0 ]
then
	echo -e "${RED}ERROR: Generate makefile failed! Please check whether you have correctly installed 'bowtie2' and 'samtools'!${BG}"
	exit 1
fi
## run analysis
echo -e "${BLUE}========== Performing analysis (4-letter mode) ==========${BG}"
echo
cd $PRG/testing_dataset/Msuite.Mode4
#time make
if [ $? != 0 ]
then
	echo -e "${RED}ERROR: analysis failed! Please check the run logs for more information.${BG}"
	exit 1
fi
cd ../../
echo -e "${BLUE}========== Analysis finished (4-letter mode) ==========${BG}"
echo
echo

########### 3-letter mode ##########
echo -e "${BLUE}Generating makefile for 3-letter alignment ...."
## build makefile
CMD="$PRG/msuite -x hg38 -1 $PRG/testing_dataset/simu.read1.fq -2 $PRG/testing_dataset/simu.read2.fq -3 -m TAPS -c 36 -q 0 -p $thread -o $PRG/testing_dataset/Msuite.Mode3"
echo -e "Command line: $CMD${BG}"
$CMD
if [ $? != 0 ]
then
	echo -e "${RED}ERROR: Generate makefile failed! Please check whether you have correctly installed 'bowtie2' and 'samtools'!${BG}"
	exit 1
fi
## run analysis
echo -e "${BLUE}========== Performing analysis (3-letter mode) ==========${BG}"
echo
cd $PRG/testing_dataset/Msuite.Mode3
time make
if [ $? != 0 ]
then
	echo -e "${RED}ERROR: analysis failed! Please check the run logs for more information.${BG}"
	exit 1
fi
cd ../../
echo -e "${BLUE}========== Analysis finished (3-letter mode) ==========${BG}"
echo
echo
echo

## prompt result
cd $PRG
echo -e "${BLUE}Done. The results for 4-letter mode are written in '$PRG/testing_dataset/Msuite.Mode4' and the HTML report is stored in '$PRG/testing_dataset/Msuite.Mode4/Msuite.report/index.html'; the results for 3-letter mode are written in '$PRG/testing_dataset/Msuite.Mode3' and the HTML report is stored in '$PRG/testing_dataset/Msuite.Mode3/Msuite.report/index.html'.${BG}"
echo
echo "For more information, please refer to README file in the package."
echo
