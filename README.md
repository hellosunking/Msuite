
# Msuite: Multi-mode DNA methylation data analysis suite
Version 1.0.3, Apr 2020<br />
Authors: Kun Sun, Lishi Li, Li Ma, Yu Zhao, Lin Deng, Huating Wang and Hao Sun<br />
Software implemented by Kun Sun \(sunkun@szbl.ac.cn\)<br />
<br />
Distributed under the
[GNU General Public License v3.0 \(GPLv3\)](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3")
for personal and academic usage only.<br />
For detailed information please read the license file under `license` directory.

---

## Installation
`Msuite` is written in `Perl` and `R` for Linux/Unix platform. To run `Msuite` you need a Linux/Unix
machine with `Bash 4 (or higher)`, `Perl 5.10 (or higher)` and `R 2.10 (or higher)` installed.

This source package contains pre-compiled executable files using `g++ v4.8.5` for Linux x86_64 system.
If you could not run the analysis normally (which is usually caused by low version of `libc++` library),
or you want to build a different version optimized for your system, you can re-compile the programs:
```
user@linux$ make clean && make
```

Note that `Msuite` depends on the following software:

* [bowtie2](https://github.com/BenLangmead/bowtie2 "bowtie2")
* [samtools](http://samtools.sourceforge.net/ "samtools")

Please install them properly and make sure that they are included in your `PATH`.
In addition, please make sure that the version of your `g++` compiler is higher than 4.8
(you can use `g++ -v` to check it).

Before running `Msuite`, genome indices must be built. To this end, we have prepared a utility named
`build.index.sh` under the `build.index` directory. To use it, you need to prepare the genome sequence
(either in one multi-fasta file or a directory containing the sequences for individual chromosomes) and
RefSeq annotation for your genome (we have included files for mm9/mm10/hg19/hg38 in this package).
You can download the annotations for other species/genome versions from the
[UCSC genome browser](http://genome.ucsc.edu/ "UCSC Genome Browser").
(The RefSeq annotation is used to profile the methylation level around Transcription Start Sites as a
quick quality control of your data.)

Then, you can build the genome indices using the following command:
```
user@linux$ build.index/build.index.sh GENOME.FA(or GENOME.DIR) REFSEQ.txt Genome.ID
```
Note that this utility will automatically incorporate the Lambda genome to build the genome indices.
`gzip` or `bzip2` compression of `GENOME.FA` and `REFSEQ.txt` are also supported, but the files must
contain the corresponding suffix (i.e., `REFSEQ.txt.gz` for `gzip` compressed and `REFSEQ.txt.bz2` for
`bzip2` compressed file). The `Genome.ID` is an identifier that you specified to name your genome and
the indices will be written to the `index` directory under the root of `Msuite`. You can add as many
genomes to `Msuite` as you need.

## Run Msuite
The main program is `msuite`. You can add its path to your `.bashrc` file under the `PATH` variable
to call it from anywhere, or you can run the following command to add it to your current session:
```
user@linux$ PATH=$PATH:$PWD
```
Call `msuite` without any parameters to see the usage (or use '-h' option):
```
########## Msuite: Multi-mode DNA methylation data analysis suite ##########

Author : Kun Sun (sunkun@szbl.ac.cn)
Version: 1.0.3 (Apr 2020)

Usage: Msuite [options] -x index -1/-U Read1.fq [ -2 Read2.fq ] -c cycle -o out.dir

Compulsory parameters:

  -1/-U Read1.fq   Specify the path to the file containing read 1
                   If your data is Paired-end, specify read 2 file using -2 option
                   If -U is used, then -2 will be ignored

                   If you have multiple files for your sample, please use ',' to separate them

  -x index         Specify the genome index
                   Please refer to README file on how to build index for Msuite

  -c cycle         Specify the sequencing cycle of your data

  -o out.dir       Specify the output directory
                   Note that your specified directory will be created if it does not exist
                   otherwise the files under that directory could get over-written


Optional parameters:

  -2 Read2.fq      Specify the path to the file containing read 2
                   Use this parameter if your data is generated in paired-end mode

                   If you have multiple files for your sample, please use ',' to separate them
                   and make sure that all the files are well paired in -1 and -2 options

  -3               Use 3-letter alignment (default)
  -4               Use 4-letter alignment
                   Note that the above two options are mutually exclusive

  -k kit           Specify the library preparation kit (default: illumina)
                   Note that only 'illumina', 'nextera' and 'bgi' are acceptable

  -m TAPS/BS       Specify the library protocol (default: BS)
                   Note that only 'TAPS' and 'BS' are acceptable

  -p threads       Specify how many threads should be used (default: 4)

  --phred33        Read cycle quality scores are in Phred33 format (default)
  --phred64        Read cycle quality scores are in Phred64 format
                   Note that the above two options are mutually exclusive

  -q score         The minimum quality score to keep the cycle (default: 20)
                   Note that 20 means 1% error rate, 30 means 0.1% error rate in Phred

                   Sometimes quality scores start from 35 ('#') in the FASTQ files,
                   in this case you could adjust '-q' option, e.g., '--phred33 -q 22'

  --minsize size   Minimum read size to be kept for alignment (default: 22)

  --minins MIN     Minimum insert size (default: 0)
  --maxins MAX     Maximum insert size (default: 1000)
                   Note that the above two options will be ignored for Single-End data

  --align-only     Stop after alignment (i.e., do not perform DNA methylation call and
                   visualization around TSS; default: not set)

  -h/--help        Show this help information and quit
  -v/--version     Show the software version and quit


Please refer to README file for more information.
```

**IMPORTANT NOTE**: If your data is generated using BS-seq protocol, you MUST use the 3-letter mode and set
`-m BS`. 4-letter mode ONLY supports processing of TAPS/5hmC-CATCH (or similar bisulfite-free assays) data.
In addition, Msuite could also directly analyze the data generated by ATAC-me or similar protocols via
setting `-k nextera`.


### Example 1
If your data is generated using TAPS protocol in 75 bp * 2 (paired-end) mode, and you want to align your
data to the hg38 reference genome in 4-letter mode, and you want to use 16 threads to speed-up the analysis,
then you can run:
```
user@linux$ msuite -1 /path/to/read1.fq -2 /path/to/read2.fq -c 75 -x hg38 \
                   -4 -m TAPS -p 16 -o /path/to/output/dir
```

### Example 2
If your data is generated using BS-seq protocol in 100 bp * 1 (single-end) mode, and you want to align your
data to the mm10 reference genome (note that you MUST use 3-letter mode here), and you want to use 32 threads
to speed-up the analysis, then you can run:
```
user@linux$ msuite -1 /path/to/read1.fq -c 100 -x mm10 \
                   -3 -m BS -p 32 -o /path/to/output/dir
```

### Example 3
If your data is generated using TAPS protocol in 75 bp * 2 (paired-end) mode, and you have 3 lanes of data,
you want to align your data to the hg19 reference genome in 4-letter mode, and you want to use 48 threads to
speed-up the analysis, then you can run:
```
user@linux$ msuite -1 /path/to/lane1.read1.fq,/path/to/lane2.read1.fq,/path/to/lane3.read1.fq \
                   -2 /path/to/lane1.read2.fq,/path/to/lane2.read2.fq,/path/to/lane3.read2.fq \
                   -c 75 -x hg19 -4 -m TAPS -p 48 -o /path/to/output/dir
```
`Msuite` will check the data and dependent programs then generate a `makefile` under `/path/to/output/dir`
('-o' option). Then you can go to `/path/to/output/dir` and run `make` to perform the analysis:
```
user@linux$ cd /path/to/output/dir; make
```
We have prepared a testing dataset under the `testing_dataset` directory. It contains *in silico* generated
reads following the TAPS (TET-assisted pyridine borane sequencing) protocol using
[SHERMAN](http://www.bioinformatics.babraham.ac.uk/projects/sherman/) software (key parameters: C-&gt;T
conversion rate: 50% for CpG loci, C-&gt;T conversion rate in CpH sites: 0.5%, sequencing error rate: 0.1%).
Note that the reads are restricted to CT/GA-rich regions (CT/GA proportion >= 80%) to illustrate the advantage
of 4-letter over 3-letter mode. We also have prepared a work shell to run `msuite` on this dataset using
both 4- and 3-letter modes:
```
user@linux$ ./run_testing_dataset.sh
```
Note that this script will automatically build the indices for hg38 genome if you have not done this before.

You can compare the performance of 4- and 3-letter analysis by inspecting the outputs, which will be written
into `testing_dataset/Msuite.Mode3/` and `testing_dataset/Msuite.Mode4/`.


## Outputs explanation
`Msuite` outputs all the results for the given region in the directory specified by `-o OUTDIR` option.
`Msuite` will write the analysis report into a HTML file named `Msuite.report/index.html`, which records the
essential statistics and visualizations for quality control, mappability, overall methylation level and
conversion rate (estimated using reads mapped to the Lamda genome).

The alignment results are recorded in the file `Msuite.final.bam` (in standard BAM format).
The CpG methylation calls are recorded in the file `Msuite.meth.call` and `Msuite.meth.bedgraph`.

You can run `make clean` in the OUTDIR to delete the intermediate files to save storage space.


## Utilities
### Mviewer
`Msuite` contains a visualization tool named `Mviewer`, adapted from the authors' previous
[BSviewer](http://sunlab.cpy.cuhk.edu.hk/BSviewer/) software. It is specially optimized to be compatiable
with `Msuite` alignment results and provides nucleotide-level, genotype-preserved DNA methylation data
visualization. For more information, please refer to README file in `Mviewer` directory.


### Others
`Msuite` also provides other utilities under the `util` directory. The `profile.meth.pl` program is designed
to summarize the methylome into bins. You can use it to prepare data for `Circos` plots.
`extract.meth.in.region` is an extended version, which is designed to extract the covered CpG sites, C-count
and T-count in the given regions (e.g., CpG islands).

The `pe_bam2bed.pl` and `se_bam2bed.pl` are designed to translate the aligned BAM file into BED format file,
and `bed2wig` is designed to translate BED file into WIG files (e.g., for coverage profiles).

The `simulation` directory contains all the scripts to reproduce our benchmark test results.

---
Please send bug reports to Kun Sun \(sunkun@szbl.ac.cn\).<br />
Msuite package is freely available at
[https://github.com/hellosunking/Msuite/](https://github.com/hellosunking/Msuite/ "Msuite @ Github").

