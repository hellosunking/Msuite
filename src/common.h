/*
 * common.h
 *
 * This header file records the constants used in Msuite
 *
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date  : Aug 2020
 *
 * */

/*
 *  Illumina adapters
 *
 *  Original adapter pair from illumina:
 *  5': ACACTCTTTCCCTACACGACGCTCTTCCGATCT
 *  3': AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
 *
 *  sequence model:
 *
 *                *read1 -->
 *  5' adapter - sequence sequence sequence - 3' adapter
 *                                <-- read2*
 *
 *  so read1 may contains 3' adapter, name it adapter_r1,
 *  read2 may contains reversed and ACGT-paired 5' adapter, name it adapter_r2
 *
 * */

#ifndef _MSUITE_COMMON_
#define _MSUITE_COMMON_

// fastq statistics structure
typedef struct {
	unsigned int A;
	unsigned int C;
	unsigned int G;
	unsigned int T;
	unsigned int N;
} fastqstat;

typedef struct {
	const char* adapter_r1;
	const char* adapter_r2;
	const unsigned int adapter_len;
	const char* adapter_index;
} adapter_info;

//illumina sequencing adapters
const char* illumina_adapter_sequence = "AGATCGGAAGAGC";
const unsigned int illumina_adapter_len = 13;	//strlen(illumina_adapter_sequence)
const char* illumina_adapter_index = "AGA";
static adapter_info illumina_adapter = {illumina_adapter_sequence, illumina_adapter_sequence,
										illumina_adapter_len, illumina_adapter_index};
//nextera sequencing adapters
const char* nextera_adapter_sequence = "CTGTCTCTTATACACATCT";
const unsigned int nextera_adapter_len = 19;	//strlen(nextera_adapter_sequence)
const char* nextera_adapter_index = "CTG";
static adapter_info nextera_adapter = {nextera_adapter_sequence, nextera_adapter_sequence,
										nextera_adapter_len, nextera_adapter_index};

//bgi sequencing adapters
const char* bgi_adapter1_sequence = "AAGTCGGAGGCCAAGCGGTC";
const char* bgi_adapter2_sequence = "AAGTCGGATCGTAGCCATGT";
const unsigned int bgi_adapter_len = 19;	//strlen(bgi_adapter_sequence)
const char* bgi_adapter_index = "AAG";
static adapter_info bgi_adapter = {bgi_adapter1_sequence, bgi_adapter2_sequence,
										bgi_adapter_len, bgi_adapter_index};

const char FILE_SEPARATOR = ',';		// separator if multiple files are provided

// seed and error configurations
const unsigned int impossible_seed = 10000;
const unsigned int MAX_READ_LENGTH = 1024;
const unsigned int MAX_CONVERSION  = 1024;

// C>T and G>A conversion related
const char LINE_NUMBER_SEPARATOR = '+';
const char CONVERSION_LOG_START  = '~';
const char CONVERSION_LOG_END    = '#';	// original version uses '&'
//const char CONVERSION_LOG_SEPARATOR = ':';	//original version uses ':'
const char CONVERSION_LOG_SEPARATOR = ';';	//':' is replaced because it's commonly used in illumina seqName
const char KEEP_QUAL_MARKER      = '|';
const char NORMAL_SEQNAME_START  = '@';

const int MAX_CHANGES_CNT  =  256;
const int MAX_SEQ_CYCLE    =  255;
const int MAX_SEQNAME_SIZE = 1024;
const int MAX_ITERM_SIZE   =   32;	// MAX sam iterm (e.g., chr, mateinfo) size
// there could be something like chr12_GA_converted, so set to 32

//const unsigned int MIN_ALIGN_SCORE = 10;	// minimum alignemnt score to keep the record
const int IMPOSSIBLE_AS_SCORE = - (1<<20);
const int DISCARD_AMB_HIT     = - (1<<15);
const int AMB_HIT_MARKER      = - (1<<10);
const int MIN_ALIGN_SCORE_AMB = 2;	// minimum alignemnt score to keep the ambigous record
const int MIN_ALIGN_SCORE_UNQ = 2;	// minimum alignemnt score to keep the unique record

const int MIN_ALIGN_SCORE = 0;		// minimum alignemnt score to keep the record (for compatible with previous versions)

// for ambigous hits, ONLY keep the one with very high score while the other hit is very poor
const int PERFECT_AMBIGOUS_SCORE = 50;	// extremely high score that will be kept
const int MIN_AMBIGOUS_SCORE = 20;		// minimum alignemnt score to keep the ambigous record
const int MAX_AMBIGOUS_HIT   = 5;		// maximum alignemnt score of the pair to keep the ambigous record

const int READS_PER_BATCH  = 1 << 18;	// process 256K reads per batch (for parallelization)
const int MAX_SAMLINE_SIZE = 1024;
const int BUFFER_SIZE_PER_BATCH_READ = 1 << 28;	// 256 MB buffer for each thread to convert FASTQ

#endif

