#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <tr1/unordered_map>
#include <stdlib.h>

using namespace std;
using namespace std::tr1;


/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Dec 2019
*/

#ifndef _MSUITE_UTIL_
#define _MSUITE_UTIL_

const unsigned int MAX_FILE_NAME = 128;

const uint64_t READS_GA = 0xf0ULL << 32;
const uint64_t READS_CT = 0x0fULL << 32;

const unsigned int MAX_SAM_LEN		= 4096;	// maximum length for fragment size and sam line
const unsigned int MAX_CpG_COVER	=  512;	// maximum CpG coverage of a fragment
const unsigned int MIN_QUAL_SCORE	=   33;	// minimum phred score for a CpG site to be considered;
                                            // Note: this parameter is not allowed to set by the user in the current version
const unsigned int MAX_MERGED_SEQ	=  512;

typedef struct {
	unsigned int lineNum;
	unsigned int score;
}fraghit;

// methylation call
typedef struct {
	unsigned short wC;	// 'C' on watson chain
	unsigned short wT;	// 'T' on watson chain
	unsigned short wZ;	// neither 'C' nor 'T'; SNPs or sequencing errors
	unsigned short cC;	// 'C' on crick chain
	unsigned short cT;	// 'T' on crick chain
	unsigned short cZ;	// neither 'C' nor 'T'; SNPs or sequencing errors
}meth;

// M-bias
typedef struct {
	unsigned int wC;
	unsigned int wT;
	unsigned int wZ;
	unsigned int cC;
	unsigned int cT;
	unsigned int cZ;
}mbias;

// load genome from multi-fasta
void loadgenome( const char * file, unordered_map<string, string> & genome );
// fix cigar
int get_readLen_from_cigar( const string &cigar );
bool fix_cigar(string &cigar, string &realSEQ, string &realQUAL, string &seq, string &qual);

// usage information for meth.call
void call_meth_usage( const char * prg );

#endif

