
#include <string>
using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Dec 2019
*/

#ifndef _MSUITE_UTIL_
#define _MSUITE_UTIL_

const unsigned int MAX_FILE_NAME = 128;

typedef struct {
	unsigned int lineNum;
	unsigned int score;
} fraghit;

const uint64_t READS_GA = 0xf0ULL << 32;
const uint64_t READS_CT = 0x0fULL << 32;

// fix cigar
int get_readLen_from_cigar( const string &cigar );

#endif

