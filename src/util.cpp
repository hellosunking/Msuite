#include <string>
#include "util.h"

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Dec 2019
*/

// fix cigar
int get_readLen_from_cigar( const string &cigar ) {
	register int i, j;
	register int size = 0;
	register int cs = cigar.size();
	register const char * p = cigar.c_str();

	for(i=0, j=0; i!=cs; ++i) {
		if( p[i] <= '9' ) {   // digital
			j *= 10;
			j += p[i] - '0';
		} else {        // MUST be M, I, or D
			if( p[i] == 'M' ) { // match or mismatch, keep
				size += j;
			} else if ( p[i] == 'I' ) { // insertion, ignore
				// do nothing
			} else if ( p[i] == 'D' ) { // deletion, add place holders
				size += j;
			} else {	// unsupported CIGAR element
				return 0;
			}
			j = 0;
		}
	}

	return size;
}


