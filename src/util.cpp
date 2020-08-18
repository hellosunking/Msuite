#include <string>
#include <iostream>
#include "util.h"

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Aug 2020
*/

// load genome from multi-fasta
void loadgenome( const char * file, unordered_map<string, string> & genome ) {
	ifstream fin( file );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << file << " !\n";
		exit(200);
	}
	cout << "Loading genome: " << file << '\n';
	string line, chr, tmp="X";	// X is for position-taking
	register unsigned int i, j;
	while( 1 ) {
		getline( fin, line );
		if( fin.eof() )break;
		j = line.length();
		if( line[0] == '>' ) {
			chr.clear();
			for( i=1; i!=j; ++i ) {
				if( line[i]==' ' || line[i]=='\t' )break;
				chr += line[i];
			}
			genome.insert( pair<string, string>(chr, tmp) );
		} else {
			genome[chr] += line;
		}
	}
	fin.close();
}

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

// update sequence and quality using CIGAR information to support indels
bool fix_cigar( string &cigar, string &realSEQ, string &realQUAL, string &seq, string &qual ) {
	register int i, j, k;
	j = 0;
	register int curr = 0;
	register int cs = cigar.size();
	realSEQ.clear();
	realQUAL.clear();
	for(i=0; i!=cs; ++i) {
		if( cigar[i] <= '9' ) {   // digital
			j *= 10;
			j += cigar[i] - '0';
		} else {	// MUST be M, I, or D
			if( cigar[i] == 'M' ) { // match or mismatch, copy seq and qual
				for(k=0; k!=j; ++k) {
					realSEQ  +=  seq[ curr+k ];
					realQUAL += qual[ curr+k ];
				}
				curr += j;
			} else if ( cigar[i] == 'I' ) { // insertion, discard this part
				curr += j;
			} else if ( cigar[i] == 'D' ) { // deletion, add place holders
				for(k=0; k!=j; ++k) {
					realSEQ  +=  'N';
					realQUAL += '\0';
				}
			} else {	// unsupported CIGAR element
				return false;
			}
			j = 0;
		}
	}
	return true;
}

void call_meth_usage( const char * prg ) {
	cerr << "\nUsage: " << prg << " <mode=SE|PE> <genome.fa> <Msuite.sam> <TAPS|BS> <cycle> <min.score> <output.prefix>\n"
		 << "\nThis program is a component of TAPSuite, designed to call CpG methylation status from SAM file.\n"
		 << "Both SE/PE data are supported; indels are also supported.\n\n"
		 << "The following files will be written for CpG sites:\n"
		 << "\tCpG.meth.call, CpG.meth.bedgraph, chr.count, and M-bias.\n\n"
		 << "If you set --CpH option, the following file will be written:\n"
		 << "\tCpH.meth.call\n\n";
}

