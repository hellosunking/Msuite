#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "util.h"

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Dec 2019
*/


int main( int argc, char *argv[] ) {
	if( argc != 6 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.info> <trim.log> <max.insert.size> <in.sam> <out.prefix>\n";
		cerr << "This program is designed to remove the duplicate reads that have the same start and end/strand.\n\n";
		return 1;
	}

	// loading info file
//	cerr << "Loading genome.info ...\n";
	map<string, map<uint64_t, fraghit> *> samRecord;
	string line, chr;
	stringstream ss;
	ifstream fin( argv[1] );
	if( fin.fail() ) {
		cerr << "Error: could not open file '" << argv[1] << "'!\n";
		exit( 1 );
	}
	while( true ) {
		getline( fin, line );
		if( fin.eof() )break;

		if( line[0] == '#' )continue;

		ss.str( line );
		ss.clear();
		ss >> chr;
//		cerr << "Adding " << chr << " ...\n";
		samRecord.insert( pair<string, map<uint64_t, fraghit>*>(chr, new map<uint64_t, fraghit>()) );
	}
	fin.close();

	// load trim.log
//	cerr << "Preparing memory ...\n";
	unsigned int readNum;
	fin.open( argv[2] );
	if( fin.fail() ) {
		cerr << "Error: could not open file '" << argv[2] << "'!\n";
		exit( 1 );
	}
	getline( fin, line );
	fin.close();
	readNum = atoi( line.c_str() );
//	cerr << "Read number: " << line << " => " << readNum << '\n';
	unsigned int * size = new unsigned int [ readNum ];	// this should be obtained from trim.log!!!
	memset( size, 0, readNum*sizeof(unsigned int) );

	// load sam file
	map<string, map<uint64_t, fraghit> *> :: iterator sam_it;
	map<string, map<uint64_t, fraghit> *> :: iterator no_such_chr = samRecord.end();
	set<unsigned int> dup;
	set<unsigned int> discard;
	map<uint64_t, fraghit> :: iterator hit_it;
	fraghit hit;

	fin.open( argv[4] );
	if( fin.fail() ) {
		cerr << "Error: could not open file '" << argv[4] << "'!\n";
		exit( 1 );
	}

	string line2, cigar1, cigar2, qual1, qual2, unk;
	unsigned int pos1, pos2, fragSize;
	uint64_t key;

	unsigned int lineNum = 0;
//	cerr << "Loading sam file ...\n";
	while( true ) {
		getline( fin, line );
		if( fin.eof() )break;
		getline( fin, line2 );
		++ lineNum;
//		if( ! (lineNum & 0x3fffff) ) {
//			cerr << '\r' << lineNum << " lines loaded.";
//		}

		//499780R1	83	chrX	14710827	42	67M	*	0	0	TCCCAATTCTAAATAGTT	HHHHHHHHHHHHHHHHH XG:Z:GA
		//499780R2	163	chrX	14710378	42	67M	*	0	0	TAACATTTCTTTAATCAC	HHHHHHHH:;:987665 XG:Z:GA

		ss.str( line );
		ss.clear();
		ss >> unk >> unk >> chr >> pos1 >> unk >> cigar1 >> unk >> unk >> unk >> unk >> qual1;
		ss.str( line2 );
		ss.clear();
		ss >> unk >> unk >> chr >> pos2 >> unk >> cigar2 >> unk >> unk >> unk >> unk >> qual2;

//		fprintf( stderr, "Line %d, chr=%s, pos1=%u (0x%x), cigar1=%s, pos2=%u (0x%x), cigar2=%s\n",
//					lineNum, chr.c_str(), pos1, pos1, cigar1.c_str(), pos2, pos2, cigar2.c_str() );

		sam_it = samRecord.find( chr );
		if( sam_it == no_such_chr ) {
			discard.insert( lineNum );
			continue;
		}

		if( line.back() == 'T' ) {	// XG:Z:CT, then pos1 < pos2, then locate the end using pos2 and cigar2
			key = pos1;
			key <<= 32;
			int readLen = get_readLen_from_cigar( cigar2 );
			fragSize = pos2 + readLen - pos1;
			key |= fragSize;
			size[ lineNum ] = fragSize;

//			fprintf( stderr, "=>left=%d, right=%d, readLen=%d, fragSize=%u, key=0x%llx\n",
//							pos1, pos2, readLen, fragSize, key );
		} else {	//XG:Z:GA, then pos1 > pos2, then locate the end using pos1 and cigar1
			key = pos2;
			key <<= 32;
			int readLen = get_readLen_from_cigar( cigar1 );
			fragSize = pos1 + readLen - pos2;
			key |= - fragSize;	// use negative values to mark the strand
			size[ lineNum ] = fragSize;

//			fprintf( stderr, "=>left=%d, right=%d, readLen=%d, fragSize=%u, key=0x%llx\n",
//							pos2, pos1, readLen, fragSize, key );
		}

		register unsigned int score = 0;
		const char * p = qual1.c_str();
		register unsigned int i = 0;
		register unsigned int readLen = qual1.size();
		for( i=0; i != readLen; ++i ) {
			score += p[i];
		}
		p = qual2.c_str();
		readLen = qual2.size();
		for( i=0; i != readLen; ++i ) {
			score += p[i];
		}

		hit_it = sam_it->second->find( key );
		if( hit_it != sam_it->second->end() ) {	// there must be a duplicate
			if( hit_it->second.score >= score ) {	// the previous one is better, mark this one as duplicate
				dup.insert( lineNum );
//				cerr << "Line " << lineNum << " meet a duplicate at line " << hit_it->second.lineNum
//					 << " and the other one is better.\n";
			} else {	// this one is better, then mark the previous one as duplicate and update the record
//				cerr << "Line " << lineNum << " meet a duplicate at line " << hit_it->second.lineNum
//					 << " but this one is better.\n";
				dup.insert( hit_it->second.lineNum );
				hit_it->second.lineNum = lineNum;
				hit_it->second.score = score;
			}
		} else {	// no such record, add this one
			hit.lineNum = lineNum;
			hit.score = score;
			sam_it->second->insert( pair<uint64_t, fraghit>( key, hit ) );
//			cerr << "Line " << lineNum << " is new.\n";
		}
	}
//	cerr << "\rDone: " << lineNum << " lines loaded, dup=" << dup.size() << ", discard=" << discard.size() << ".\n";

//	cerr << "Writing output ...\n";
	// prepare output file
	char * outfile = new char [ MAX_FILE_NAME ];
	sprintf( outfile, "%s.sam", argv[5] );
	ofstream fout( outfile );
	if( fout.fail() ) {
		cerr << "Error: could not write output SAM file!\n";
		exit( 1 );
	}
	// rewind sam file
	fin.clear();
	fin.seekg( ios_base::beg );
	lineNum = 0;
	set<unsigned int> :: iterator non_dup = dup.end();
	set<unsigned int> :: iterator non_discard = discard.end();
	unsigned int unique=0;
	while( true ) {
		getline( fin, line );
		if( fin.eof() )break;
		getline( fin, line2 );
		++ lineNum;
//		if( ! (lineNum & 0x3fffff) ) {
//			cerr << '\r' << lineNum << " lines loaded.";
//		}

		if( discard.find(lineNum)==non_discard && dup.find(lineNum)==non_dup ) {
			fout << line << '\n' << line2 << '\n';
			++ unique;
		} else {
//			cerr << "Line " << lineNum << " is marked as duplicate.\n";
		}
	}
	fin.close();
	fout.close();
//	cerr << "\rDone: " << unique << " lines written.\n";

//	cerr << "Writing log ...\n";
	// write log
	sprintf( outfile, "%s.log", argv[5] );
	fout.open( outfile );
	if( fout.fail() ) {
		cerr << "Error: could not write output log file!\n";
		exit( 1 );
	}
	fout << "All\t"       << lineNum        << '\n'
		 << "Unique\t"    << unique         << '\n'
		 << "Duplicate\t" << dup.size()     << '\n'
		 << "Discard\t"   << discard.size() << '\n';
	fout.close();

//	cerr << "Writing size distribution ...\n";
	// write size distribution
	unsigned int max_size = atoi( argv[3] );
	if( max_size == 0 ) {
		cerr << "Error: Invalid maximum insert size!\n";
		exit( 1 );
	}
	max_size ++;
//	cerr << "Max insert: " << max_size << '\n';
	unsigned int * sizeCnt = new unsigned int [ max_size ];
	for( register unsigned int i=0; i!=max_size; ++i )
		sizeCnt[i] = 0;
	for( register unsigned int i=1; i<=lineNum; ++i ) {
		if( discard.find(i)==non_discard && dup.find(i)==non_dup ) {
//			cerr << "Add size " << size[i] << '\n';
			if( size[i] < max_size ) {
				sizeCnt[ size[i] ] ++;
			} else {	// the CIGAR makes the fragment longer than MAX_INSERT_SIZE !!!
				cerr << "WARNING: Line " << i << " has an unacceptable size (" << size[i]
						<< ") and will be discarded!\n";
			}
		}
	}
	sprintf( outfile, "%s.size.dist", argv[5] );
	fout.open( outfile );
	if( fout.fail() ) {
		cerr << "Error: could not write output size.dist file!\n";
		exit( 1 );
	}
	fout << "Size\tCount\tProportion\n";
	for( register unsigned int i=1; i!=max_size; ++i ) {
		if( sizeCnt[i] != 0 ) {
			fout << i << '\t' << sizeCnt[i] << '\t' << sizeCnt[i]*100.0/unique << '\n';
		}
	}
	fout.close();

//	cerr << "Freeing memory ...\n";
	delete [] size;
	delete [] sizeCnt;
	delete [] outfile;
	for( sam_it=samRecord.begin(); sam_it!=no_such_chr; ++sam_it ) {
		delete sam_it->second;
	}

	return 0;
}


