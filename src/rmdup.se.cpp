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
	if( argc != 5 ) {
		cerr << "\nUsage: " << argv[0] << " <chr.info> <trim.log> <in.sam> <out.prefix>\n";
		cerr << "This program is designed to remove the duplicate reads that have the same start and end/strand.\n\n";
		return 1;
	}

	// loading info file
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
//
	// load sam file
	map<string, map<uint64_t, fraghit> *> :: iterator sam_it;
	map<string, map<uint64_t, fraghit> *> :: iterator no_such_chr = samRecord.end();
	set<unsigned int> dup;
	set<unsigned int> discard;
	map<uint64_t, fraghit> :: iterator hit_it;
	fraghit hit;

	fin.open( argv[3] );
	if( fin.fail() ) {
		cerr << "Error: could not open file '" << argv[3] << "'!\n";
		exit( 1 );
	}

	string cigar, qual, unk;
	unsigned int pos;
	uint64_t key;

	unsigned int lineNum = 0;
	while( true ) {
		getline( fin, line );
		if( fin.eof() ) break;
		++ lineNum;

		//499780R1	83	chrX	14710827	42	67M	*	0	0	TCCCAATTCTAAATAGTT	HHHHHHHHHHHHHHHHH XG:Z:GA
		//499780R2	163	chrX	14710378	42	67M	*	0	0	TAACATTTCTTTAATCAC	HHHHHHHH:;:987665 XG:Z:GA

		ss.str( line );
		ss.clear();
		ss >> unk >> unk >> chr >> pos >> unk >> cigar >> unk >> unk >> unk >> unk >> qual;

//		fprintf( stderr, "Line %d, chr=%s, pos1=%u (0x%x), cigar1=%s, pos2=%u (0x%x), cigar2=%s\n",
//					lineNum, chr.c_str(), pos1, pos1, cigar1.c_str(), pos2, pos2, cigar2.c_str() );

		sam_it = samRecord.find( chr );
		if( sam_it == no_such_chr ) {
			discard.insert( lineNum );
			continue;
		}

		key = pos;
		if( line.back() == 'T' ) {	// XG:Z:CT
			key |= READS_CT;
		} else {	//XG:Z:GA
			key |= READS_GA;
		}

		register unsigned int score = 0;
		const char * p = qual.c_str();
		register unsigned int i = 0;
		register unsigned int readLen = qual.size();
		for( i=0; i != readLen; ++i ) {
			score += p[i];
		}
//		fprintf( stderr, "line %d => key=0x%llx, score=%u\n",
//					lineNum, key, score );

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

	// prepare output file
	string outpre = argv[4];
	outpre += ".sam";
	ofstream fout( outpre.c_str() );
	if( fout.fail() ) {
		cerr << "Error: could not write output SAM ile!\n";
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
		++ lineNum;
		//499780R1	83	chrX	14710827	42	67M	*	0	0	TCCCAATTCTAAATAGTT	HHHHHHHHHHHHHHHHH XG:Z:GA
		//499780R2	163	chrX	14710378	42	67M	*	0	0	TAACATTTCTTTAATCAC	HHHHHHHH:;:987665 XG:Z:GA

		if( discard.find(lineNum)==non_discard && dup.find(lineNum)==non_dup ) {
			fout << line << '\n';
			++ unique;
		} else {
//			cerr << "Line " << lineNum << " is marked as duplicate.\n";
		}
	}
	fin.close();
	fout.close();

	// write log
	outpre = argv[4];
	outpre += ".log";
	fout.open( outpre.c_str() );
	fout << "All\t"       << lineNum         << '\n'
		 << "Unique\t"    << unique          << '\n'
		 << "Duplicate\t" << dup.size()      << '\n'
		 << "Discard\t"   << discard.size()  << '\n';
	fout.close();

	return 0;
}


