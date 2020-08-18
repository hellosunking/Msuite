#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <tr1/unordered_map>
#include <stdlib.h>
#include "util.h"

using namespace std;
using namespace std::tr1;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Aug 2020
 * In this version, M-bias data is provided
*/

bool TAPS;
unsigned int MIN_ALIGN_SCORE;	// minimum alignment score for a read to be considered
unsigned int cycle;				// sequencing cycle

// function declarations, the implementation is at the end of this file
// some functions are implemented in util.cpp
void deal_SE_CpG( const char *gfile, const char *samfile, const char *output );
void deal_PE_CpG( const char *gfile, const char *samfile, const char *output );

void callmeth_CpG( string &realSEQ, string &realQUAL, string &chr, int pos, bool strand,
				unordered_map<string, string> :: iterator &git, map<string, map<int, meth>*> &methcall, mbias *mb );
void write_methcall_CpG( map<string, map<int, meth>*> & mc, map<string, int> & chrcount,
				unordered_map<string, string> & g, const char *outfile );

int main( int argc, char *argv[] ) {
	if( argc != 8 ) {
		call_meth_usage( argv[0] );
		return 2;
	}

	string mode = argv[4];
	if( mode=="TAPS" || mode=="taps" ) {
		TAPS = true;
	} else if( mode=="BS" || mode=="bs" ) {
		TAPS = false;
	} else {
		cerr << "Error: Unknown protocol! Must be TAPS or BS!\n";
		exit( 3 );
	}

	cycle = atoi( argv[5] );
	if( cycle == 0 ) {
		cerr << "Error: Invalid cycle!\n";
		exit( 4 );
	}
	MIN_ALIGN_SCORE = atoi( argv[6] );

	mode = argv[1];
	if( mode=="SE" || mode=="se" ) {	// for SE data, only need to calculate target 1
		deal_SE_CpG( argv[2], argv[3], argv[7] );
	} else if( mode=="PE" || mode=="pe" ) {
		deal_PE_CpG( argv[2], argv[3], argv[7] );
	} else {
		cerr << "Error: Unknown mode! Must be PE or SE!\n";
		exit( 5 );
	}

	return 0;
}

// process SE data
void deal_SE_CpG( const char *gfile, const char *samfile, const char *output) {
	// load genome
	unordered_map<string, string> g;
	loadgenome( gfile, g );
	unordered_map<string, string> :: iterator git;
	unordered_map<string, string> :: iterator no_such_chr = g.end();

	map<string, map<int, meth>*> methcall;
	map<string, int> chrcount;
	for( git=g.begin(); git!=g.end(); ++git ) {
		methcall.insert( pair<string, map<int, meth>*>(git->first, new map<int, meth>()) );
		chrcount.insert( pair<string, int>(git->first, 0) );
	}

	mbias *mb = new mbias[ MAX_SAM_LEN ];
	for( register int i=0; i!=MAX_SAM_LEN; ++i ) {
		mb[i].wC=0;mb[i].wT=0;mb[i].wZ=0;
		mb[i].cC=0;mb[i].cT=0;mb[i].cZ=0;
	}

	// open sam file
	ifstream fin( samfile );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << samfile << " to read!\n";
		exit(200);
	}
	cout << "Loading alignment " << samfile << " in SE mode ...\n";
	unsigned int count = 0;
	string line, seqName, chr, cigar, seq, qual;
	int flag;
	string mateinfo, matepos, dist;   //fields that are ignored; all the sequence are converted to WATSON chain
	register unsigned int pos, score;
	stringstream ss;
	string realSEQ, realQUAL;   //these are CIGAR-processed seq and qual
	line.resize( MAX_SAM_LEN );
	bool strand;
	// load sam file
	while( true ) {
		getline( fin, line );
		if( fin.eof() ) break;
		//14_R1	83	chr9	73301642	42	36M	=	73301399	-279	TCCTTCTCTCCCTC	GHHHHHHHHHH	XG:Z:GA

		ss.clear();
		ss.str( line );
		ss >> seqName >> flag >> chr >> pos >> score >> cigar >> mateinfo >> matepos >> dist >> seq >> qual;
		
		if( score < MIN_ALIGN_SCORE ) {
			//cerr << "Discard " << seqName << " due to poor alignment score.\n";
			continue;
		}
		// determine whether the alignemnt is on watson chain or crick chain using the XG:Z: tag
		// which is ALWAYS at the end of read1 for TAPSaligner and Bismark
		if( line.back() == 'T' ) {	// XG:Z:CT => watson
			strand = true;
		} else {
			strand = false;
		}

		git = g.find(chr);
		if( git == no_such_chr ) continue;   // there is NO such chromosome in the genome!!!

		// process the CIGAR, handle the indels
		if( ! fix_cigar(cigar, realSEQ, realQUAL, seq, qual ) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar << ") at line " << line << "!\n";
			continue;
		}

		// call CpG methylation
		callmeth_CpG( realSEQ, realQUAL, chr, pos, strand, git, methcall, mb );
		// chr count
		chrcount.find( chr )->second ++;

		// report progress for every 4 million reads
		++ count;
//		if( ! (count & 0x003fffff) )
//			cout << '\r' << count << " lines loaded.";
	}
	fin.close();
	cout << '\r' << "Done: " << count << " lines loaded.\n";

	cout << "Writing methylation call ...\n";
	write_methcall_CpG( methcall, chrcount, g, output );

	// write M-bias data
	string outfile = output;
	outfile += ".R1.mbias";
	ofstream fmbias( outfile.c_str() );
	if( fmbias.fail() ) {
		cerr << "ERROR: write output M-bias failed.\n";
		exit(20);
	}
	fmbias << "Cycle\twC\twT\twZ\tcC\tcT\tcZ\n";
	for( unsigned int i=0, j=cycle-1; i!=cycle; ++i, --j ) {
		fmbias << i+1 << '\t' << mb[i].wC << '\t' << mb[i].wT << '\t' << mb[i].wZ << '\t'
				<< mb[j].cC << '\t' << mb[j].cT << '\t' << mb[j].cZ << '\n';
	}
	fmbias.close();
	delete [] mb;
}

/////////////////////////////////////////////////////////////////////////////////////////
void deal_PE_CpG( const char *gfile, const char *samfile, const char *output) {
	// load genome
	unordered_map<string, string> g;
	loadgenome( gfile, g );
	unordered_map<string, string> :: iterator git;
	unordered_map<string, string> :: iterator no_such_chr = g.end();

	map<string, map<int, meth>*> methcall;
	map<string, int> chrcount;

	for( git=g.begin(); git!=g.end(); ++git ) {
		methcall.insert( pair<string, map<int, meth>*>(git->first, new map<int, meth>()) );
		chrcount.insert( pair<string, int>(git->first, 0) );
	}
	map<int, meth>* mp;
	map<int, meth> :: iterator mit;
	meth m;

	mbias *mb1 = new mbias[ MAX_SAM_LEN ];
	mbias *mb2 = new mbias[ MAX_SAM_LEN ];
	mbias *mb3 = new mbias[ MAX_SAM_LEN ];	// for fragments with overlap; currently ignored
	for( register int i=0; i!=MAX_SAM_LEN; ++i ) {
		mb1[i].wC = 0;mb1[i].wT = 0;mb1[i].wZ = 0;
		mb1[i].cC = 0;mb1[i].cT = 0;mb1[i].cZ = 0;
		mb2[i].wC = 0;mb2[i].wT = 0;mb2[i].wZ = 0;
		mb2[i].cC = 0;mb2[i].cT = 0;mb2[i].cZ = 0;
		mb3[i].wC = 0;mb3[i].wT = 0;mb3[i].wZ = 0;
		mb3[i].cC = 0;mb3[i].cT = 0;mb3[i].cZ = 0;
	}

	// open sam file
	ifstream fin( samfile );
	if( fin.fail() ) {
		cerr << "Error file: cannot open " << samfile << " to read!\n";
		exit(200);
	}
	cout << "Loading alignment " << samfile << " in PE mode ...\n";

	unsigned int count = 0;
	string line1, line2, seqName, chr, cigar1, seq1, qual1, cigar2, seq2, qual2;
	register int flag;
	string mateinfo, matepos, dist;   //fields that are ignored; all the sequence are converted to WATSON chain
	register unsigned int pos1, pos2, score;
	stringstream ss;
	string realSEQ1, realQUAL1, realSEQ2, realQUAL2;   //these are CIGAR-processed seq and qual
	string mSEQ, mQUAL; //merged sequence and quality if read1 and read2 has overlap
	line1.resize( MAX_SAM_LEN );
	line2.resize( MAX_SAM_LEN );
	mSEQ.resize( MAX_MERGED_SEQ );
	mQUAL.resize( MAX_MERGED_SEQ );
	
	int p1, p2;
	string *r1, *r2, *q1, *q2;
	bool strand;

	// load sam file
	while( true ) {
		getline( fin, line1 );
		if( fin.eof() ) break;
		getline( fin, line2 );
		//14_R1	83	chr9	73301642	42	36M	=	73301399	-279	TCCTCCTTCTCTCCCTC	HHHHHHHHH	XG:Z:CT
		//14_R2	163	chr9	73301399	42	36M	=	73301642	279	TTTATTTTGATCCTGTA	DDCBA@?>=<;986420.

		ss.clear();
		ss.str( line1 );
		ss >> seqName >> flag >> chr >> pos1 >> score >> cigar1 >> mateinfo >> matepos >> dist >> seq1 >> qual1;
		//cerr << seqName << '\n';

		if( score < MIN_ALIGN_SCORE ) {
			//cerr << "Discard " << seqName << " due to poor alignment score.\n";
			continue;
		}

		// determine whether the alignemnt is on watson chain or crick chain using the XG:Z: tag
		// which is ALWAYS at the end of read1 for TAPSaligner and Bismark
		if( line1.back() == 'T' ) {  // XG:Z:CT => watson
			strand = true;
		} else {	// XG:Z:GA => crick
			strand = false;
		}

		git = g.find( chr );
		if( git == no_such_chr ) continue;   // there is NO such chromosome in the genome!!!

		ss.clear();
		ss.str( line2 );
		ss >> seqName >> flag >> chr >> pos2 >> score >> cigar2 >> mateinfo >> matepos >> dist >> seq2 >> qual2;

		// process CIGAR 1, handle the indels
		realSEQ1.clear();
		realQUAL1.clear();
		if( ! fix_cigar( cigar1, realSEQ1, realQUAL1, seq1, qual1 ) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar1 << ") in " << seqName << "!\n";
			continue;
		}
		// process CIGAR 2, handle the indels
		realSEQ2.clear();
		realQUAL2.clear();
		if( ! fix_cigar( cigar2, realSEQ2, realQUAL2, seq2, qual2 ) ) {
			cerr << "ERROR: Unsupported CIGAR (" << cigar2 << ") in " << seqName << "!\n";
			continue;
		}

		if( pos1 <= pos2 ) {
			p1 = pos1; r1 = & realSEQ1; q1 = & realQUAL1;
			p2 = pos2; r2 = & realSEQ2; q2 = & realQUAL2;
		} else {
			p1 = pos2; r1 = & realSEQ2; q1 = & realQUAL2;
			p2 = pos1; r2 = & realSEQ1; q2 = & realQUAL1;
		}

		if( p1 + r1->size() <= p2 ) { //there is NO overlap
			callmeth_CpG( realSEQ1, realQUAL1, chr, pos1, strand, git, methcall, mb1 );
			callmeth_CpG( realSEQ2, realQUAL2, chr, pos2, strand, git, methcall, mb2 );
		} else {	// there is overlap in read 1 and read 2
			//cerr << "Found overlap in " << seqName << '\n';
			if( p2+r2->size() >= p1+r1->size() ) {	// most case
				mSEQ.clear();
				mQUAL.clear();
				int len = p2 + r2->size() - p1;
				//cerr << "len\t" << len << "\tp1=" << p1 << "\tp2=" << p2 << '\n';
				register int rs = r1->size();
				register int offset = p2 - p1;
				register unsigned int k;
				for( k=0; k != offset; ++k ) {	// read1 only
					mSEQ  += r1->at(k);		// can also use substr
					mQUAL += q1->at(k);
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";

				for( ; k != rs; ++k ) {	// overlapped region, peak the one with higher quality
					if( q1->at(k) >= q2->at(k-offset) ) {
						mSEQ  += r1->at(k);
						mQUAL += q1->at(k);
					} else {
						mSEQ  += r2->at(k-offset);
						mQUAL += q2->at(k-offset);
					}
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";

				for( ; k != len; ++k ) {	//read2 only
					mSEQ  += r2->at(k-offset);
					mQUAL += q2->at(k-offset);
				}
				//cerr << "mSEQ\t" << mSEQ << "\n";
				callmeth_CpG( mSEQ, mQUAL, chr, p1, strand, git, methcall, mb3 );
			} else {	// rare case that R1 completely contains R2 => use R1 directly
				callmeth_CpG( realSEQ1, realQUAL1, chr, pos1, strand, git, methcall, mb3 );
			}
		}

		chrcount.find( chr )->second ++;
		++ count;
//		if( ! (count & 0x003fffff) )
//			cout << '\r' << count << " lines loaded.";
	}
	fin.close();
	cout << '\r' << "Done: " << count << " lines loaded.\n";

	write_methcall_CpG( methcall, chrcount, g, output );

	// write M-bias data
	string outfile = output;
	outfile += ".R1.mbias";
	ofstream fmbias( outfile.c_str() );
	if( fmbias.fail() ) {
		cerr << "ERROR: write output M-bias failed.\n";
		exit(20);
	}
	fmbias << "Cycle\twC\twT\twZ\tcC\tcT\tcZ\n";
	for( unsigned int i=0, j=cycle-1; i!=cycle; ++i, --j ) {
		fmbias << i+1 << '\t' << mb1[i].wC << '\t' << mb1[i].wT << '\t' << mb1[i].wZ << '\t'
				<< mb1[j].cC << '\t' << mb1[j].cT << '\t' << mb1[j].cZ << '\n';
	}
	fmbias.close();

	outfile = output;
	outfile += ".R2.mbias";
	fmbias.open( outfile.c_str() );
	if( fmbias.fail() ) {
		cerr << "ERROR: write output M-bias failed.\n";
		exit(20);
	}
	fmbias << "Cycle\twC\twT\twZ\tcC\tcT\tcZ\n";
	for( unsigned int j=0, i=cycle-1; j!=cycle; ++j, --i ) {
		fmbias << j+1 << '\t' << mb2[i].wC << '\t' << mb2[i].wT << '\t' << mb2[i].wZ << '\t'
				<< mb2[j].cC << '\t' << mb2[j].cT << '\t' << mb2[j].cZ << '\n';
	}
	fmbias.close();
	delete [] mb1;
	delete [] mb2;
	delete [] mb3;
}

// call meth from sequence
void callmeth_CpG( string &realSEQ, string &realQUAL, string &chr, int pos, bool strand,
				unordered_map<string, string> :: iterator &git, map<string, map<int, meth>*> &methcall, mbias *mb ) {
	map<int, meth>* mp;
	map<int, meth> :: iterator mit;
	meth m;
	unsigned int rs = realSEQ.size();
	char c1, c2, c3;
	unsigned int i = 0;
	unsigned int j = pos + i;
	for( ; i!=rs; ++i, ++j) {
//		if( realQUAL[i] < MIN_QUAL_SCORE )
//			continue;
		if( strand ) {	// watson strand
			if( j == git->second.size() - 1 )
				continue;

			c1 = git->second[j];
			c2 = git->second[j+1];

			if( (c1!='C' && c1!='c') || (c2!='G' && c2!='g') )	// not a CpG site
				continue;

			mp  = methcall.find( chr )->second;
			mit = mp->find( j );
			if( mit == mp->end() ) {	// no record, insert one
				if( realSEQ[i] == 'C' ) {
					m.wC=1; m.wT=0; m.wZ=0; m.cC=0; m.cT=0; m.cZ=0;
					mb[i].wC ++;
				} else if( realSEQ[i] == 'T' ) {
					m.wC=0; m.wT=1; m.wZ=0; m.cC=0; m.cT=0; m.cZ=0;
					mb[i].wT ++;
				} else {
					m.wC=0; m.wT=0; m.wZ=1; m.cC=0; m.cT=0; m.cZ=0;
					mb[i].wZ ++;
				}
				mp->insert( pair<int, meth>(j, m) );
			} else {
				if( realSEQ[i] == 'C' ) {
					mit->second.wC ++;
					mb[i].wC ++;
				} else if( realSEQ[i] == 'T' ) {
					mit->second.wT ++;
					mb[i].wT ++;
				} else {
					mit->second.wZ ++;
					mb[i].wZ ++;
				}
			}
		} else {	// check G in CpG for crick strand reads
			if( j == 0 )
				continue;

			c1 = git->second[j-1];
			c2 = git->second[j];

			if( (c1!='C' && c1!='c') || (c2!='G' && c2!='g') )	// not a CpG site
				continue;

			mp  = methcall.find( chr )->second;
			mit = mp->find( j-1 );
			if( mit == mp->end() ) {	// no record, insert one
				if( realSEQ[i] == 'G' ) {
					m.cC=1; m.cT=0; m.cZ=0; m.wC=0; m.wT=0; m.wZ=0;
					mb[i].cC ++;
				} else if( realSEQ[i] == 'A' ) {
					m.cC=0; m.cT=1; m.cZ=0; m.wC=0; m.wT=0; m.wZ=0;
					mb[i].cT ++;
				} else {
					m.cC=0; m.cT=0; m.cZ=1; m.wC=0; m.wT=0; m.wZ=0;
					mb[i].cZ ++;
				}
				mp->insert( pair<int, meth>(j-1, m) );
			} else {	// there is such a record
				if( realSEQ[i] == 'G' ) {
					mit->second.cC ++;
					mb[i].cC ++;
				} else if( realSEQ[i] == 'A' ) {
					mit->second.cT ++;
					mb[i].cT ++;
				} else {
					mit->second.cZ ++;
					mb[i].cZ ++;
				}
			}
		}
	}
}

// write meth call into file
void write_methcall_CpG( map<string, map<int, meth>*> & mc, map<string, int> & chrcount,
						unordered_map<string, string> & g, const char *outpre ) {
	string outfile = outpre;
	outfile += ".CpG.meth.call";
	ofstream fcpg( outfile.c_str() );
	if( fcpg.fail() ) {
		cerr << "ERROR: write output CpG meth call failed.\n";
		exit(20);
	}

	outfile = outpre;
	outfile += ".CpG.meth.bedgraph";
	ofstream fbed( outfile.c_str() );
	if( fbed.fail() ) {
		cerr << "ERROR: write output bedgraph failed.\n";
		exit(21);
	}

	outfile = outpre;
	outfile += ".CpG.meth.log";
	ofstream flog( outfile.c_str() );
	if( flog.fail() ) {
		cerr << "ERROR: write output failed.\n";
		exit(21);
	}

	unordered_map<string, string> :: iterator git;
	map<string, map<int, meth>*> :: iterator mit;	// methcall iterator
	map<int, meth> :: iterator cit;	// methcall iterator for each chromosome
	map<string, int> :: iterator chrit;

	meth m;

	fcpg << "#chr\tLocus\tTotal\twC\twT\twOther\tContext\tcC\tcT\tcOther\n";
	flog << "#chr\tNo.Reads\tCpG.wC\tCpG.wT\tCpG.cC\tCpG.cT\n";
	for( mit=mc.begin(); mit!=mc.end(); ++mit ) {
		git = g.find( mit->first );
		int totalWC=0, totalWT=0, totalCC=0, totalCT=0;	//total C, T
		int total_CpG_WC=0, total_CpG_WT=0, total_CpG_CC=0, total_CpG_CT=0;	//total C, T on CpG sites

		for( cit=mit->second->begin(); cit!=mit->second->end(); ++cit ) {
			int i = cit->first;
			m = cit->second;
			unsigned int Valid = m.wC+m.wT+m.cC+m.cT;

			if( (git->second[i]=='C' || git->second[i]=='c') && (git->second[i+1]=='G' || git->second[i+1]=='g') ) {	// CpG sites
				fcpg << mit->first << '\t' << i << '\t' << Valid+m.wZ+m.cZ << '\t'
					 << m.wC << '\t' << m.wT << '\t' << m.wZ << '\t'
					 << git->second[i-1] << git->second[i] << git->second[i+1] << git->second[i+2] << '\t'
					 << m.cC << '\t' << m.cT << '\t' << m.cZ << '\n';

				if( Valid != 0 ) {
					float md;
					if( TAPS ) {
						md = (m.wT+m.cT)*100.0/Valid;
					} else {
						md = (m.wC+m.cC)*100.0/Valid;
					}
					fbed << mit->first << '\t' << i-1 << '\t' << i << '\t' << md << '\n';
				}
				total_CpG_WC += m.wC;
				total_CpG_WT += m.wT;
				total_CpG_CC += m.cC;
				total_CpG_CT += m.cT;
			}
		}

		chrit = chrcount.find( mit->first );
		flog << mit->first << '\t' << chrit->second << '\t'
			 << total_CpG_WC << '\t' << total_CpG_WT << '\t'
			 << total_CpG_CC << '\t' << total_CpG_CT << '\n';

		delete mit->second;
	}
	fcpg.close();
	flog.close();
	fbed.close();
}

