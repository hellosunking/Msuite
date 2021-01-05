#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <string>
#include <omp.h>
#include <unistd.h>
#include "common.h"

using namespace std;

/*
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package
 * Date: Dec 2019
 *
 * In Ktrim, the conversion log format is placed BEFORE the raw seqName. Here is the rule:
 *   read 1:
 *       @LINE_NUMBER+S|xx;xx$  => for reads with endC and C>T changes, 'S' is its quality score
 *       @LINE_NUMBER+xx;xx$    => for reads without endC
 *       @LINE_NUMBER+$         => for reads without endC and conversions
 *   read 2:
 *       @S|xx;xx$  => for reads with frontG and G>A changes, 'S' is its quality score
 *       @xx;xx$    => for reads without frontG
 *       @$         => for reads without frontG and conversions
 *
 * This could fasten the de-convert procedure since the whole seqName is untouched
 */

int main( int argc, char *argv[] ) {
	if( argc < 5 ) {
		cerr << "\nUsage: " << argv[0] << " <CG2TG.sam> <CG2CA.sam> <trim.log> <output.sam> [thread=1]\n"
			 << "\nThis program is part of Msuite, designed to generate the final alignment file.\n\n"
			 << "This program will convert 'T' back to 'C' in the alignment file.\n"
			 << "For ambigous reads, only those with unique best hits and have good scores are kept.\n"
			 << "Align score cutoff for ambigous reads: " << MIN_ALIGN_SCORE_AMB << '\n'
			 << "Align score cutoff for non-ambigous reads: " << MIN_ALIGN_SCORE_UNQ << "\n\n"
			 << "Multi-thread is supported, 4-8 threads are recommended.\n\n";
		//cerr << "Rescue mode is ON.\n\n";

		return 2;
	}

    int thread = 1;
	if( argc > 5 ) {
		thread = atoi( argv[5] );
		if( thread < 0 ) {  // use all threads
			cerr << "Warning: thread is set to 0! I will use all threads instead (maximum 8).\n";
			thread = omp_get_max_threads();
			if( thread > 8 )
				thread = 8;
		}
	}

	// load Ktrim.log, get the read number
	FILE* ktrimlog = fopen( argv[3], "r" );
	if( ktrimlog == NULL ) {
		cerr << "Error: cannot open file " << argv[3] << "!\n";
		exit(10);
	}
	unsigned int readNum;
	int fscanfRet = fscanf( ktrimlog, "%d", &readNum );
	fclose( ktrimlog );
	//printf( "Line number: %d!\n", readNum ); 

    //cerr << "Requesting memory ...\n";
	++ readNum;
	int *hits = new int[ readNum ];
	// initialization
	for(register int i=0; i!=readNum; ++i)
		hits[i] = IMPOSSIBLE_AS_SCORE;

	// prepare files
	ifstream CG2TG( argv[1] );
	if( CG2TG.fail() ) {
		cerr << "Error: cannot open file " << argv[1] << "!\n";
		exit(11);
	}
	ifstream CG2CA( argv[2] );
	if( CG2CA.fail() ) {
		cerr << "Error: cannot open file " << argv[2] << "!\n";
		CG2TG.close();
		exit(12);
	}
	// check output file, delete it if exists
	if( access(argv[4], F_OK) == 0 ) { // file exists, delete it!
		cerr << "WARNING: output file exists! It will be overwritten!\n";
		if( remove(argv[4]) != 0 ) {
			cerr << "FATAL ERROR: Could not delete file " << argv[4] << "!\n";
			return 100;
		}
	}
	FILE * outsam = fopen( argv[4], "a" );	//this file must be EMPTY or NULL before calling this program!!!
	if( outsam == NULL ) {
		cerr << "Error: cannot open file " << argv[4] << " to write!\n";
		CG2TG.close();
		CG2CA.close();
		exit(13);
	}

	// prepare Hits arrary
//	cout << "Preparing HITS array ...\n";
	// load align.CG2TG.sam, record the hit information
	string *R1 = new string[ READS_PER_BATCH ];
	while( true ) {
		unsigned int loaded = 0;
		while( true ) {
			getline( CG2TG, R1[ loaded ] );
			if( CG2TG.eof() )break;

			++ loaded;
			if( loaded == READS_PER_BATCH )
				break;
		}
//		cerr << loaded << " lines loaded\n";
		if( loaded == 0 ) break;

		omp_set_num_threads( thread );
		// get line number
		#pragma omp parallel
		{
			unsigned int tn = omp_get_thread_num();
			unsigned int start = loaded * tn / thread;
			unsigned int end   = loaded * (tn+1) / thread;
			const char * p;
//          cerr << "Thread " << tn << ": s=" << start << ", e=" << end << '\n';

			// Bowtie2 records the mismatch-score in the AS:i: tag and next match (if any) in XS:i: tag
			// we will extract AS:i: tag in both read1 and read2
			// 
			// 10753_chr1:226691056-226691358_R1    99  chr1    226691056   35  50M =   226691309   303
			// GTCTTTACAATTTGGCATGTTTTTGCAGTGGCTGGTACCAGTTGTTCCTT  HHHHHHHHHHHHHHGGGGGFFFFEEDDCCBBA@@?>=<<;:8765321/.
			// AS:i:0  XS:i:-6 XN:i:0  XM:i:0  XO:i:0      XG:i:0  NM:i:0  MD:Z:50 YS:i:0  YT:Z:CP
			// 10753_chr1:226691056-226691358_R2    147 chr1    226691309   35  50M =   226691056   -303
			// TTCTGCTGAGAGATCTGCTGTTAGTCTGATGGGTTTCCCTTTGTGGGTAA  ./1235678:;<<=>?@@ABBCCDDEEFFFFGGGGGHHHHHHHHHHHHHH
			// AS:i:0  XS:i:-6 XN:i:0  XM:i:0  XO:i:0      XG:i:0  NM:i:0  MD:Z:50 YS:i:0  YT:Z:CP
			//
			// Mostly AS-scores are 0 or negative
			for( unsigned int ii=start; ii!=end; ++ii ) {
				p = R1[ii].c_str();
                register int i;
				register unsigned int j = 0;
				for( i=0; p[i]!=LINE_NUMBER_SEPARATOR; ++i ) {
					j <<= 4;
					if( p[i]<='9' ) {	//0-9
						j += p[i] - '0';
					} else {	//a-f
						j += p[i] - 87;	// 'a'-10
					}
				}

				// skip the data before AS:i: tag (which is on column 12 (1-based)
				// therefore there are 11 TABS before it)
				register int k = 0;
				++ i;
				while( true ) {
					if( p[i] == '\t' ) {
						++ k;
						if( k == 11 )
							break;
					}
					++ i;
				}
				int AStag1value = 0;
				bool AStag1flag;
				i += 6; // move from '\t' to AS:i:
				if( p[i] == '-' ) {
					AStag1flag = true;
					++ i;
				} else {
					AStag1flag = false;
				}
				while( p[i] != '\t' ) {
					AStag1value *= 10;
					AStag1value += p[i] - '0';
					++ i;
				}
				if( AStag1flag )
				  AStag1value = - AStag1value;

				hits[j] = AStag1value;
//				cerr << "j=" << j << ", score=" << hits[j] << '\n';
			}
		}
		if( CG2TG.eof() )break;
	}
	//cerr << "Done.\n";

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// process CG2CA.sam file
	// In this file, read1 is ALWAYS on crick chain and read2 is always on WATSON chain
//	cout << "Processing CG2CA ...\n";
	int *cntCA = new int [ thread ];
	for( unsigned int i=0; i!=thread; ++i )
		cntCA[i] = 0;

	while( true ) {
		unsigned int loaded = 0;
		while( true ) {
			getline( CG2CA, R1[ loaded ] );
			if( CG2CA.eof() )break;

			++ loaded;
			if( loaded == READS_PER_BATCH ) break;
		}
		//cerr << loaded << "lines loaded\n";

		// start task
		omp_set_num_threads( thread );
		#pragma omp parallel
		{
			unsigned int tn = omp_get_thread_num();
			unsigned int sindex = loaded * tn / thread;
			unsigned int eindex = loaded * (tn+1) / thread;

			const char *p;	// *IDstr, *CIGARstr;
			register unsigned int IDstart, IDend;
			register bool endC, frontG, posStrand;	// indicators: "C" at the end, "G" at the front, and WATSON strand
			register char Qend;	// quality score for the "C"
			stringstream ss;
			register char * seqName = (char *) malloc( MAX_SEQNAME_SIZE );
			char chr[  MAX_ITERM_SIZE ];
			char flag[ MAX_ITERM_SIZE ];
			string cigar, seq, qual;
			register int pos, score;
			// the following items will not be updated, therefore string is OK
			char scoreSTR[ MAX_ITERM_SIZE ];
			char mateinfo[ MAX_ITERM_SIZE ];
			register int len, bias;

			for( unsigned int ii=sindex; ii!=eindex; ++ii ) {
                // deal read 1: should be on CRICK chain
                ss.clear();
                ss.str( R1[ii] );
                // extract seqName to score first; score is needed in case there is another hit on CT2TG
                // this fragment could be discarded, so no need to extract others here
                ss >> seqName >> flag >> chr >> pos >> score;

                // deal seqName
				register int i, k;
                register unsigned int j = 0;
                for( i=0; seqName[i] != LINE_NUMBER_SEPARATOR; ++i ) {
                    j <<= 4;	// the number is in HEX
                    if( seqName[i] <= '9' ) {	//0-9
                        j += seqName[i] - '0';
                    } else {	//a-f
                        j += seqName[i] - 87;	// 'a'-10
                    }
                }

				// check whether it is an ambigous hit
				if( hits[j] != IMPOSSIBLE_AS_SCORE ) {	// this is an ambigous hit
					register int ASindex = i + 1;
					p = R1[ii].c_str();
					k = 0;
					while( true ) {
						if( p[ASindex] == '\t' ) {
							++ k;
							if( k == 11 )
							  break;
						}
						++ ASindex;
					}
					int AStagvalue = 0;
					bool AStagflag;
					ASindex += 6;
					if( p[ASindex] == '-' ) {
						AStagflag  = true;
						++ ASindex;
					} else {
						AStagflag  = false;
					}
					while( p[ASindex] != '\t' ) {
						AStagvalue *= 10;
						AStagvalue += p[ASindex] - '0';
						++ ASindex;
					}
					if( AStagflag ) AStagvalue = - AStagvalue;

					ASindex = AStagvalue;	// AS-score of this hit

//					cerr << "Ambigous: j=" << j << ", score=" << ASindex << '\n';

					if( ASindex > hits[j] ) {	// this one is the unique best hit
						hits[j] = DISCARD_AMB_HIT;
						if( score < MIN_ALIGN_SCORE_AMB )
							continue;
						score >>= 1;	// lower the score
					} else if ( ASindex == hits[j] ) {	// discard both (non-unique best hits)
						hits[j] = DISCARD_AMB_HIT;
						continue;
					} else {	// the other one is the unique best hit
						hits[j] = AMB_HIT_MARKER;
						continue;
					}
				} else {	// non-ambigous
					if( score < MIN_ALIGN_SCORE_UNQ ) {	// too poor quality, discard
						// TODO: if both R1 and R2 do not have 2nd hits (i.e., no XS tag), still keep it
						// if XS tag exists, it will be always after the AS tag, i.e., on column 13
						//cerr << "Checking " << seqName << '\n';
						register int XSindex = i + 1;
						p = R1[ii].c_str();
						k = 0;
						while( true ) {
							if( p[XSindex] == '\t' ) {
								++ k;
								if( k == 12 )
								  break;
							}
							++ XSindex;
						}
						if( p[XSindex+1]=='X' && p[XSindex+2]=='S' )	// there is a XS index, DISCARD
							continue;
						//cerr << "UNIQUE_POOR\t" << seqName << '\n';
						//continue;

						// here this read will be kept !!!
						//cerr << "Rescued: " << seqName << '\n';
					}
				}

                //this fragment will be kept, extract other information
                ss >> cigar >> mateinfo >> mateinfo >> mateinfo >> seq >> qual;

                // look for the conversion log start
                if( seqName[i+2] == KEEP_QUAL_MARKER ) {	// marker for endC
                    Qend = seqName[i+1];
                    endC = true;
                    i += 3;
                } else {
                    endC = false;
                    ++ i;
                }
                if( seqName[i] != CONVERSION_LOG_END ) {	// there are C>T conversions
                    // process the conversion log
                    len = seq.size() - 1;
                    // the counting is from the right to the left and the missing 'C' (if exists)
                    // should be the leftmost, therefore it DOES not affect calculating the index
                    j = 0;
                    for( ; seqName[i] != CONVERSION_LOG_END; ++i ) {
                        if ( seqName[i] == CONVERSION_LOG_SEPARATOR ) {	// one change meet
                            seq[len-j] = 'G';
                            j = 0;
                        } else {
                            j <<= 4;
                            if( seqName[i] <= '9' ) {	// 0-9
                                j += seqName[i] - '0';
                            } else {	// a-f
                                j += seqName[i] - 87;	// 'a'-10
                            }
                        }
                    }
                    seq[len-j] = 'G';
                }
                IDstart = i + 1;

                // now seqName is processed, then deal with pos and CIGAR;
                // matepos is NOT considered here!!! Dist is always not affected, so ignore it
                if( endC ) {
                    // for CG2CA, read1 is always on crick chain; so if there is an endC,
                    // add '1M' at the beginning of the CIGAR and update POS
                    -- pos;
                    j = 0;
                    p = cigar.c_str();
                    for(i=0; p[i] != 'M'; ++i) {
                        j *= 10;
                        j += p[i] - '0';
                    }
                    ++ j;	// this is to add the '1M' at the beginning of the CIGAR
					// the XG:Z:GA is to mark that this fragment is aligned to the crick chain
					// this information is used in meth.caller and is consistent with Bismark
					// this mark is ONLY added to read1
                    fprintf( outsam, "%s\t%s\t%s\t%d\t%d\t%d%s\t*\t0\t0\tG%s\t%c%s\tXG:Z:GA\n",
								seqName+IDstart, flag, chr, pos, score, j, p+i,
								seq.c_str(), Qend, qual.c_str() );
                } else {	// do not need to update CIGAR and pos
                    fprintf( outsam, "%s\t%s\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tXG:Z:GA\n",
								seqName+IDstart, flag, chr, pos, score, cigar.c_str(),
								seq.c_str(), qual.c_str() );
                }

                ++ cntCA[tn];
            }
            free( seqName );
        }

        if( CG2CA.eof() )break;
	}
    unsigned int CAall = 0;
    for( unsigned int i=0; i!=thread; ++i )
        CAall += cntCA[i];

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// process CG2TG.sam
	// In this file, read1 is ALWAYS on WATSON chain and read2 is always on CRICK chain
//	cout << "Processing CG2TG ...\n";
	CG2TG.clear();
	CG2TG.seekg( ios_base::beg );
    int *cntGT = new int [ thread ];
    for( unsigned int i=0; i!=thread; ++i )
        cntGT[i] = 0;
    while( true ) {
        unsigned int loaded = 0;
		while( true ) {
			getline( CG2TG, R1[ loaded ] );
			if( CG2TG.eof() )break;

			++ loaded;
			if( loaded == READS_PER_BATCH )
				break;
		}
        //cerr << loaded << "lines loaded\n";

        // start task
		omp_set_num_threads( thread );
		#pragma omp parallel
		{
            unsigned int tn = omp_get_thread_num();
			unsigned int sindex = loaded * tn / thread;
			unsigned int eindex = loaded * (tn+1) / thread;

			const char *p;	//*IDstr, *CIGARstr;
			register unsigned int IDstart, IDend;
			register bool endC, frontG, posStrand;	// indicators: "C" at the end, "G" at the front, and WATSON strand
			register char Qend;	// quality score for the "C"
			stringstream ss;
			register char * seqName = (char *) malloc( MAX_SEQNAME_SIZE );
			char chr[  MAX_ITERM_SIZE ];
			char flag[ MAX_ITERM_SIZE ];
			string cigar, seq, qual;
			register int pos, score;
			//the following items will not be updated, therefore string is OK
			char scoreSTR[ MAX_ITERM_SIZE ];
			char mateinfo[ MAX_ITERM_SIZE ];
			register int len, bias;

            for( unsigned int ii=sindex; ii!=eindex; ++ii ) {
                // deal read 1: should be on WATSON chain
                ss.clear();
                ss.str( R1[ ii ] );
                ss >> seqName >> flag >> chr >> pos >> score;

                // deal seqName
				register int i, k;
                register unsigned int j = 0;
                for( i=0; seqName[i] != LINE_NUMBER_SEPARATOR; ++i ) {
                    j <<= 4;	// the number is in HEX
                    if( seqName[i]<='9' ) {	//0-9; seqName > '0' is always true
                        j += seqName[i] - '0';
                    } else {	//a-f
                        j += seqName[i] - 87;	// 'a'-10
                    }
				}

				// check ambigous
				if( hits[j] == DISCARD_AMB_HIT ) { // ambigous marked and discarded by CG2CA
					continue;
				} else if ( hits[j] == AMB_HIT_MARKER ) {	// ambigous marked BUT kept by CG2CA
					if( score < MIN_ALIGN_SCORE_AMB )
						continue;
					score >>= 1;	// lower the score
				} else {	// unique hit
					if( score < MIN_ALIGN_SCORE_UNQ ) {
						// TODO: if R1 does not have 2nd hits (i.e., no XS tag), still keep it
						// if XS tag exists, it will be always after the AS tag, i.e., on column 13
						// It seems that this will introduce lots of errors!!!
						//cerr << "Checking " << seqName << '\n';
						register int XSindex = i + 1;
						p = R1[ii].c_str();
						k = 0;
						while( true ) {
							if( p[XSindex] == '\t' ) {
								++ k;
								if( k == 12 )
								  break;
							}
							++ XSindex;
						}
						if( p[XSindex+1]=='X' && p[XSindex+2]=='S' )	// there is a XS index, DISCARD
							continue;

						//cerr << "UNIQUE_POOR\t" << seqName << '\n';
						//continue;
						//now this hit will be kept!!!
						//cerr << "Rescued: " << seqName << '\n';
					}
				}

                // extract the other sections in the SAM record
                ss >> cigar >> mateinfo >> mateinfo >> mateinfo >> seq >> qual;

                // check endC marker
                if( seqName[i+2] == KEEP_QUAL_MARKER ) {	// there is a C at the end
                    Qend = seqName[i+1];
                    endC = true;
                    i += 3;
                } else {	// no C at the end
                    endC = false;
                    ++ i;
                }

                // process the conversion log
                if( seqName[i] != CONVERSION_LOG_END ) {	// there are changes in this read
                    j = 0;
                    for( ; seqName[i] != CONVERSION_LOG_END; ++i ) {
                        if ( seqName[i] == CONVERSION_LOG_SEPARATOR ) {	// one change met
                            seq[j] = 'C';	//convert back to 'C'; read 1 is on WATSON
                            j = 0;
                        } else {
                            j <<= 4;
                            if( seqName[i]<='9' ) {	//0-9, seqName > '0' is always true
                                j += seqName[i] - '0';
                            } else {	//a-f
                                j += seqName[i] - 87;	// 'a'-10
                            }
                        }
                    }
                    seq[j] = 'C';	//convert back to C
                }
                IDstart = i + 1;

                // now seqName is processed, then deal with pos and CIGAR;
                // matepos is NOT considered here!!! Dist is always not affected, so ignore it
                if( endC ) {
                    // for CG2TG, read1 is always on WATSON chain, so add 1M to the end of CIGAR; do not update POS
                    // CIGAR: xM[yID]zM
                    len = cigar.size();
                    p = cigar.c_str();
                    i = len - 3;	// len-1 is 'M', len-2 MUST be a digital
                    while( i >= 0 ) {
                        if( p[i] > '9' ) break;	// it is not a digital (should be I/D/), stop here
                        -- i;
                    }
                    ++ i;	// now it point to the first digital of the LAST match segment; could be 0 (i.e., cigar is xxM only)
                    j = 0;
                    for( k=i; p[k] != 'M'; ++k ) {
                        j *= 10;
                        j += p[k] - '0';
                    }
                    ++ j;	// this is to add the 1M at the end of CIGAR
                    cigar.resize( i );
                    /*fout << seqName+IDstart << '\t' << flag << '\t' << chr << '\t' << pos << '\t' << scoreSTR
                        << '\t' << cigar << j << 'M' << '\t' << mateinfo << '\t' << matepos << '\t' << dist
                        << '\t' << seq << 'C' << '\t' << qual << Qend << endl;*/
					// the XG:Z:CT is to mark that this fragment is aligned to the watson chain
					// this information is used in meth.caller and is the same as Bismark
					// this mark is ONLY added to read1
                    fprintf( outsam, "%s\t%s\t%s\t%d\t%d\t%s%dM\t*\t0\t0\t%sC\t%s%c\tXG:Z:CT\n",
								seqName+IDstart, flag, chr, pos, score, cigar.c_str(), j,
								seq.c_str(), qual.c_str(), Qend );
                } else {	// do not need to update CIGAR
                    /*fout << seqName+IDstart << '\t' << flag << '\t' << chr << '\t' << pos << '\t' << scoreSTR
                        << '\t' << cigar << '\t' << mateinfo << '\t' << matepos << '\t' << dist
                        << '\t' << seq << '\t' << qual << endl;*/
                    fprintf( outsam, "%s\t%s\t%s\t%d\t%d\t%s\t*\t0\t0\t%s\t%s\tXG:Z:CT\n",
								seqName+IDstart, flag, chr, pos, score, cigar.c_str(),
								seq.c_str(), qual.c_str() );
                }

                ++ cntGT[tn];
            }
            free( seqName );
        }
        if( CG2TG.eof() )break;
	}
    unsigned int GTall = 0;
    for( unsigned int i=0; i!=thread; ++i )
        GTall += cntGT[i];

	cout << "WATSON\t" << GTall << '\n'
		 << "CRICK\t"  << CAall << '\n';

	CG2CA.close();
	CG2TG.close();
    fclose( outsam );

	delete [] hits;
	delete [] R1;
    delete [] cntCA;
    delete [] cntGT;

	return 0;
}

