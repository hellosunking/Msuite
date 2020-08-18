#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <omp.h>
#include "common.h"

using namespace std;

/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * This program is part of the Msuite package, adapted from Ktrim
 * Date: Jan 2020
 *
 * The conversion log's rule:
 *   read 1:
 *	   @LINE_NUMBER+S|xx;xx$ => for reads with endC and C>T changes, 'S' is its quality score
 *	   @LINE_NUMBER+xx;xx$	 => for reads without endC
 *	   @LINE_NUMBER+$		 => for reads without endC and conversions
 *   read 2:
 *	   @S|xx;xx$ => for reads with frontG and G>A changes, 'S' is its quality score
 *	   @xx;xx$	 => for reads without frontG
 *	   @$		 => for reads without frontG and conversions
 *
 * TODO: Use LARGE buffer to store the modified reads per batch then write
 *       it to the output files after the multi-thread trimming/conversion
**/

/*
 * use dynamic max_mismatch as the covered size can range from 3 to a large number such as 50,
 * so use static values (e.g., 4) is not good
*/
bool check_mismatch_dynamic_SE( const string & s, unsigned int pos, const adapter_info* ai ) {
	register unsigned int mis=0;
	register unsigned int i, len;
	len = s.length() - pos;
	if( len > ai->adapter_len )
	  len = ai->adapter_len;
	register unsigned int max_mismatch_dynamic = len >> 2;
	if( (max_mismatch_dynamic<<2) != len )
	  ++ max_mismatch_dynamic;
	const char * p = s.c_str();
	for( i=0; i!=len; ++i ) {
		if( p[pos+i] != ai->adapter_r1[i] ) {
			++ mis;
			if( mis > max_mismatch_dynamic )
			  return false;
		}
	}

	return true;
}


int main( int argc, char *argv[] ) {
	if( argc < 5 ) {
		cerr << "\nUsage: " << argv[0] << " <r1.fq> <r2.fq=placeholder> <cycle> <out.prefix> "
			 << "[mode=0|3|4] [thread=1] [min.length=36] [min.quality=53] [libraryKit=illumina]\n\n"

			 << "This program is part of Msuite and is designed to do fastq statistics, quality-trimming,\n"
			 << "adapter-trimming and C->T/G->A conversions for Paired-End reads generated by illumina sequencers.\n\n"

			 << "Run modes:\n"
			 << "  0. Do not perform any conversions. Usable for normal DNA/RNA-seq alignment.\n"
			 << "  3. Do C>T (G>A) for all C (G) sites in read1 (read2). Usable for 3-letter BS-seq data.\n"
			 << "  4. Do C>T (G>A) for C (G) in CpG sites in read1 (read2). Usable for 4-letter TAPS data.\n\n"

			 << "Default parameters:\n"
			 << "  mode: 0\n"
			 << "  min.length: 36\n"
			 << "  min.quality: 53 (33+20 for phred33('!') scoring system)\n\n"

			 << "Other commonly used Phred scoring systems are 35('#') and 64('@').\n"
			 << "You may need to set this parameter manually based on your data.\n\n"
			 << "By default, illumina adapters are used:\n\n"
			 << "  5': ACACTCTTTCCCTACACGACGCTCTTCCGATCT\n"
			 << "  3': AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG\n\n"
			 << "However, Nextera adapters are also supported:\n\n"\
			 << "  5': CTGTCTCTTATACACATCT\n"
			 << "  3': CTGTCTCTTATACACATCT\n\n"
			 << "Sequencing model (built-in):\n\n"
			 << "              *read1 -->\n"
			 << "  5'adapter - NNNNNNNNNNsequenceNNNNNNNNNN - 3'adapter\n"
			 << "                                <-- read2*\n\n";

		return 2;
	}
	// for SE, the r2.fq will be ignored (could be set to /dev/null ro anything else)
	// but to provide an identical interface to preprocessor.pe, its position is kept

	unsigned int  mode = 0;
	unsigned int  thread = 1;
	unsigned int  min_length = 36;
	unsigned char quality = 53;
	const char *libraryKit;
	const adapter_info* ai;

	unsigned int  cycle = atoi( argv[3] );
	if( cycle == 0 ) {
		cerr << "Error: Unacceptable cycle!\n";
		exit(1);
	}
	if( argc > 5 ) {
		mode = atoi( argv[5] );
		if( argc > 6) {
			thread = atoi( argv[6] );
			if( argc > 7 ) {
				min_length = atoi( argv[7] );
				if( argc > 8 ) {
					quality = (unsigned char) atoi( argv[8] );
					if( argc > 9 ) {
						libraryKit = argv[9];
					} else {
						libraryKit = "illumina";
					}
				}
			}
		}
	}
	if( mode != 0 && mode !=3 && mode != 4 ) {
		cerr << "Error: invalid run mode! Must be 0, 3, or 4!\n";
		return 100;
	}
	if( thread == 0 ) {
		cerr << "Warning: thread is set to 0! I will use all threads instead.\n";
		thread = omp_get_max_threads();
	}
	if( min_length == 0 ) {
		cerr << "Error: invalid min_length! Must be a positive number!\n";
		return 101;
	}
	if( quality == 0 ) {
		cerr << "Error: invalid quality! Must be a positive number!\n";
		return 102;
	}
	if( strcmp(libraryKit, "illumina")==0 || strcmp(libraryKit, "Illumina")==0 ) {
		ai = &illumina_adapter;
	} else if ( strcmp(libraryKit, "nextera")==0 || strcmp(libraryKit, "Nextera")==0 ) {
		ai = &nextera_adapter;
	} else if ( strcmp(libraryKit, "bgi")==0 || strcmp(libraryKit, "Nextera")==0 ) {
		ai = &bgi_adapter;
	} else {
		cerr << "Error: invalid library kit! Currently only supports illumina and nextera!\n";
		return 103;
	}

	ifstream fq1( argv[1] );
	if( fq1.fail() ) {
		cerr << "Error: open fastq file failed!\n";
		fq1.close();
		return 2;
	}
	string base = argv[4];
	base += ".R1.fq";
	ofstream fout1( base.c_str() );
	if( fout1.fail() ) {
		cerr << "Error: write file failed!\n";
		fq1.close();
		return 3;
	}
	// set HEX format number output
	fout1.setf(ios::hex, ios::basefield);

	string *id1   = new string [READS_PER_BATCH];
	string *seq1  = new string [READS_PER_BATCH];
	string *qual1 = new string [READS_PER_BATCH];
	string unk;

	register unsigned int line = 1;
	int *dropped	  = new int [thread];
	int *real_adapter = new int [thread];
	int *tail_adapter = new int [thread];

	fastqstat * AllR1stat = new fastqstat[ cycle ];
	memset( AllR1stat, 0, cycle*sizeof(fastqstat) );

	// buffer for storing the modified reads per thread
	char ** buffer1 = new char * [thread];
	int  * b1stored = new int	 [thread];
	fastqstat **R1stat = new fastqstat * [thread];
	for(unsigned int i=0; i!=thread; ++i) {
		buffer1[i] = new char[ BUFFER_SIZE_PER_BATCH_READ ];
		R1stat[i]  = new fastqstat [ cycle ];
		dropped[i] = 0;
		real_adapter[i] = 0;
		tail_adapter[i] = 0;
	}

	cerr << "Loading files ...\n";
	while( true ) {
		// get fastq reads
		unsigned int loaded = 0;
		while( true ) {
			getline( fq1, id1  [ loaded ] );
			if( fq1.eof() )break;
			getline( fq1, seq1 [ loaded ] );
			getline( fq1, unk );
			getline( fq1, qual1[ loaded ] );

			if( id1[loaded].size() > cycle ) {
				id1[loaded].resize( cycle );
				qual1[loaded].resize( cycle );
			}

			++ loaded;
			if( loaded == READS_PER_BATCH )
				break;
		}
		//cerr << "loaded: " << loaded << '\n';
		if( loaded == 0 )
			break;

		// start parallalization
		omp_set_num_threads( thread );
		#pragma omp parallel
		{
			unsigned int tn = omp_get_thread_num();
			unsigned int start = loaded * tn / thread;
			unsigned int end   = loaded * (tn+1) / thread;

			// normalization
			b1stored[tn] = 0;
			memset( R1stat[tn], 0, cycle*sizeof(fastqstat) );
		
			string conversionLog;
			register int i, j;
			register unsigned int last_seed;
			vector<unsigned int> seed;
			vector<unsigned int> :: iterator it;
			const char *p, *q;
			char *conversion = new char [MAX_CONVERSION];
			char numstr[10]; // enough to hold all numbers up to 99,999,999 plus ':'

			for( unsigned int ii=start; ii!=end; ++ii ) {
				// fqstatistics
				p = seq1[ii].c_str();
				j = seq1[ii].size();
				if( j > cycle )
					j = cycle;
				for( i=0; i!=j; ++i ) {
					switch ( p[i] ) {
						case 'a':
						case 'A': R1stat[tn][i].A ++; break;
						case 'c':
						case 'C': R1stat[tn][i].C ++; break;
						case 'g':
						case 'G': R1stat[tn][i].G ++; break;
						case 't':
						case 'T': R1stat[tn][i].T ++; break;
						default : R1stat[tn][i].N ++; break;
					}
				}

				// quality control
				p = qual1[ii].c_str();
				for( i=qual1[ii].length()-1; i; --i ) {
					if( p[i] >= quality ) break;
				}
				++ i;
				if( i < min_length ) { // not long enough
					++ dropped[ tn ];
					continue;
				}
				seq1[ii].resize(  i );
				qual1[ii].resize( i );

				// looking for seed target, 1 mismatch is allowed for these 2 seeds
				// which means seq1 and seq2 at least should take 1 perfect seed match
				seed.clear();
				for( i=0; (i=seq1[ii].find(ai->adapter_index, i)) != string::npos; ++i )
					seed.push_back( i );

				last_seed = impossible_seed;	// a position which cannot be in seed
				for( it=seed.begin(); it!=seed.end(); ++it ) {
					if( check_mismatch_dynamic_SE(seq1[ii], *it, ai) )
						break;
				}
				if( it != seed.end() ) {	// adapter found
					++ real_adapter[tn];
					if( *it >= min_length )	{
						seq1[ii].resize(  *it );
						qual1[ii].resize( *it );
					} else {	// drop this read as its length is not enough
						++ dropped[tn];
						continue;
					}
				} else {	// seed not found, now check the tail 2 or 1, if perfect match, drop these 2
					i = seq1[ii].length() - 2;
					p = seq1[ii].c_str();
					if( p[i]==ai->adapter_r1[0] && p[i+1]==ai->adapter_r1[1] ) {
						if( i < min_length ) {
							++ dropped[tn];
							continue;
						}
						seq1[ii].resize(  i );
						qual1[ii].resize( i );

						++ tail_adapter[tn];
/*//maybe it is not that good to check tail-1?
					} else {	// tail 2 is not good, check tail 1
						++ i;
						if( p[i] == ai->adapter_r1[0] ) {
							if( i < min_length ) {
								++ dropped[tn];
								continue;
							}
							seq1[ii].resize(  i );
							qual1[ii].resize( i );

							++ tail_adapter[tn];
						}
*/
					}
				}

				//check if there is any white space in the IDs; if so, remove all the data after the whitespace
				j = id1[ii].size();
				p = id1[ii].c_str();
				for( i=1; i!=j; ++i ) {
					if( p[i]==' ' || p[i]=='\t' ) {	// white space, then trim ID
						id1[ii].resize( i );
						break;
					}
				}

				// do C->T and G->A conversion
				if( mode == 0 ) {	// no need to do conversion
					b1stored[tn] += sprintf( buffer1[tn]+b1stored[tn], "%s\n%s\n+\n%s\n",
								id1[ii].c_str(), seq1[ii].c_str(), qual1[ii].c_str() );
				} else if( mode == 3 ) {	// in this implementation, id1 and id2 are different!!!
					// modify id1 to add line number (to facilitate removing ambigous step)
					// in mode 3, there is NO endC and frontG issues
					id1[ii][0] = CONVERSION_LOG_END;
					j = seq1[ii].size();	// seq1 and seq2 are of the same size
					conversionLog = LINE_NUMBER_SEPARATOR;
					for( i=0; i!=j; ++i ) {
						if( seq1[ii][i] == 'C' ) {
							seq1[ii][i] = 'T';
							sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
							conversionLog += numstr;
						}
					}
					if( conversionLog.back() == CONVERSION_LOG_SEPARATOR )
						conversionLog.pop_back();
					/*fout1 << NORMAL_SEQNAME_START << line << conversionLog << id1 << '\n'
							<< seq1 << "\n+\n" << qual1 << '\n';*/
					b1stored[tn] += sprintf( buffer1[tn]+b1stored[tn], "%c%x%s%s\n%s\n+\n%s\n",
												NORMAL_SEQNAME_START, line+ii, conversionLog.c_str(),
												id1[ii].c_str(), seq1[ii].c_str(), qual1[ii].c_str() );
				} else if ( mode == 4 ) {	// this is the major task for EMaligner
					// modify id1 to add line number (to facilitate the removing ambigous step)
					// check seq1 for C>T conversion
					id1[ii][0] = CONVERSION_LOG_END;
					conversionLog = LINE_NUMBER_SEPARATOR;
					j = seq1[ii].size()-1;
					if( seq1[ii].back() == 'C' ) { //ther is a 'C' and the end, discard it (but record its Quality score);
						//otherwise it may introduce a mismatch in alignment
						conversionLog += qual1[ii].back();
						conversionLog += KEEP_QUAL_MARKER;
						seq1[ii].pop_back();
						qual1[ii].pop_back();
					}
					for( i=0; i!=j; ++i ) {
						if( seq1[ii][i]=='C' && seq1[ii][i+1]=='G' ) {
							seq1[ii][i] = 'T';
							sprintf( numstr, "%x%c", i, CONVERSION_LOG_SEPARATOR );
							conversionLog += numstr;
						}
					}
					if( conversionLog.back() == CONVERSION_LOG_SEPARATOR )
						conversionLog.pop_back();
					/*fout1 << NORMAL_SEQNAME_START << line << conversionLog << id1 << '\n'
							<< seq1 << "\n+\n" << qual1 << '\n';*/
					b1stored[tn] += sprintf( buffer1[tn]+b1stored[tn], "%c%x%s%s\n%s\n+\n%s\n",
								NORMAL_SEQNAME_START, line+ii, conversionLog.c_str(),
								id1[ii].c_str(), seq1[ii].c_str(), qual1[ii].c_str() );
					// format for ID1:
					// if there is a C at the end
					//	@line_number '+' x| C1;C2;C3$ raw_seq_name
					//	the | is the marker for the existence of tail 'C' and 'x' is its quality score
					//	if exists, | is ALWAYS two bytes after '+' (use this to test its existence)
					// if there is No C at the end
					//	@line_number '+' C1&C2&C3$ raw_seq_name
					//
					// All the numbers in line_number and C1,C2,C3... are HEX
				}
			}
		}	// parallel body
		// write output and update fastq statistics
		for( unsigned int i=0; i!=thread; ++i ) {
			fout1 << buffer1[i];

			for( unsigned int j=0; j!=cycle; ++j ) {
				AllR1stat[j].A += R1stat[i][j].A;
				AllR1stat[j].C += R1stat[i][j].C;
				AllR1stat[j].G += R1stat[i][j].G;
				AllR1stat[j].T += R1stat[i][j].T;
				AllR1stat[j].N += R1stat[i][j].N;
			}
		}
		line += loaded;
		cerr << '\r' << line-1 << " reads loaded";

		if( fq1.eof() )break;
	}
	fq1.close();
	fout1.close();

	cerr << "\rDone: " << line-1 << " lines processed.\n";

	// write trim.log
	ofstream fout( "Msuite.trim.log" );
	if( fout.fail() ) { 
		cerr << "Error: cannot write log file!\n";
		return 4;
	}
	int dropped_all=0, real_all=0, tail_all=0;
	for( unsigned int i=0; i!=thread; ++i ) {
		dropped_all += dropped[i];
		real_all += real_adapter[i];
		tail_all += tail_adapter[i];
	}
	fout << line << '\n'	// total
		 << "Dropped : " << dropped_all << '\n'
		 << "Aadaptor: " << real_all	<< '\n'
		 << "Tail Hit: " << tail_all	<< '\n';
	fout.close();

	// write fqstatistics
	fout.open( "R1.fqstat" );
	if( fout.fail() ) {
		cerr << "Error: cannot write R1.fqstat file!\n";
		return 5;
	}
	fout << "Cycle\tA\tC\tG\tT\tN\n";
	for( unsigned int j=0; j!=cycle; ++j ) {
		fout << j+1 << '\t' << AllR1stat[j].A << '\t' << AllR1stat[j].C << '\t'<< AllR1stat[j].G
				<< '\t' << AllR1stat[j].T << '\t'<< AllR1stat[j].N << '\n';
	}
	fout.close();

	//free memory
	for(unsigned int i=0; i!=thread; ++i) {
		delete buffer1[i];
		delete R1stat[i];
	}
	delete [] buffer1;
	delete [] R1stat;
	delete [] AllR1stat;

	return 0;
}

