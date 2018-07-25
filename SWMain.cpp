#include <iostream>
#include <string.h>
#include <time.h>
#include "SmithWatermanGotoh.h"
//#include "BandedSmithWaterman.h"
#include "parameter_parser.h"
#include "references.h"
#include "fastq_reader.h"

//include kokkos


using namespace std;

inline void PrintAlignment(const string& name, const string& seq,
                    const string& cigar, const Alignment& al) {
 cout << "read_name: " << name
      << endl
      << "read_seq: " << seq
      << endl
      << "max score: " << al.sw_score 
      << ", begin_ref: " << al.ref_begin << ", begin_read: " << al.query_begin
      << ", end_ref: " << al.ref_end << ", query_end: " << al.query_end
      << endl
      << "cigar: " << cigar
      << endl;
}

int main(int argc, char* argv[]) {
  
  //initialize kokkos
  Kokkos::initialize();
  
  Parameters param;
  ParseArgumentsOrDie(argc, argv, &param);

  References refs;
  refs.LoadReferences(param.fasta.c_str());
  const int refs_count = refs.GetReferenceCount();

  FastqReader fastq;
  fastq.Open(param.fastq.c_str());

  string readname, *sequences, *quals, cigarSW;
  int length = 0, readsize;
  Alignment alignment;
  clock_t start, end;
  unsigned int max_sequence_length, max_reference_length;

  // determine max sequence length
  {
    FastqReader fastqtmp;
    fastqtmp.Open(param.fastq.c_str());
    while (fastqtmp.LoadNextRead(&readname, &sequence, &qual)) {
      const unsigned int sequence_length = sequence.size();
      if (sequence_length > max_sequence_length) {
	max_sequence_length = sequence_length;
      }
    }
  }
  
  cout << "max sequence length: " << max_sequence_length << "\n";

  // max reference length
  for (int i = 0; i < refs_count; ++i) {
    const char* pReference = refs.GetReferenceSequence(i, &length);
    if (length > max_reference_length) {
      max_reference_length = length;
    }
  }

  cout << "max reference length: " << max_reference_length << "\n";
  
  //scope ensures proper deletion of sw object
  {
    CSmithWatermanGotoh sw(param.match, 0-param.mismatch, param.open_gap, param.extend_gap, 
			   max_reference_length, max_sequence_length);

    //start clock
    start = clock();

    //do batches
    while (fastq.LoadNextBatch(&readname, &sequences, &quals, &readsize, param.batchsize)) {
      for (int j = 0; j < param.batchsize; ++j){
        for (int i = 0; i < refs_count; ++i) {
          const char* pReference = refs.GetReferenceSequence(i, &length);
          sw.Align(&alignment, cigarSW, pReference, length, sequences[j].c_str(), sequences[j].size());
          PrintAlignment(readname, sequences[j], cigarSW, alignment);
        }
      }
    }
    //do remainder
    for (int j = 0; j < readsize; ++j){
      for (int i = 0; i < refs_count; ++i) {
        const char* pReference = refs.GetReferenceSequence(i, &length);
        sw.Align(&alignment, cigarSW, pReference, length, sequences[j].c_str(), sequences[j].size());
        PrintAlignment(readname, sequences[j], cigarSW, alignment);
      }
    }
    
    //end clock
    end = clock();
  }
  
  float cpu_time = (static_cast<float>(end - start)) / static_cast<float>(CLOCKS_PER_SEC);
  fprintf(stdout, "CPU time: %f seconds\n", cpu_time);

	//=====================================================
	// defind the hash region
	// first.first:   reference begin
	// first.second:  reference end
	// second.first:  query begin
	// second.second: query end
	//=====================================================
	
	//pair< pair<unsigned int, unsigned int>, pair<unsigned int, unsigned int> > hr;
	//hr.first.first   = 5;
	//hr.first.second  = 13;
	//hr.second.first  = 0;
	//hr.second.second = 8;

	// for 76 bp reads, we expect as much as 12 mismatches - however this does not
	// translate to a bandwidth of 12 * 2 + 1 since most of these will be
	// substitution errors
	//const unsigned char bandwidth = 11;

	// initialize
	//const char* pReference = "ATGGCGGGGATCGGGACACTCGCCGGTGCGGGTACCCTA";
	//const char* pQuery     =      "GGGGATCGGGACACTCGCTCTCCGGTGCGGGTA";
	
	//const unsigned int referenceLen = strlen(pReference);
	//const unsigned int queryLen     = strlen(pQuery);

	// create a new Smith-Waterman alignment object
	//CSmithWatermanGotoh sw(10.0f, -9.0f, 15.0f, 6.66f);
	//CBandedSmithWaterman bsw(10.0f, -9.0f, 15.0f, 6.66f, bandwidth);

	//   referenceBegin, referenceEnd
	//unsigned int referenceSW, referenceBSW;
	//string cigarSW, cigarBSW;
	//sw.Align(referenceSW, cigarSW, pReference, referenceLen, pQuery, queryLen);
	//bsw.Align(referenceBSW, cigarBSW, pReference, referenceLen, pQuery, queryLen, hr);

	//printf("Smith-Waterman\n");
	//printf("reference:    %s %3u\n", cigarSW.c_str(), referenceSW);
	//printf("Banded Smith-Waterman\n");
	//printf("reference:    %s %3u\n", cigarBSW.c_str(), referenceBSW);
  
  //finalize kokkos
  Kokkos::finalize();
  
	return 0;
}
