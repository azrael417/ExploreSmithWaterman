#include <iostream>
#include <string.h>
#include <time.h>
#include "SmithWatermanGotoh.h"
//#include "BandedSmithWaterman.h"
#include "parameter_parser.h"
#include "references.h"
#include "fastq_reader.h"

//include kokkos
#include "types.h"

using namespace std;

inline void PrintAlignment(const string& name, const string& seq, const Alignment& al) {
 
 cout << "read_name: " << name
      << '\n'
      << "read_seq: " << seq
      << '\n'
      << "max score: " << al.sw_score 
      << ", begin_ref: " << al.ref_begin << ", begin_read: " << al.query_begin
      << ", end_ref: " << al.ref_end << ", query_end: " << al.query_end
      << '\n'
      << "cigar: " << al.cigarCount << al.cigarChar
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

  int length = 0, num_total_aligns;
  clock_t start, end;
  unsigned int max_sequence_length=0, max_reference_length=0;

  // determine max sequence length and counts
  int sequence_count=0;
  {
    FastqReader fastqtmp;
    string readname, sequence, qual;
    fastqtmp.Open(param.fastq.c_str());
    while (fastqtmp.LoadNextRead(&readname, &sequence, &qual)) {
      const unsigned int sequence_length = sequence.size();
      if (sequence_length > max_sequence_length) {
	      max_sequence_length = sequence_length;
      }
      sequence_count++;
    }
  }
  cout << "Number of sequences: " << sequence_count << endl;
  cout << "Max sequence length: " << max_sequence_length << endl;

  // max reference length
  max_reference_length = refs.GetMaxSequenceLength();
  cout << "Number of references: " << refs_count << endl;
  cout << "Max reference length: " << max_reference_length << "\n";
  
  //handle negative num batches
  if(param.num_batches <= 0) param.num_batches = sequence_count;
  if(param.num_tasks <= 0) param.num_tasks = 1;
  
  //scope ensures proper deletion of sw object
  {
    CSmithWatermanGotoh sw(param.match, 0-param.mismatch, param.open_gap, param.extend_gap, 
			   max_reference_length, max_sequence_length, static_cast<unsigned int>(param.num_tasks));

    //allocate alignment buffer
    View2D<Alignment> alignments = View2D<Alignment>("alignments", refs_count, param.batchsize);

    //start clock
    start = clock();

    //do batches
    num_total_aligns = 0;
    int num_batches=0;
    while ( fastq.LoadNextBatch(param.batchsize) & (num_batches<param.num_batches) ) {
      
      //perform alignment
      auto sequences = fastq.GetSequences();
      auto sequences_length = fastq.GetSequenceLengths();
      
      //define partition launch kernel
      auto master = [&]( int partition_id, int num_partitions ) {
        
        //perform tasking over batch dim:
        int jblock = param.batchsize / num_partitions;
        int jb = jblock * partition_id;
        int je = std::min(jb+jblock, param.batchsize);

        //do the loop
        for (int j=jb; j < je; ++j){
          for (int i=0; i < refs_count; ++i){
            int tmplen;
            auto tmpread = refs.GetReferenceSequence(i, &tmplen);
            sw.Align(partition_id, alignments(i, j), tmpread, tmplen, Kokkos::subview(sequences, j, Kokkos::ALL), sequences_length(j));
          }
        }
      };
      
      //dispatch
      Kokkos::OpenMP::partition_master( master, param.num_tasks );
      
      //print alignment
      for (int j = 0; j < param.batchsize; ++j){
        string tmpseq, tmpread;
        ViewToString(tmpseq, fastq.GetSequence(j));
        ViewToString(tmpread, fastq.GetRead(j));
        
        for (int i = 0; i < refs_count; ++i) {
          PrintAlignment(tmpread, tmpseq, alignments(i,j));
          num_total_aligns++;
        }
      }
      num_batches++;
    }
    //do remainder if necessary
    int remaindercount = fastq.GetReadsize();
    if(remaindercount != param.batchsize){
      //perform alignment
      for (int j = 0; j < remaindercount; ++j){
        for (int i = 0; i < refs_count; ++i) {
          auto tmpread = refs.GetReferenceSequence(i, &length);
          sw.Align(0, alignments(i, j), tmpread, length, fastq.GetSequence(j), fastq.GetSequenceLength(j));
        }
      }
      //print alignment
      for (int j = 0; j < remaindercount; ++j){
        string tmpseq, tmpread;
        ViewToString(tmpseq, fastq.GetSequence(j));
        ViewToString(tmpread, fastq.GetRead(j));
  
        for (int i = 0; i < refs_count; ++i) {
          PrintAlignment(tmpread, tmpseq, alignments(i,j));
          num_total_aligns++;
        }
      }
    }
    
    //end clock
    end = clock();
  }
  
  float cpu_time = (static_cast<float>(end - start)) / static_cast<float>(CLOCKS_PER_SEC);
  fprintf(stdout, "CPU time: %f seconds\n", cpu_time);
  fprintf(stdout, "Number of total alignments: %i\n", num_total_aligns);
  fprintf(stdout, "Alignments per second: %f\n", (float)num_total_aligns / cpu_time);

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
