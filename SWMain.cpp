#include <iostream>
#include <string.h>
#include <time.h>
#include "SmithWatermanGotoh.h"
#include "BandedSmithWaterman.h"
#include "parameter_parser.h"
#include "references.h"
#include "fastq_reader.h"

using namespace std;

int main(int argc, char* argv[]) {
  Parameters param;
  ParseArgumentsOrDie(argc, argv, &param);

  References refs;
  refs.LoadReferences(param.fasta.c_str());
  const int refs_count = refs.GetReferenceCount();

  FastqReader fastq;
  fastq.Open(param.fastq.c_str());

  string readname, sequence, qual, cigarSW;
  int length = 0;
  unsigned int referenceSW;
  CSmithWatermanGotoh sw(10.0f, -9.0f, 15.0f, 6.66f);
  clock_t start = clock();

  while (fastq.LoadNextRead(&readname, &sequence, &qual)) {
    for (int i = 0; i < refs_count; ++i) {
      const char* pReference = refs.GetReferenceSequence(i, &length);
      sw.Align(referenceSW, cigarSW, pReference, length, sequence.c_str(), sequence.size());

    }
  }

  clock_t end = clock();
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
	return 0;
}
