#include <iostream>
#include <string.h>
#include "SmithWatermanGotoh.h"
#include "BandedSmithWaterman.h"
#include "parameter_parser.h"
#include "references.h"

using namespace std;

int main(int argc, char* argv[]) {
  Parameters param;
  ParseArgumentsOrDie(argc, argv, &param);

  References refs;
  refs.LoadReferences(param.fasta.c_str());

	//=====================================================
	// defind the hash region
	// first.first:   reference begin
	// first.second:  reference end
	// second.first:  query begin
	// second.second: query end
	//=====================================================
	
	pair< pair<unsigned int, unsigned int>, pair<unsigned int, unsigned int> > hr;
	hr.first.first   = 5;
	hr.first.second  = 13;
	hr.second.first  = 0;
	hr.second.second = 8;

	//=====================================================

	// for 76 bp reads, we expect as much as 12 mismatches - however this does not
	// translate to a bandwidth of 12 * 2 + 1 since most of these will be
	// substitution errors
	const unsigned char bandwidth = 11;

	// initialize
	const char* pReference = "ATGGCGGGGATCGGGACACTCGCCGGTGCGGGTACCCTA";
	const char* pQuery     =      "GGGGATCGGGACACTCGCTCTCCGGTGCGGGTA";
	
	const unsigned int referenceLen = strlen(pReference);
	const unsigned int queryLen     = strlen(pQuery);

	// ==============================================================================================
	// benchmarking reference on koi.bc.edu when NUM_ITERATIONS = 38000 on 76 bp read (1 try):
	// CPU time: 23.920 s, wall time: 24.012 s (1582.6 alignments/s)
	// ==============================================================================================
	//const unsigned int NUM_ITERATIONS = 38000;
	//unsigned int NUM_ITERATIONS = 1;

	// create a new Smith-Waterman alignment object
	CSmithWatermanGotoh sw(10.0f, -9.0f, 15.0f, 6.66f);
	CBandedSmithWaterman bsw(10.0f, -9.0f, 15.0f, 6.66f, bandwidth);

	// start timing the algorithm
	//CBenchmark bench;
	//bench.Start();

	// perform NUM_ITERATIONS alignments
	//Alignment bswAl;
	//Alignment swAl;
	//   referenceBegin, referenceEnd
	unsigned int referenceSW, referenceBSW;
	string cigarSW, cigarBSW;
	//for(unsigned int i = 0; i < NUM_ITERATIONS; i++) {
	  sw.Align(referenceSW, cigarSW, pReference, referenceLen, pQuery, queryLen);
	  bsw.Align(referenceBSW, cigarBSW, pReference, referenceLen, pQuery, queryLen, hr);
	//}

	// stop timing the algorithm
	//bench.Stop();
	
	// calculate the alignments per second
	//double elapsedWallTime = bench.GetElapsedWallTime();
	//double alignmentsPerSecond = (double)NUM_ITERATIONS / elapsedWallTime;
	
	// show our results
	//printf("%d\t%d\n", al.ReferenceBegin,al.QueryBegin);

	printf("Smith-Waterman\n");
	printf("reference:    %s %3u\n", cigarSW.c_str(), referenceSW);
	printf("Banded Smith-Waterman\n");
	printf("reference:    %s %3u\n", cigarBSW.c_str(), referenceBSW);
	/*  
	printf("Smith-Waterman\n");
	printf("reference:    %s %3u %3u\n", swAl.Reference.CData(), swAl.ReferenceBegin, swAl.ReferenceEnd);
        printf("query:        %s %3u %3u\n", swAl.Query.CData(), swAl.QueryBegin, swAl.QueryEnd);
        printf("mismatches:   %u\n", swAl.NumMismatches);
	printf("\n");	
	printf("Banded Smith-Waterman\n");
	printf("reference:    %s %3u %3u\n", bswAl.Reference.CData(), bswAl.ReferenceBegin, bswAl.ReferenceEnd);
        printf("query:        %s %3u %3u\n", bswAl.Query.CData(), bswAl.QueryBegin, bswAl.QueryEnd);
        printf("mismatches:   %u\n", bswAl.NumMismatches);
	*/
	//printf("alignments/s: %.1f\n\n", alignmentsPerSecond);

	//bench.DisplayTime("BandedSmithWaterman");

	return 0;
}
