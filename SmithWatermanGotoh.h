#ifndef SMITHWATERMANGOTOH_H_
#define SMITHWATERMANGOTOH_H_

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <memory>
#include "Mosaik.h"
#include <sstream>
#include <string>

//kokkosification
#include <kokkos/core/src/Kokkos_Core.hpp>

using namespace std;

#define MOSAIK_NUM_NUCLEOTIDES 26
#define GAP '-'

struct Alignment {
  uint64_t ref_begin;
  uint64_t ref_end;
  uint64_t query_begin;
  uint64_t query_end;
  float sw_score;
};

class CSmithWatermanGotoh {
public:
	// constructor
        CSmithWatermanGotoh(float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty, 
		            unsigned int max_reference_length, unsigned int max_sequence_length);
	// destructor
	~CSmithWatermanGotoh(void);
	// aligns the query sequence to the reference using the Smith Waterman Gotoh algorithm
	void Align(Alignment* alignment, string& cigarAl, const char* s1, const unsigned int s1Length, const char* s2, const unsigned int& s2Length);
	// enables homo-polymer scoring
	void EnableHomoPolymerGapPenalty(float hpGapOpenPenalty);
private:
	// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
	void CreateScoringMatrix(void);
	// deletes the simple scoring matrix
	void DestroyScoringMatrix(void);
	// corrects the homopolymer gap order for forward alignments
	void CorrectHomopolymerGapOrder(const unsigned int numBases, const unsigned int numMismatches);
	// returns the maximum floating point number
	static inline float MaxFloats(const float& a, const float& b, const float& c);
	// our simple scoring matrix
	Kokkos::View<float**> mScoringMatrix; //size [MOSAIK_NUM_NUCLEOTIDES][MOSAIK_NUM_NUCLEOTIDES];
	// keep track of maximum initialized sizes
	uint64_t mCurrentMatrixSize;
	uint64_t mCurrentAnchorSize;
	uint64_t mCurrentQuerySize;
	uint64_t mCurrentAQSumSize;
	// define our traceback directions
	// N.B. This used to be defined as an enum, but gcc doesn't like being told
	// which storage class to use
	const static char Directions_STOP;
	const static char Directions_LEFT;
	const static char Directions_DIAGONAL;
	const static char Directions_UP;
	// define scoring constants
	const float mMatchScore;
	const float mMismatchScore;
	const float mGapOpenPenalty;
	const float mGapExtendPenalty;
	// store the backtrace pointers
	char* mPointers;
	// store the vertical gap sizes - assuming gaps are not longer than 32768 bases long
	short* mSizesOfVerticalGaps;
	// store the horizontal gap sizes - assuming gaps are not longer than 32768 bases long
	short* mSizesOfHorizontalGaps;	
	// score if xi aligns to a gap after yi
	float* mQueryGapScores;
	// best score of alignment x1...xi to y1...yi
	float* mBestScores;
	// our reversed alignment
	char* mReversedAnchor;
	char* mReversedQuery;
	// define static constants
	static const float FLOAT_NEGATIVE_INFINITY;
	// toggles the use of the homo-polymer gap open penalty
	bool mUseHomoPolymerGapOpenPenalty;
	// specifies the homo-polymer gap open penalty
	float mHomoPolymerGapOpenPenalty;
};

// returns the maximum floating point number
inline float CSmithWatermanGotoh::MaxFloats(const float& a, const float& b, const float& c) {
	float max = 0.0f;
	if(a > max) max = a;
	if(b > max) max = b;
	if(c > max) max = c;
	return max;
}

#endif // SMITHWATERMANGOTOH_H_
