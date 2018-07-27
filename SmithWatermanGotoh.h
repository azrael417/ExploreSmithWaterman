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
#include "types.h"

using namespace std;

#define MOSAIK_NUM_NUCLEOTIDES 26
#define GAP '-'

#define DIRECTIONS_STOP 0
#define DIRECTIONS_LEFT 1
#define DIRECTIONS_DIAGONAL 2
#define DIRECTIONS_UP 3

#define FLOAT_NEGATIVE_INFINITY -1e+30

struct Alignment {
  uint64_t ref_begin;
  uint64_t ref_end;
  uint64_t query_begin;
  uint64_t query_end;
  uint64_t cigarCount;
  char cigarChar;
  float sw_score;
};

class CSmithWatermanGotoh {
public:
	// constructor
  CSmithWatermanGotoh(float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty, 
		            unsigned int max_reference_length, unsigned int max_sequence_length, unsigned int num_tasks);
	// destructor
	~CSmithWatermanGotoh(void);
	// aligns the query sequence to the reference using the Smith Waterman Gotoh algorithm
	KOKKOS_INLINE_FUNCTION void Align(const unsigned int teamid, Alignment& alignment, View1D<char> s1, const unsigned int& s1Length, View1D<char> s2, const unsigned int& s2Length) const;
	// enables homo-polymer scoring
	void EnableHomoPolymerGapPenalty(float hpGapOpenPenalty);
private:
	// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
	void CreateScoringMatrix(void);
  //allocate temp buffers
  void InitArrays();
  void InitArrays(unsigned int max_reference_length, unsigned int max_sequence_length, unsigned int num_tasks=0);
	// corrects the homopolymer gap order for forward alignments
	KOKKOS_INLINE_FUNCTION void CorrectHomopolymerGapOrder(const unsigned int teamid, const unsigned int numBases, const unsigned int numMismatches) const;
	// returns the maximum floating point number
	static inline float MaxFloats(const float& a, const float& b, const float& c);
	// our simple scoring matrix
	View2D<float> mScoringMatrix; //size [MOSAIK_NUM_NUCLEOTIDES][MOSAIK_NUM_NUCLEOTIDES];
	// keep track of maximum initialized sizes
	uint64_t mCurrentMatrixSize;
	uint64_t mCurrentAnchorSize;
	uint64_t mCurrentQuerySize;
	uint64_t mCurrentAQSumSize;
  uint64_t mNumTasks;
	// define our traceback directions
	// N.B. This used to be defined as an enum, but gcc doesn't like being told
	// which storage class to use
	//const char Directions_STOP;
	//const char Directions_LEFT;
	//const char Directions_DIAGONAL;
	//const char Directions_UP;
	// define scoring constants
	const float mMatchScore;
	const float mMismatchScore;
	const float mGapOpenPenalty;
	const float mGapExtendPenalty;
	// store the backtrace pointers
	View2D<char, Kokkos::LayoutRight> mPointers;
	// store the vertical gap sizes - assuming gaps are not longer than 32768 bases long
	View2D<short, Kokkos::LayoutRight> mSizesOfVerticalGaps;
	// store the horizontal gap sizes - assuming gaps are not longer than 32768 bases long
	View2D<short, Kokkos::LayoutRight> mSizesOfHorizontalGaps;	
	// score if xi aligns to a gap after yi
	View2D<float, Kokkos::LayoutRight> mQueryGapScores;
	// best score of alignment x1...xi to y1...yi
	View2D<float, Kokkos::LayoutRight> mBestScores;
	// our reversed alignment
	View2D<char, Kokkos::LayoutRight> mReversedAnchor;
	View2D<char, Kokkos::LayoutRight> mReversedQuery;
	// define static constants
	//static const float FLOAT_NEGATIVE_INFINITY;
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


// aligns the query sequence to the reference using the Smith Waterman Gotoh algorithm
KOKKOS_INLINE_FUNCTION void CSmithWatermanGotoh::Align(const unsigned int teamid, Alignment& alignment, View1D<char> s1, const unsigned int& s1Length, View1D<char> s2, const unsigned int& s2Length) const {

	if((s1Length == 0) || (s2Length == 0)) {
		cout << "ERROR: Found a read with a zero length." << endl;
		exit(1);
	}

	uint64_t referenceLen      = s1Length + 1;
	uint64_t queryLen          = s2Length + 1;
	uint64_t sequenceSumLength = s1Length + s2Length;
	// set current matrix size
	uint64_t matrix_size = referenceLen * queryLen;
	//mCurrentMatrixSize = matrix_size;

	// initialize the traceback matrix to STOP
  for(uint64_t i = 0; i < queryLen; ++i) mPointers(teamid, i) = 0;
	//memset((char*)mPointers, 0, SIZEOF_CHAR * queryLen);
	for(uint64_t i = 1; i < referenceLen; ++i) {
	  mPointers(teamid, i * queryLen) = 0;
	}

	// initialize the gap matrices to 1
	//uninitialized_fill(mSizesOfVerticalGaps, mSizesOfVerticalGaps + mCurrentMatrixSize, 1);
	//uninitialized_fill(mSizesOfHorizontalGaps, mSizesOfHorizontalGaps + mCurrentMatrixSize, 1);
  for(uint64_t i=0; i < matrix_size; ++i){
    mSizesOfVerticalGaps(teamid, i) = 1;
    mSizesOfHorizontalGaps(teamid, i) = 1;
  }

	//
	// sequence lengths
	//
	//mCurrentQuerySize = s2Length;
	//mCurrentAQSumSize = sequenceSumLength;

	// initialize the gap score and score vectors
  //uninitialized_fill(mQueryGapScores, mQueryGapScores + queryLen, FLOAT_NEGATIVE_INFINITY);
  //memset((char*)mBestScores, 0, SIZEOF_FLOAT * queryLen);
  for(uint64_t i=0; i < queryLen; ++i){
    mQueryGapScores(teamid, i) = FLOAT_NEGATIVE_INFINITY;
    mBestScores(teamid, i) = 0.;
  }

	float similarityScore, totalSimilarityScore, bestScoreDiagonal;
	float queryGapExtendScore, queryGapOpenScore;
	float referenceGapExtendScore, referenceGapOpenScore, currentAnchorGapScore;

	uint64_t BestColumn = 0;
	uint64_t BestRow    = 0;
	float BestScore         = FLOAT_NEGATIVE_INFINITY;

	for(uint64_t i = 1, k = queryLen; i < referenceLen; i++, k += queryLen) {

		currentAnchorGapScore = FLOAT_NEGATIVE_INFINITY;
		bestScoreDiagonal = mBestScores(teamid, 0);

		for(uint64_t j = 1, l = k + 1; j < queryLen; j++, l++) {

			// calculate our similarity score
			similarityScore = mScoringMatrix(s1(i - 1) - 'A', s2(j - 1) - 'A');

			// fill the matrices
			totalSimilarityScore = bestScoreDiagonal + similarityScore;

			//cout << "i: " << i << ", j: " << j << ", totalSimilarityScore: " << totalSimilarityScore << endl;

			queryGapExtendScore = mQueryGapScores(teamid, j) - mGapExtendPenalty;
			queryGapOpenScore   = mBestScores(teamid, j) - mGapOpenPenalty;

			// compute the homo-polymer gap score if enabled
			if(mUseHomoPolymerGapOpenPenalty)
				if((j > 1) && (s2(j - 1) == s2(j - 2)))
					queryGapOpenScore = mBestScores(teamid, j) - mHomoPolymerGapOpenPenalty;

			if(queryGapExtendScore > queryGapOpenScore) {
				mQueryGapScores(teamid, j) = queryGapExtendScore;
				mSizesOfVerticalGaps(teamid, l) = (short)(mSizesOfVerticalGaps(teamid, l - queryLen) + 1);
			} else mQueryGapScores(teamid, j) = queryGapOpenScore;

			referenceGapExtendScore = currentAnchorGapScore - mGapExtendPenalty;
			referenceGapOpenScore   = mBestScores(teamid, j - 1) - mGapOpenPenalty;

			// compute the homo-polymer gap score if enabled
			if(mUseHomoPolymerGapOpenPenalty)
				if((i > 1) && (s1(i - 1) == s1(i - 2)))
					referenceGapOpenScore = mBestScores(teamid, j - 1) - mHomoPolymerGapOpenPenalty;

			if(referenceGapExtendScore > referenceGapOpenScore) {
				currentAnchorGapScore = referenceGapExtendScore;
				mSizesOfHorizontalGaps(teamid, l) = (short)(mSizesOfHorizontalGaps(teamid, l - 1) + 1);
			} else currentAnchorGapScore = referenceGapOpenScore;

			bestScoreDiagonal = mBestScores(teamid, j);
			mBestScores(teamid, j) = MaxFloats(totalSimilarityScore, mQueryGapScores(teamid, j), currentAnchorGapScore);
			

			// determine the traceback direction
			// diagonal (445364713) > stop (238960195) > up (214378647) > left (166504495)
			if(mBestScores(teamid, j) == 0)                         mPointers(teamid, l) = DIRECTIONS_STOP;
			else if(mBestScores(teamid, j) == totalSimilarityScore) mPointers(teamid, l) = DIRECTIONS_DIAGONAL;
			else if(mBestScores(teamid, j) == mQueryGapScores(teamid, j))   mPointers(teamid, l) = DIRECTIONS_UP;
			else                                            mPointers(teamid, l) = DIRECTIONS_LEFT;

			// set the traceback start at the current cell i, j and score
			if(mBestScores(teamid, j) > BestScore) {
				BestRow    = i;
				BestColumn = j;
				BestScore  = mBestScores(teamid, j);
			}
		}
	}


	//
	// traceback
	//

	alignment.sw_score = BestScore;
	// aligned sequences
	int gappedAnchorLen  = 0;   // length of sequence #1 after alignment
	int gappedQueryLen   = 0;   // length of sequence #2 after alignment
	int numMismatches    = 0;   // the mismatched nucleotide count

	char c1, c2;

	uint64_t ci = BestRow;
	uint64_t cj = BestColumn;
	uint64_t ck = ci * queryLen;

	// traceback flag
	bool keepProcessing = true;

	while(keepProcessing) {

		// diagonal (445364713) > stop (238960195) > up (214378647) > left (166504495)
		switch(mPointers(teamid, ck + cj)) {

			case DIRECTIONS_DIAGONAL:
				c1 = s1(--ci);
				c2 = s2(--cj);
				ck -= queryLen;

				mReversedAnchor(teamid, gappedAnchorLen++) = c1;
				mReversedQuery(teamid, gappedQueryLen++)   = c2;

				// increment our mismatch counter
				if(mScoringMatrix(c1 - 'A', c2 - 'A') == mMismatchScore) numMismatches++;	
				break;

			case DIRECTIONS_STOP:
				keepProcessing = false;
				break;

			case DIRECTIONS_UP:
				for(uint64_t l = 0, len = mSizesOfVerticalGaps(teamid, ck + cj); l < len; l++) {
					mReversedAnchor(teamid, gappedAnchorLen++) = s1(--ci);
					mReversedQuery(teamid, gappedQueryLen++)   = GAP;
					ck -= queryLen;
					numMismatches++;
				}
				break;

			case DIRECTIONS_LEFT:
				for(uint64_t l = 0, len = mSizesOfHorizontalGaps(teamid, ck + cj); l < len; l++) {
					mReversedAnchor(teamid, gappedAnchorLen++) = GAP;
					mReversedQuery(teamid, gappedQueryLen++)   = s2(--cj);
					numMismatches++;
				}
				break;
		}
	}

	// define the reference and query sequences
	mReversedAnchor(teamid, gappedAnchorLen) = 0;
	mReversedQuery(teamid, gappedQueryLen)   = 0;

	// catch sequences with different lengths
	if(gappedAnchorLen != gappedQueryLen) {
		cout << "ERROR: The aligned sequences have different lengths after Smith-Waterman-Gotoh algorithm." << endl;
		exit(1);
	}

	// reverse the strings and assign them to our alignment structure
  //reverse(mReversedAnchor, mReversedAnchor + gappedAnchorLen);
	//reverse(mReversedQuery,  mReversedQuery  + gappedQueryLen);
  for(uint64_t r=0; r < int(gappedAnchorLen/2); ++r){
    char last = mReversedAnchor(teamid, gappedAnchorLen - r - 1);
    mReversedAnchor(teamid, gappedAnchorLen - r - 1) = mReversedAnchor(teamid, r);
    mReversedAnchor(teamid, r) = last;
  }
  for(uint64_t r=0; r < int(gappedQueryLen/2); ++r){
    char last = mReversedQuery(teamid, gappedQueryLen - r - 1);
    mReversedQuery(teamid, gappedQueryLen - r - 1) = mReversedQuery(teamid, r);
    mReversedQuery(teamid, r) = last;
  }

/*
	cerr << mReversedAnchor << endl;
	cerr << mReversedQuery << endl;
	cerr << endl;
*/

	//alignment.Reference = mReversedAnchor;
	//alignment.Query     = mReversedQuery;

	// set the reference endpoints
	//alignment.ReferenceBegin = ci;
	//alignment.ReferenceEnd   = BestRow - 1;
	alignment.ref_begin   = ci;
	alignment.ref_end     = BestRow - 1;
	alignment.query_begin = cj;
	alignment.query_end   = BestColumn - 1;

	// set the query endpoints
	/*  
	if(alignment.IsReverseComplement) {
		alignment.QueryBegin = s2Length - BestColumn;
		alignment.QueryEnd   = s2Length - cj - 1;
		// alignment.QueryLength= alignment.QueryBegin - alignment.QueryEnd + 1;
	} else {
		alignment.QueryBegin = cj;
		alignment.QueryEnd   = BestColumn - 1;
		// alignment.QueryLength= alignment.QueryEnd - alignment.QueryBegin + 1;
	}
	*/

	// set the query length and number of mismatches
	//alignment.QueryLength = alignment.QueryEnd - alignment.QueryBegin + 1;
	//alignment.NumMismatches  = numMismatches;
  
  string tststring;
  ViewToString(tststring, Kokkos::subview(mReversedAnchor, teamid, Kokkos::ALL));
	uint64_t alLength = strlen(tststring.c_str());
	uint64_t m = 0, d = 0, i = 0;
	bool dashRegion = false;
	//ostringstream oCigar (ostringstream::out);
	for ( unsigned int j = 0; j < alLength; j++ ) {
		// m
		if ( ( mReversedAnchor(teamid, j) != GAP ) && ( mReversedQuery(teamid, j) != GAP ) ) {
			if ( dashRegion ) {
				if ( d != 0 ){
          alignment.cigarCount = d;
          alignment.cigarChar = 'D';
          //oCigar << d << 'D';
        }
				else{
          alignment.cigarCount = i;
          alignment.cigarChar = 'I';
          //oCigar << i << 'I';
        }
			}
			dashRegion = false;
			m++;
			d = 0;
			i = 0;
		}
		else {
			if ( !dashRegion ){
        alignment.cigarChar = 'M';
				//oCigar << m << 'M';
      }
			dashRegion = true;
			m = 0;
			if ( mReversedAnchor(teamid, j) == GAP ) {
				if ( d != 0 ){
          alignment.cigarCount = d;
          alignment.cigarChar = 'D';
          //oCigar << d << 'D';
        }
				i++;
				d = 0;
			}
			else {
				if ( i != 0 ){
          alignment.cigarCount = i;
          alignment.cigarChar = 'I';
          //oCigar << i << 'I';
        }
				d++;
				i = 0;
			}
		}
	}
	if      ( m != 0 ){
    alignment.cigarCount = m;
    alignment.cigarChar = 'M';
    //oCigar << m << 'M';
  }
	else if ( d != 0 ){
    alignment.cigarCount = d;
    alignment.cigarChar = 'D';
    //oCigar << d << 'D';
  }
	else if ( i != 0 ){
    alignment.cigarCount = i;
    alignment.cigarChar = 'I';
    //oCigar << i << 'I';
  }
	
	// fix the gap order
	CorrectHomopolymerGapOrder(teamid, alLength, numMismatches);
}


// corrects the homopolymer gap order for forward alignments
KOKKOS_INLINE_FUNCTION void CSmithWatermanGotoh::CorrectHomopolymerGapOrder(const unsigned int teamid, const unsigned int numBases, const unsigned int numMismatches) const{

	// this is only required for alignments with mismatches
	//if(al.NumMismatches == 0) return;
	if ( numMismatches == 0 ) return;

	// localize the alignment data
	//char* pReference = al.Reference.Data();
	//char* pQuery     = al.Query.Data();
	//const unsigned int numBases = al.Reference.Length();
	View2D<char, Kokkos::LayoutRight>* pReference = const_cast<View2D<char, Kokkos::LayoutRight>*>(&mReversedAnchor);
	View2D<char, Kokkos::LayoutRight>* pQuery     = const_cast<View2D<char, Kokkos::LayoutRight>*>(&mReversedQuery);

	// initialize
	bool hasReferenceGap = false, hasQueryGap = false;
	View2D<char, Kokkos::LayoutRight>* pNonGapSeq = NULL;
	View2D<char, Kokkos::LayoutRight>* pGapSeq    = NULL;
	char nonGapBase  = 'J';

	// identify gapped regions
	for(unsigned int i = 0; i < numBases; i++) {

		// check if the current position is gapped
		hasReferenceGap = false;
		hasQueryGap     = false;

		if((*pReference)(teamid, i) == GAP) {
			hasReferenceGap = true;
			pNonGapSeq      = pQuery;
			pGapSeq         = pReference;
			nonGapBase      = (*pQuery)(teamid, i);
		}

		if((*pQuery)(teamid, i) == GAP) {
			hasQueryGap = true;
			pNonGapSeq  = pReference;
			pGapSeq     = pQuery;
			nonGapBase  = (*pReference)(teamid, i);
		}

		// continue if we don't have any gaps
		if(!hasReferenceGap && !hasQueryGap) continue;

		// sanity check
		if(hasReferenceGap && hasQueryGap) {
			printf("ERROR: Found a gap in both the reference sequence and query sequence.\n");
			exit(1);
		}

		// find the non-gapped length (forward)
		unsigned short numGappedBases = 0;
		unsigned short nonGapLength   = 0;
		unsigned short testPos = i;
		while(testPos < numBases) {

			const char gs  = (*pGapSeq)(teamid, testPos);
			const char ngs = (*pNonGapSeq)(teamid, testPos);

			bool isPartofHomopolymer = false;
			if(((gs == nonGapBase) || (gs == GAP)) && (ngs == nonGapBase)) isPartofHomopolymer = true;
			if(!isPartofHomopolymer) break;

			if(gs == GAP) numGappedBases++;
			else nonGapLength++;
			testPos++;
		}

		// fix the gap order
		if(numGappedBases != 0) {
			
      for(unsigned int kk=0; kk<nonGapLength; kk++) (*pGapSeq)(teamid, i + kk) = nonGapBase;
      for(unsigned int kk=0; kk<numGappedBases; kk++) (*pGapSeq)(teamid, i + nonGapLength + kk) = GAP;
      
      //char* pCurrentSequence = pGapSeq + i;
			//memset(pCurrentSequence, nonGapBase, nonGapLength);
			//pCurrentSequence += nonGapLength;
			//memset(pCurrentSequence, GAP, numGappedBases);
		}

		// increment
		i += numGappedBases + nonGapLength - 1;
	}
}


#endif // SMITHWATERMANGOTOH_H_
