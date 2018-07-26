#include "SmithWatermanGotoh.h"

#include <stdint.h>

//const float CSmithWatermanGotoh::FLOAT_NEGATIVE_INFINITY = (float)-1e+30;

//const char CSmithWatermanGotoh::Directions_STOP     = 0;
//const char CSmithWatermanGotoh::Directions_LEFT     = 1;
//const char CSmithWatermanGotoh::Directions_DIAGONAL = 2;
//const char CSmithWatermanGotoh::Directions_UP       = 3;

CSmithWatermanGotoh::CSmithWatermanGotoh(float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty,
					 unsigned int max_reference_length, unsigned int max_sequence_length) 
: mCurrentMatrixSize((max_reference_length + 1) * (max_sequence_length + 1))
, mCurrentAnchorSize(0)
, mCurrentQuerySize(max_sequence_length + 1)
, mCurrentAQSumSize(max_sequence_length + max_reference_length)
, mMatchScore(matchScore)
, mMismatchScore(mismatchScore)
, mGapOpenPenalty(gapOpenPenalty)
, mGapExtendPenalty(gapExtendPenalty)
, mPointers(NULL)
, mSizesOfVerticalGaps(NULL)
, mSizesOfHorizontalGaps(NULL)
, mQueryGapScores(NULL)
, mBestScores(NULL)
, mReversedAnchor(NULL)
, mReversedQuery(NULL)
, mUseHomoPolymerGapOpenPenalty(false)
{
	CreateScoringMatrix();
  InitArrays();
}


void CSmithWatermanGotoh::InitArrays(){
  
  //create views
  mPointers              = View1D<char>("mPointers", mCurrentMatrixSize);
  mSizesOfVerticalGaps   = View1D<short>("mSizesOfVerticalGaps", mCurrentMatrixSize);
  mSizesOfHorizontalGaps = View1D<short>("mSizesOfHorizontalGaps", mCurrentMatrixSize);
  mQueryGapScores        = View1D<float>("mQueryGapScores", mCurrentQuerySize + 1);
  mBestScores            = View1D<float>("mBestScores", mCurrentQuerySize + 1);
  mReversedAnchor        = View1D<char>("mReversedAnchor", mCurrentAQSumSize + 1);	// reversed sequence #1
  mReversedQuery         = View1D<char>("mReversedQuery", mCurrentAQSumSize + 1);	// reversed sequence #2
  
}


void CSmithWatermanGotoh::InitArrays(unsigned int max_reference_length, unsigned int max_sequence_length){
  mCurrentMatrixSize = (max_reference_length + 1) * (max_sequence_length + 1);
  mCurrentQuerySize = max_sequence_length + 1;
  mCurrentAQSumSize = max_sequence_length + max_reference_length;
  InitArrays();
}


CSmithWatermanGotoh::~CSmithWatermanGotoh(void) {}


// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
void CSmithWatermanGotoh::CreateScoringMatrix(void) {

	unsigned int nIndex = 13;
	unsigned int xIndex = 23;
  
  //allocate memory
  mScoringMatrix = View2D<float>("mScoringMatrix", MOSAIK_NUM_NUCLEOTIDES, MOSAIK_NUM_NUCLEOTIDES);
  
	// define the N score to be 1/4 of the span between mismatch and match
	//const short nScore = mMismatchScore + (short)(((mMatchScore - mMismatchScore) / 4.0) + 0.5);

	// calculate the scoring matrix
	for(unsigned char i = 0; i < MOSAIK_NUM_NUCLEOTIDES; i++) {
    
		for(unsigned char j = 0; j < MOSAIK_NUM_NUCLEOTIDES; j++) {

			// N.B. matching N to everything (while conceptually correct) leads to some
			// bad alignments, lets make N be a mismatch instead.

			// add the matches or mismatches to the hashtable (N is a mismatch)
			if((i == nIndex) || (j == nIndex)) mScoringMatrix(i,j) = mMismatchScore;
			else if((i == xIndex) || (j == xIndex)) mScoringMatrix(i,j) = mMismatchScore;
			else if(i == j) mScoringMatrix(i,j) = mMatchScore;
			else mScoringMatrix(i,j) = mMismatchScore;
		}
	}

	// add ambiguity codes
	mScoringMatrix('M' - 'A', 'A' - 'A') = mMatchScore;	// M - A
	mScoringMatrix('A' - 'A', 'M' - 'A') = mMatchScore;
	mScoringMatrix('M' - 'A', 'C' - 'A') = mMatchScore; // M - C
	mScoringMatrix('C' - 'A', 'M' - 'A') = mMatchScore;

	mScoringMatrix('R' - 'A', 'A' - 'A') = mMatchScore;	// R - A
	mScoringMatrix('A' - 'A', 'R' - 'A') = mMatchScore;
	mScoringMatrix('R' - 'A', 'G' - 'A') = mMatchScore; // R - G
	mScoringMatrix('G' - 'A', 'R' - 'A') = mMatchScore;

	mScoringMatrix('W' - 'A', 'A' - 'A') = mMatchScore;	// W - A
	mScoringMatrix('A' - 'A', 'W' - 'A') = mMatchScore;
	mScoringMatrix('W' - 'A', 'T' - 'A') = mMatchScore; // W - T
	mScoringMatrix('T' - 'A', 'W' - 'A') = mMatchScore;

	mScoringMatrix('S' - 'A', 'C' - 'A') = mMatchScore;	// S - C
	mScoringMatrix('C' - 'A', 'S' - 'A') = mMatchScore;
	mScoringMatrix('S' - 'A', 'G' - 'A') = mMatchScore; // S - G
	mScoringMatrix('G' - 'A', 'S' - 'A') = mMatchScore;

	mScoringMatrix('Y' - 'A', 'C' - 'A') = mMatchScore;	// Y - C
	mScoringMatrix('C' - 'A', 'Y' - 'A') = mMatchScore;
	mScoringMatrix('Y' - 'A', 'T' - 'A') = mMatchScore; // Y - T
	mScoringMatrix('T' - 'A', 'Y' - 'A') = mMatchScore;

	mScoringMatrix('K' - 'A', 'G' - 'A') = mMatchScore;	// K - G
	mScoringMatrix('G' - 'A', 'K' - 'A') = mMatchScore;
	mScoringMatrix('K' - 'A', 'T' - 'A') = mMatchScore; // K - T
	mScoringMatrix('T' - 'A', 'K' - 'A') = mMatchScore;

	mScoringMatrix('V' - 'A', 'A' - 'A') = mMatchScore;	// V - A
	mScoringMatrix('A' - 'A', 'V' - 'A') = mMatchScore;
	mScoringMatrix('V' - 'A', 'C' - 'A') = mMatchScore; // V - C
	mScoringMatrix('C' - 'A', 'V' - 'A') = mMatchScore;
	mScoringMatrix('V' - 'A', 'G' - 'A') = mMatchScore; // V - G
	mScoringMatrix('G' - 'A', 'V' - 'A') = mMatchScore;

	mScoringMatrix('H' - 'A', 'A' - 'A') = mMatchScore;	// H - A
	mScoringMatrix('A' - 'A', 'H' - 'A') = mMatchScore;
	mScoringMatrix('H' - 'A', 'C' - 'A') = mMatchScore; // H - C
	mScoringMatrix('C' - 'A', 'H' - 'A') = mMatchScore;
	mScoringMatrix('H' - 'A', 'T' - 'A') = mMatchScore; // H - T
	mScoringMatrix('T' - 'A', 'H' - 'A') = mMatchScore;

	mScoringMatrix('D' - 'A', 'A' - 'A') = mMatchScore;	// D - A
	mScoringMatrix('A' - 'A', 'D' - 'A') = mMatchScore;
	mScoringMatrix('D' - 'A', 'G' - 'A') = mMatchScore; // D - G
	mScoringMatrix('G' - 'A', 'D' - 'A') = mMatchScore;
	mScoringMatrix('D' - 'A', 'T' - 'A') = mMatchScore; // D - T
	mScoringMatrix('T' - 'A', 'D' - 'A') = mMatchScore;

	mScoringMatrix('B' - 'A', 'C' - 'A') = mMatchScore;	// B - C
	mScoringMatrix('C' - 'A', 'B' - 'A') = mMatchScore;
	mScoringMatrix('B' - 'A', 'G' - 'A') = mMatchScore; // B - G
	mScoringMatrix('G' - 'A', 'B' - 'A') = mMatchScore;
	mScoringMatrix('B' - 'A', 'T' - 'A') = mMatchScore; // B - T
	mScoringMatrix('T' - 'A', 'B' - 'A') = mMatchScore;
}

// enables homo-polymer scoring
void CSmithWatermanGotoh::EnableHomoPolymerGapPenalty(float hpGapOpenPenalty) {
	mUseHomoPolymerGapOpenPenalty = true;
	mHomoPolymerGapOpenPenalty    = hpGapOpenPenalty;
}

