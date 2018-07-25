#include "SmithWatermanGotoh.h"

#include <stdint.h>

const float CSmithWatermanGotoh::FLOAT_NEGATIVE_INFINITY = (float)-1e+30;

const char CSmithWatermanGotoh::Directions_STOP     = 0;
const char CSmithWatermanGotoh::Directions_LEFT     = 1;
const char CSmithWatermanGotoh::Directions_DIAGONAL = 2;
const char CSmithWatermanGotoh::Directions_UP       = 3;

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
  mPointers              = Kokkos::View<char*>("mPointers", mCurrentMatrixSize);
  mSizesOfVerticalGaps   = Kokkos::View<short*>("mSizesOfVerticalGaps", mCurrentMatrixSize);
  mSizesOfHorizontalGaps = Kokkos::View<short*>("mSizesOfHorizontalGaps", mCurrentMatrixSize);
  mQueryGapScores        = Kokkos::View<float*>("mQueryGapScores", mCurrentQuerySize + 1);
  mBestScores            = Kokkos::View<float*>("mBestScores", mCurrentQuerySize + 1);
  mReversedAnchor        = Kokkos::View<char*>("mReversedAnchor", mCurrentAQSumSize + 1);	// reversed sequence #1
  mReversedQuery         = Kokkos::View<char*>("mReversedQuery", mCurrentAQSumSize + 1);	// reversed sequence #2
  
}


void CSmithWatermanGotoh::InitArrays(unsigned int max_reference_length, unsigned int max_sequence_length){
  mCurrentMatrixSize = (max_reference_length + 1) * (max_sequence_length + 1);
  mCurrentQuerySize = max_sequence_length + 1;
  mCurrentAQSumSize = max_sequence_length + max_reference_length;
  InitArrays();
}


CSmithWatermanGotoh::~CSmithWatermanGotoh(void) {}

// string to view and view to string converter
void ViewToString(string& result, Kokkos::View<char*> view){
  result.clear();
  for(uint64_t l=0; l<view.extent(0); l++){
      result+=view(l);
  }
}

//convert string to view:
void StringToView(Kokkos::View<char*> view, const string input){
  for(uint64_t l=0; l<input.size(); l++){
    view(l) = input[l];
  }
  if(input.size() < view.extent(0)) view(input.size())=0;
  //for(uint64_t l=input.size(); l<view.extent(0); l++){
  //  view(l) = 0;
  //}
}


// aligns the query sequence to the reference using the Smith Waterman Gotoh algorithm
void CSmithWatermanGotoh::Align(Alignment& alignment, const char* s1, const unsigned int s1Length, const char* s2, const unsigned int& s2Length) {

	if((s1Length == 0) || (s2Length == 0)) {
		cout << "ERROR: Found a read with a zero length." << endl;
		exit(1);
	}

	uint64_t referenceLen      = s1Length + 1;
	uint64_t queryLen          = s2Length + 1;
	uint64_t sequenceSumLength = s1Length + s2Length;
	// set current matrix size
	uint64_t matrix_size = referenceLen * queryLen;
	mCurrentMatrixSize = matrix_size;

	// initialize the traceback matrix to STOP
  for(uint64_t i = 0; i < queryLen; ++i) mPointers(i) = 0;
	//memset((char*)mPointers, 0, SIZEOF_CHAR * queryLen);
	for(uint64_t i = 1; i < referenceLen; ++i) {
	  mPointers(i * queryLen) = 0;
	}

	// initialize the gap matrices to 1
	//uninitialized_fill(mSizesOfVerticalGaps, mSizesOfVerticalGaps + mCurrentMatrixSize, 1);
	//uninitialized_fill(mSizesOfHorizontalGaps, mSizesOfHorizontalGaps + mCurrentMatrixSize, 1);
  for(uint64_t i=0; i < mCurrentMatrixSize; ++i){
    mSizesOfVerticalGaps(i) = 1;
    mSizesOfHorizontalGaps(i) = 1;
  }

	//
	// sequence lengths
	//
	mCurrentQuerySize = s2Length;
	mCurrentAQSumSize = sequenceSumLength;

	// initialize the gap score and score vectors
  //uninitialized_fill(mQueryGapScores, mQueryGapScores + queryLen, FLOAT_NEGATIVE_INFINITY);
  //memset((char*)mBestScores, 0, SIZEOF_FLOAT * queryLen);
  for(uint64_t i=0; i < queryLen; ++i){
    mQueryGapScores(i) = FLOAT_NEGATIVE_INFINITY;
    mBestScores(i) = 0.;
  }

	float similarityScore, totalSimilarityScore, bestScoreDiagonal;
	float queryGapExtendScore, queryGapOpenScore;
	float referenceGapExtendScore, referenceGapOpenScore, currentAnchorGapScore;

	uint64_t BestColumn = 0;
	uint64_t BestRow    = 0;
	float BestScore         = FLOAT_NEGATIVE_INFINITY;

	for(uint64_t i = 1, k = queryLen; i < referenceLen; i++, k += queryLen) {

		currentAnchorGapScore = FLOAT_NEGATIVE_INFINITY;
		bestScoreDiagonal = mBestScores(0);

		for(uint64_t j = 1, l = k + 1; j < queryLen; j++, l++) {

			// calculate our similarity score
			similarityScore = mScoringMatrix(s1[i - 1] - 'A', s2[j - 1] - 'A');

			// fill the matrices
			totalSimilarityScore = bestScoreDiagonal + similarityScore;

			//cout << "i: " << i << ", j: " << j << ", totalSimilarityScore: " << totalSimilarityScore << endl;

			queryGapExtendScore = mQueryGapScores(j) - mGapExtendPenalty;
			queryGapOpenScore   = mBestScores(j) - mGapOpenPenalty;

			// compute the homo-polymer gap score if enabled
			if(mUseHomoPolymerGapOpenPenalty)
				if((j > 1) && (s2[j - 1] == s2[j - 2]))
					queryGapOpenScore = mBestScores(j) - mHomoPolymerGapOpenPenalty;

			if(queryGapExtendScore > queryGapOpenScore) {
				mQueryGapScores(j) = queryGapExtendScore;
				mSizesOfVerticalGaps(l) = (short)(mSizesOfVerticalGaps(l - queryLen) + 1);
			} else mQueryGapScores(j) = queryGapOpenScore;

			referenceGapExtendScore = currentAnchorGapScore - mGapExtendPenalty;
			referenceGapOpenScore   = mBestScores(j - 1) - mGapOpenPenalty;

			// compute the homo-polymer gap score if enabled
			if(mUseHomoPolymerGapOpenPenalty)
				if((i > 1) && (s1[i - 1] == s1[i - 2]))
					referenceGapOpenScore = mBestScores(j - 1) - mHomoPolymerGapOpenPenalty;

			if(referenceGapExtendScore > referenceGapOpenScore) {
				currentAnchorGapScore = referenceGapExtendScore;
				mSizesOfHorizontalGaps(l) = (short)(mSizesOfHorizontalGaps(l - 1) + 1);
			} else currentAnchorGapScore = referenceGapOpenScore;

			bestScoreDiagonal = mBestScores(j);
			mBestScores(j) = MaxFloats(totalSimilarityScore, mQueryGapScores(j), currentAnchorGapScore);
			

			// determine the traceback direction
			// diagonal (445364713) > stop (238960195) > up (214378647) > left (166504495)
			if(mBestScores(j) == 0)                         mPointers(l) = Directions_STOP;
			else if(mBestScores(j) == totalSimilarityScore) mPointers(l) = Directions_DIAGONAL;
			else if(mBestScores(j) == mQueryGapScores(j))   mPointers(l) = Directions_UP;
			else                                            mPointers(l) = Directions_LEFT;

			// set the traceback start at the current cell i, j and score
			if(mBestScores(j) > BestScore) {
				BestRow    = i;
				BestColumn = j;
				BestScore  = mBestScores(j);
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
		switch(mPointers[ck + cj]) {

			case Directions_DIAGONAL:
				c1 = s1[--ci];
				c2 = s2[--cj];
				ck -= queryLen;

				mReversedAnchor(gappedAnchorLen++) = c1;
				mReversedQuery(gappedQueryLen++)   = c2;

				// increment our mismatch counter
				if(mScoringMatrix(c1 - 'A', c2 - 'A') == mMismatchScore) numMismatches++;	
				break;

			case Directions_STOP:
				keepProcessing = false;
				break;

			case Directions_UP:
				for(uint64_t l = 0, len = mSizesOfVerticalGaps(ck + cj); l < len; l++) {
					mReversedAnchor(gappedAnchorLen++) = s1[--ci];
					mReversedQuery(gappedQueryLen++)   = GAP;
					ck -= queryLen;
					numMismatches++;
				}
				break;

			case Directions_LEFT:
				for(uint64_t l = 0, len = mSizesOfHorizontalGaps(ck + cj); l < len; l++) {
					mReversedAnchor(gappedAnchorLen++) = GAP;
					mReversedQuery(gappedQueryLen++)   = s2[--cj];
					numMismatches++;
				}
				break;
		}
	}

	// define the reference and query sequences
	mReversedAnchor(gappedAnchorLen) = 0;
	mReversedQuery(gappedQueryLen)   = 0;

	// catch sequences with different lengths
	if(gappedAnchorLen != gappedQueryLen) {
		cout << "ERROR: The aligned sequences have different lengths after Smith-Waterman-Gotoh algorithm." << endl;
		exit(1);
	}

	// reverse the strings and assign them to our alignment structure
  //reverse(mReversedAnchor, mReversedAnchor + gappedAnchorLen);
	//reverse(mReversedQuery,  mReversedQuery  + gappedQueryLen);
  for(uint64_t r=0; r < int(gappedAnchorLen/2); ++r){
    char last = mReversedAnchor(gappedAnchorLen - r - 1);
    mReversedAnchor(gappedAnchorLen - r - 1) = mReversedAnchor(r);
    mReversedAnchor(r) = last;
  }
  for(uint64_t r=0; r < int(gappedQueryLen/2); ++r){
    char last = mReversedQuery(gappedQueryLen - r - 1);
    mReversedQuery(gappedQueryLen - r - 1) = mReversedQuery(r);
    mReversedQuery(r) = last;
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
  ViewToString(tststring, mReversedAnchor);
	uint64_t alLength = strlen(tststring.c_str());
	uint64_t m = 0, d = 0, i = 0;
	bool dashRegion = false;
	//ostringstream oCigar (ostringstream::out);
	for ( unsigned int j = 0; j < alLength; j++ ) {
		// m
		if ( ( mReversedAnchor(j) != GAP ) && ( mReversedQuery(j) != GAP ) ) {
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
			if ( mReversedAnchor(j) == GAP ) {
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
	CorrectHomopolymerGapOrder(alLength, numMismatches);
}

// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
void CSmithWatermanGotoh::CreateScoringMatrix(void) {

	unsigned int nIndex = 13;
	unsigned int xIndex = 23;
  
  //allocate memory
  mScoringMatrix = Kokkos::View<float**>("mScoringMatrix", MOSAIK_NUM_NUCLEOTIDES, MOSAIK_NUM_NUCLEOTIDES);
  
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

// corrects the homopolymer gap order for forward alignments
void CSmithWatermanGotoh::CorrectHomopolymerGapOrder(const unsigned int numBases, const unsigned int numMismatches) {

	// this is only required for alignments with mismatches
	//if(al.NumMismatches == 0) return;
	if ( numMismatches == 0 ) return;

	// localize the alignment data
	//char* pReference = al.Reference.Data();
	//char* pQuery     = al.Query.Data();
	//const unsigned int numBases = al.Reference.Length();
	Kokkos::View<char*>* pReference = &mReversedAnchor;
	Kokkos::View<char*>* pQuery     = &mReversedQuery;

	// initialize
	bool hasReferenceGap = false, hasQueryGap = false;
	Kokkos::View<char*>* pNonGapSeq = NULL;
	Kokkos::View<char*>* pGapSeq    = NULL;
	char nonGapBase  = 'J';

	// identify gapped regions
	for(unsigned int i = 0; i < numBases; i++) {

		// check if the current position is gapped
		hasReferenceGap = false;
		hasQueryGap     = false;

		if((*pReference)(i) == GAP) {
			hasReferenceGap = true;
			pNonGapSeq      = pQuery;
			pGapSeq         = pReference;
			nonGapBase      = (*pQuery)(i);
		}

		if((*pQuery)(i) == GAP) {
			hasQueryGap = true;
			pNonGapSeq  = pReference;
			pGapSeq     = pQuery;
			nonGapBase  = (*pReference)(i);
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

			const char gs  = (*pGapSeq)(testPos);
			const char ngs = (*pNonGapSeq)(testPos);

			bool isPartofHomopolymer = false;
			if(((gs == nonGapBase) || (gs == GAP)) && (ngs == nonGapBase)) isPartofHomopolymer = true;
			if(!isPartofHomopolymer) break;

			if(gs == GAP) numGappedBases++;
			else nonGapLength++;
			testPos++;
		}

		// fix the gap order
		if(numGappedBases != 0) {
			
      for(unsigned int kk=0; kk<nonGapLength; kk++) (*pGapSeq)(i + kk) = nonGapBase;
      for(unsigned int kk=0; kk<numGappedBases; kk++) (*pGapSeq)(i + nonGapLength + kk) = GAP;
      
      //char* pCurrentSequence = pGapSeq + i;
			//memset(pCurrentSequence, nonGapBase, nonGapLength);
			//pCurrentSequence += nonGapLength;
			//memset(pCurrentSequence, GAP, numGappedBases);
		}

		// increment
		i += numGappedBases + nonGapLength - 1;
	}
}
