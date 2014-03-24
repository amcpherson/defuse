/*
 *  SplitAlignment.h
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#ifndef SPLITALIGNMENT_H_
#define SPLITALIGNMENT_H_

#include "AlignmentStream.h"
#include "ExonRegions.h"
#include "FastaIndex.h"
#include "ReadStream.h"
#include "ReadIndex.h"
#include <boost/unordered_set.hpp>

#include <vector>

using namespace std;
using namespace boost;
 
class SplitAlignment
{
public:
	typedef unordered_map<int,SplitAlignment> SplitAlignmentMap;
	typedef unordered_map<int,SplitAlignment>::iterator SplitAlignmentMapIter;
	typedef unordered_map<int,SplitAlignment>::const_iterator SplitAlignmentMapConstIter;

	bool Initialize(const LocationVec& alignPair, const FastaIndex& reference, 
					const ExonRegions& exonRegions, double fragmentLengthMean, 
					double fragmentLengthStdDev, int maxReadLength, int minReadLength);
	
	static void WriteCandidateRefSeqs(ostream& out, SplitAlignmentMap& splitAlignments);
	static void WriteCandidateMateRegions(ostream& out, SplitAlignmentMap& splitAlignments);
	
	static void ReadCandidateRefSeqs(istream& in, SplitAlignmentMap& splitAlignments);
	static void ReadCandidateMateRegions(istream& in, SplitAlignmentMap& splitAlignments);
	
	static bool FindCandidates(AlignmentStream* alignments, SplitAlignmentMap& splitAlignments);
	
	static void ReadCandidateSequences(IReadStream* readStream, SplitAlignmentMap& splitAlignments);
	static void ReadCandidateSequences(const ReadIndex& readIndex, SplitAlignmentMap& splitAlignments);

	bool Align(bool generateAlignmentText);

	static void WriteAlignments(ostream& out, SplitAlignmentMap& splitAlignments);
	static void ReadAlignments(istream& in, SplitAlignmentMap& splitAlignments);
	
	bool Evaluate();
	
	static void WriteSequences(ostream& out, SplitAlignmentMap& splitAlignments);
	static void WriteBreaks(ostream& out, SplitAlignmentMap& splitAlignments);
	static void WriteAlignText(ostream& out, SplitAlignmentMap& splitAlignments);
	
private:	
	inline void CalculateBreakRegion(int minReadLength, int maxReadLength, int maxFragmentLength, int alignStart, 
									 int alignEnd, int strand, int& breakStart, int& breakLength);
	inline void CalculateSplitMateRegion(int minReadLength, int maxReadLength, int minFragmentLength, 
										 int maxFragmentLength, int breakRegionStart, int breakRegionLength, 
										 int strand, int& mateRegionStart, int& mateRegionEnd);

	string mAlignRefName[2];
	int mAlignStrand[2];

	int mSplitAlignSeqStart[2];
	int mSplitAlignSeqLength[2];

	int mSplitSeqStrand[2];
	
	string mSplitAlignSeq[2];
	string mSplitRemainderSeq[2];
	
	LocationVec mMateRegions[2];
	
	IntegerVec mCandidateReadID;
	IntegerVec mCandidateRevComp;
	
	StringVec mCandidateSequence;
	
	IntegerVec mAlignmentReadID;
	IntegerPairVec mAlignmentRefSplit;
	IntegerPairVec mAlignmentReadSplit;
	IntegerVec mAlignmentScore;
	StringVec mAlignmentText;
	
	IntegerPair mBestSplit;
	IntegerVec mBreakPos;
	string mSequence;
	int mSplitReadCount;
	double mSplitPosAvg;
	double mSplitMinAvg;	
};

#endif

