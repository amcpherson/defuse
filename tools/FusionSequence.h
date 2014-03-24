/*
 *  FusionSequence.h
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#ifndef FUSIONSEQUENCE_H_
#define FUSIONSEQUENCE_H_

#include "AlignmentIndex.h"
#include "ExonRegions.h"
#include "FastaIndex.h"
#include "ReadIndex.h"

#include <vector>

using namespace std;
 
class FusionSequence
{
public:
	FusionSequence(const AlignmentIndex& discordant, const AlignmentIndex& anchored, const ReadIndex& reads, const FastaIndex& reference,
				   const ExonRegions& exonRegions, double fragmentLengthMean, double fragmentLengthStdDev, int maxReadLength, int minReadLength);
	
	struct Location
	{
		string refName;
		int strand;
		int start;
		int end;
	};
	
	typedef vector<Location> LocationVec;
	
	bool Calculate(const LocationVec& alignPair, string& sequence, int& breakPos1, int& breakPos2, int& interLen, DoubleVec& stats);
				   
private:
	void FindBreaks(string& sequence, int& breakPos1, int& breakPos2, int& interLen, DoubleVec& stats);
	bool FindBreakDeNovo(string& sequence, int& breakPos1, int& breakPos2, int& interLen, int& minCount) const;
	bool FindBreakSplitReads(string& sequence, int& breakPos1, int& breakPos2, int& splitReadCount, double& splitPosAvg, double& splitMinAvg) const;
	
	void AddAnchoredReads(const string& transcript, int strand, bool revComp, int start, int end);
	void AddSpanningReads(const LocationVec& alignPair);
	
	void CalculateBreakRegion(int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength);
	void CalculateSplitMateRegion(int breakRegionStart, int breakRegionLength, int strand, int& mateRegionStart, int& mateRegionEnd);

	typedef unordered_map<int,string> SequenceMap;
	typedef unordered_map<int,string>::iterator SequenceMapIter;
	typedef unordered_map<int,string>::const_iterator SequenceMapConstIter;
	
	const static int mNumBreakPadding;

	const AlignmentIndex& mDiscordant;
	const AlignmentIndex& mAnchored;
	const ReadIndex& mReads;
	const FastaIndex& mReference;
	const ExonRegions& mExonRegions;

	double mFragmentLengthMean;
	double mFragmentLengthStdDev;
	int mMinFragmentLength;
	int mMaxFragmentLength;
	int mMinReadLength;
	int mMaxReadLength;
		
	string mOutputDirectory;
	
	int mAlignStrand[2];

	Region mSplitAlignSeqPos[2];
	Region mSplitRemainderSeqPos[2];
	Region mDeNovoAlignSeqPos[2];

	string mSplitAlignSeq[2];
	string mSplitRemainderSeq[2];
	string mDeNovoAlignSeq[2];

	SequenceMap mSpanningSequences[2];
	SequenceMap mCandidateSequences;

	string mClusterBreakSequences;
	int mSplitReadCount;
	double mSplitPosAvg;
	double mSplitMinAvg;
};

#endif
