/*
 *  DeNovoFusion.h
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#ifndef DENOVOFUSION_H_
#define DENOVOFUSION_H_

#include "AlignmentIndex.h"
#include "ExonRegions.h"
#include "FastaIndex.h"
#include "ReadIndex.h"

#include <vector>

using namespace std;
 
class DeNovoFusion
{
public:
	DeNovoFusion(const AlignmentIndex& discordant, const AlignmentIndex& anchored, const ReadIndex& reads, const FastaIndex& reference,
	             const ExonRegions& exonRegions, double fragmentLengthMean, double fragmentLengthStdDev, int maxReadLength, int minReadLength,
				 ostream& logFile);
	
	bool Calculate(const LocationVec& alignPair, IntegerVec& breakPos, string& sequence, int& interLen, int& minCount);
	
private:	
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
	
	ostream& mLogFile;
	
	int mDeNovoAlignSeqStart[2];
	int mDeNovoAlignSeqLength[2];

	int mDeNovoSeqStrand[2];
	
	string mDeNovoAlignSeq[2];
	
	SequenceMap mCandidateSequences;
	SequenceMap mSpanningSequences[2];
};

#endif
