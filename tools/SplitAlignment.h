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
#include "SplitReadAligner.h"
#include <boost/unordered_set.hpp>

#include <vector>
#include <map>

using namespace std;
using namespace boost;

struct SplitReadInfo;
struct SplitAlignment;
struct BreakPrediction;

class SplitAlignmentTask
{
public:
	SplitAlignmentTask() : GenerateAlignmentText(false) {}

	bool Initialize(int id, const LocationVec& alignPair, const FastaIndex& reference,
			const ExonRegions& exonRegions, double fragmentLengthMean, 
			double fragmentLengthStdDev, int maxReadLength, int minReadLength);

	void CalculateBreakRegion(int minReadLength, int maxReadLength, int maxFragmentLength,
			int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength);

	vector<SplitAlignment> Align(SplitReadAligner& aligner, const SplitReadInfo& readInfo,
			const string& readSeq);

	SplitAlignment ReAlign(SplitReadAligner& aligner, const ReadIndex& readIndex,
			const SplitAlignment& alignment);
	vector<SplitAlignment> ReAlign(SplitReadAligner& aligner, const ReadIndex& readIndex,
			const vector<SplitAlignment>& alignments);

	BreakPrediction Evaluate(const vector<SplitAlignment>& alignments) const;

	void WriteAlignText(ostream& out, const vector<SplitAlignment>& alignments);

    bool GenerateAlignmentText;

	int mFusionID;

	string mAlignRefName[2];
	int mAlignStrand[2];

	int mSplitAlignSeqStart[2];
	int mSplitAlignSeqLength[2];

	int mSplitSeqStrand[2];
	
	string mSplitAlignSeq[2];
	string mSplitRemainderSeq[2];
	
	LocationVec mMateRegions[2];
};

class BinnedLocations
{
public:
	BinnedLocations(int binSpacing);

	void Add(int id, const Location& location);
	void Overlapping(const RawAlignment& alignment, unordered_set<int>& ids) const;
	
private:	
	int mBinSpacing;
	unordered_map<string,unordered_map<int,IntegerVec> > mBinned[2];
	vector<int> mIDs;
	vector<Region> mRegions;
};

class SplitReadRealigner
{
public:
	SplitReadRealigner();

	void AddTask(const SplitAlignmentTask& alignTasks);
	void AddReads(IReadStream* readStream);
	void DoAlignment(AlignmentStream* mateAlignments, ostream& alignmentsOut);

	SplitReadAligner& GetAligner() { return mAligner; }

private:
	unordered_map<int,SplitAlignmentTask> mAlignTasks;
	unordered_map<int,string> mReads;
	BinnedLocations mBinnedMateRegions;
	SplitReadAligner mAligner;
};

struct SplitReadInfo
{
	ReadID readID;
	bool revComp;
};

struct SplitAlignment
{
	int fusionID;
	SplitReadInfo readInfo;

	pair<int,int> refSplit;
	pair<int,int> readSplit;

	int score;

	string text;

	void WriteAlignment(ostream& out);
	static vector<SplitAlignment> ReadSortedAlignments(istream& in);
};

struct BreakPrediction
{
	int fusionID;

	string alignRefName[2];
	int alignStrand[2];
	int breakPos[2];

	string sequence;

	int splitReadCount;
	double splitPosAvg;
	double splitMinAvg;

	vector<SplitAlignment> alignments;

	void WriteSequence(ostream& out);
	void WriteBreak(ostream& out);
	void WriteAlignments(ostream& out);
};

unordered_map<int,SplitAlignmentTask> CreateTasks(const string& referenceFasta,
		const string& exonRegionsFilename, double fragmentLengthMean, double fragmentLengthStdDev,
		int minReadLength, int maxReadLength, const LocationVecMap& fusionRegions);

#endif

