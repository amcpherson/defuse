/*
 *  FusionSequence.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "FusionSequence.h"
#include "GenomeIndex.h"
#include "PathFinder.h"
#include "SplitReadAligner.h"
#include "DebugCheck.h"

#include <map>
#include <list>
#include <fstream>
#include <sstream>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

const int FusionSequence::mNumBreakPadding = 10;

FusionSequence::FusionSequence(const AlignmentIndex& discordant, const AlignmentIndex& anchored, const ReadIndex& reads, const FastaIndex& reference,
							   const ExonRegions& exonRegions, double fragmentLengthMean, double fragmentLengthStdDev, int minReadLength, int maxReadLength)
: mDiscordant(discordant), mAnchored(anchored), mReads(reads), mReference(reference), mFragmentLengthMean(fragmentLengthMean),
  mExonRegions(exonRegions), mFragmentLengthStdDev(fragmentLengthStdDev), mMinReadLength(minReadLength), mMaxReadLength(maxReadLength)
{
	mMinFragmentLength = mFragmentLengthMean - 3 * mFragmentLengthStdDev;
	mMaxFragmentLength = mFragmentLengthMean + 3 * mFragmentLengthStdDev;
}

bool FusionSequence::Calculate(const LocationVec& alignPair, string& sequence, int& breakPos1, int& breakPos2, int& interLen, DoubleVec& stats)
{
	if (alignPair.size() != 2)
	{
		cerr << "Error: Incorrect input for FusionSequence::Calculate()" << endl;
		return false;
	}
	
	AddSpanningReads(alignPair);
		
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		const string& alignRefName = alignPair[clusterEnd].refName;
		int alignStrand = alignPair[clusterEnd].strand;
		int alignStart = alignPair[clusterEnd].start;
		int alignEnd = alignPair[clusterEnd].end;
		
		int refSeqStrand = (clusterEnd == 0) ? alignStrand : OtherStrand(alignStrand);
		bool revCompReads = (clusterEnd == 0) ? true : false;
		
		int breakRegionStart;
		int breakRegionLength;
		CalculateBreakRegion(alignStart, alignEnd, alignStrand, breakRegionStart, breakRegionLength);
		
		mAlignStrand[clusterEnd] = alignStrand;
		
		// Assumption: the break region we get from the paired end analysis will not involve a read that pushes more than 50%
		// into the breakpoint.  This is aligner dependent.
		if (alignStrand == PlusStrand)
		{
			mSplitAlignSeqPos[clusterEnd].start = breakRegionStart - (int)(1.0 * mMaxReadLength);
			mSplitAlignSeqPos[clusterEnd].end = breakRegionStart + breakRegionLength - 1;

			if (alignStart < mSplitRemainderSeqPos[clusterEnd].start)
			{
				mSplitRemainderSeqPos[clusterEnd].start = alignStart;
				mSplitRemainderSeqPos[clusterEnd].end = mSplitRemainderSeqPos[clusterEnd].start - 1;
			}
			else
			{
				mSplitRemainderSeqPos[clusterEnd].start = 0;
				mSplitRemainderSeqPos[clusterEnd].end = -1;
			}

			mDeNovoAlignSeqPos[clusterEnd].start = alignStart;
			mDeNovoAlignSeqPos[clusterEnd].end = breakRegionStart + breakRegionLength - 1;
		}
		else
		{
			mSplitAlignSeqPos[clusterEnd].start = breakRegionStart - breakRegionLength + 1;
			mSplitAlignSeqPos[clusterEnd].end = breakRegionStart + (int)(1.0 * mMaxReadLength);

			if (alignEnd > mSplitAlignSeqPos[clusterEnd].end)
			{
				mSplitRemainderSeqPos[clusterEnd].start = mSplitAlignSeqPos[clusterEnd].end + 1;
				mSplitRemainderSeqPos[clusterEnd].end = alignEnd;
			}
			else
			{
				mSplitRemainderSeqPos[clusterEnd].start = 0;
				mSplitRemainderSeqPos[clusterEnd].end = -1;
			}

			mDeNovoAlignSeqPos[clusterEnd].start = breakRegionStart - breakRegionLength + 1;
			mDeNovoAlignSeqPos[clusterEnd].end = alignEnd;
		}
		
		// Retrieve split reference alignment sequence
		mReference.Get(alignRefName, refSeqStrand, mSplitAlignSeqPos[clusterEnd].start, mSplitAlignSeqPos[clusterEnd].end, mSplitAlignSeq[clusterEnd]);
		
		// Retrieve split reference alignment remainder sequence
		mReference.Get(alignRefName, refSeqStrand, mSplitRemainderSeqPos[clusterEnd].start, mSplitRemainderSeqPos[clusterEnd].end, mSplitRemainderSeq[clusterEnd]);
		
		// Retrieve denovo reference alignment sequence
		mReference.Get(alignRefName, refSeqStrand, mDeNovoAlignSeqPos[clusterEnd].start, mDeNovoAlignSeqPos[clusterEnd].end, mDeNovoAlignSeq[clusterEnd]);
				
		// Find gene for this transcript
		const string& gene = mExonRegions.GetTranscriptGene(alignRefName);
		
		// Find all transcripts of the same gene
		const StringVec& geneTranscripts = mExonRegions.GetGeneTranscripts(gene);
		
		// Iterate through all transcripts of this gene
		for (StringVecConstIter geneTranscriptIter = geneTranscripts.begin(); geneTranscriptIter != geneTranscripts.end(); geneTranscriptIter++)
		{
			const string& remapTranscript = *geneTranscriptIter;
			
			// Remap break region start to this transcript if possible
			int remappedBreakRegionStart;
			if (!mExonRegions.RemapTranscriptToTranscript(alignRefName, breakRegionStart, remapTranscript, remappedBreakRegionStart))
			{
				continue;
			}
			
			DebugCheck(remapTranscript != alignRefName || breakRegionStart == remappedBreakRegionStart);					
			
			// Calculate the anchored mate alignment region in this transcript
			int mateRegionStart;
			int mateRegionEnd;
			CalculateSplitMateRegion(remappedBreakRegionStart, breakRegionLength, alignStrand, mateRegionStart, mateRegionEnd);

			AddAnchoredReads(remapTranscript, alignStrand, revCompReads, mateRegionStart, mateRegionEnd);
		}
	}
	
	FindBreaks(sequence, breakPos1, breakPos2, interLen, stats);
}

void FusionSequence::AddAnchoredReads(const string& transcript, int strand, bool revComp, int start, int end)
{
	CompAlignVec candidateMates;
	
	NameIndex referenceNames;
	mDiscordant.Find(transcript, strand, start, end, referenceNames, candidateMates);
	mAnchored.Find(transcript, strand, start, end, referenceNames, candidateMates);
	
	for (int candidateIndex = 0; candidateIndex < candidateMates.size(); candidateIndex++)
	{
		string readSeqs[2];
		mReads.Find(candidateMates[candidateIndex].readID.fragmentIndex, readSeqs[0], readSeqs[1]);
		
		ReadID candidateReadID;
		candidateReadID.fragmentIndex = candidateMates[candidateIndex].readID.fragmentIndex;
		candidateReadID.readEnd = (candidateMates[candidateIndex].readID.readEnd == 0) ? 1 : 0;
		
		string candidateSeq = readSeqs[candidateReadID.readEnd];
		
		if (revComp)
		{
			ReverseComplement(candidateSeq);
		}
		
		mCandidateSequences[candidateReadID.id] = candidateSeq;
	}
}

void FusionSequence::AddSpanningReads(const LocationVec& alignPair)
{
	typedef unordered_set<int> IntegerSet;
	typedef unordered_map<int,IntegerSet> IntSetMap;
	typedef unordered_map<int,CompactAlignment> CompAlignMap;
	
	IntegerSet fragments;
	IntSetMap spanningReadEnd[2];
	CompAlignMap spanningAlignments[2];
	
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		const string& alignRefName = alignPair[clusterEnd].refName;
		int alignStrand = alignPair[clusterEnd].strand;
		int alignStart = alignPair[clusterEnd].start;
		int alignEnd = alignPair[clusterEnd].end;

		NameIndex referenceNames;
		CompAlignVec alignments;
		
		mDiscordant.Find(alignRefName, alignStrand, alignStart, alignEnd, referenceNames, alignments);
	
		for (CompAlignVecConstIter alignIter = alignments.begin(); alignIter != alignments.end(); alignIter++)
		{
			fragments.insert(alignIter->readID.fragmentIndex);
			spanningReadEnd[clusterEnd][alignIter->readID.fragmentIndex].insert(alignIter->readID.readEnd);
			spanningAlignments[clusterEnd].insert(make_pair(alignIter->readID.fragmentIndex,*alignIter));			
		}
	}
	
	for (unordered_set<int>::const_iterator fragmentIter = fragments.begin(); fragmentIter != fragments.end(); fragmentIter++)
	{
		int fragmentIndex = *fragmentIter;
		
		if (spanningReadEnd[0][fragmentIndex].size() != 1 || spanningReadEnd[1][fragmentIndex].size() != 1)
		{
			continue;
		}
		
		string readSeqs[2];
		mReads.Find(fragmentIndex, readSeqs[0], readSeqs[1]);
		
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			ReadID spanningReadID;
			spanningReadID.fragmentIndex = fragmentIndex;
			spanningReadID.readEnd = *(spanningReadEnd[clusterEnd][fragmentIndex].begin());
			
			mSpanningSequences[clusterEnd][spanningReadID.id] = readSeqs[spanningReadID.readEnd];			
		}
	}
}

void FusionSequence::FindBreaks(string& sequence, int& breakPos1, int& breakPos2, int& interLen, DoubleVec& stats)
{
	// Generate output
	string splitAlignOutputFilename = "split_alignments.txt";
	ofstream splitAlignOutputFile(splitAlignOutputFilename.c_str());

	if (!splitAlignOutputFile.good())
	{
		cerr << "Error: Unable to open " << splitAlignOutputFilename << endl;
		exit(1);
	}

	string deNovoBreakSequence;
	int deNovoBreakPos1;
	int deNovoBreakPos2;
	int deNovoInterLen;
	int deNovoMinCount;

	bool deNovoResult = FindBreakDeNovo(deNovoBreakSequence, deNovoBreakPos1, deNovoBreakPos2, deNovoInterLen, deNovoMinCount);

	string splitBreakSequence;
	int splitBreakPos1;
	int splitBreakPos2;
	int splitReadCount;
	double splitPosAvg;
	double splitMinAvg;
	double splitFragmentMean;
	
	bool splitResult = false;
	if (mCandidateSequences.size() < 20000)
	{
		splitResult = FindBreakSplitReads(splitBreakSequence, splitBreakPos1, splitBreakPos2, splitReadCount, splitPosAvg, splitMinAvg);
	}
	
	// Select the prediction with the sample fragment length mean closest to the overall fragment length mean
	double deNovoFragmentMeanDiff = abs(mFragmentLengthMean - deNovoFragmentMean);
	double splitFragmentMeanDiff = abs(mFragmentLengthMean - splitFragmentMean);
	
	if (deNovoResult && (!splitResult || deNovoFragmentMeanDiff < splitFragmentMeanDiff))
	{
		mClusterBreakSequences = deNovoBreakSequence;
		mSplitReadCount = deNovoMinCount;
		mSplitPosAvg = 1.0;
		mSplitMinAvg = 1.0;
	}
	else if (splitResult)
	{
		mClusterBreakSequences = splitBreakSequence;
		mSplitReadCount = splitReadCount;
	}
	else
	{
		mClusterBreakSequences = "";
		mSplitReadCount = 0;
		
		cout << "\tFailed to find break" << endl;
	}
}

bool FusionSequence::FindBreakDeNovo(string& sequence, int& breakPos1, int& breakPos2, int& interLen, int& minCount) const
{
	int k = 21;
	int maxLength = 2 * mMaxFragmentLength;
	int reconstructOffset = k / 2;
	
	const int matchScore = 2;
	const int mismatchScore = -1;
	const int gapScore = -2;
	const int minAnchor = k;
	
	if (mCandidateSequences.size() == 0)
	{
		return false;
	}
	
	// Retreive candidate sequences
	SequenceVec breakSequences;
	IntegerVec sequenceEnds[2];
	for (SequenceMapConstIter sequenceIter = mCandidateSequences.begin(); sequenceIter != mCandidateSequences.end(); sequenceIter++)
	{
		int readID = sequenceIter->first;
		const string& readSeq = sequenceIter->second;
		
		breakSequences.push_back(readSeq);
		sequenceEnds[0].push_back(0);
		sequenceEnds[1].push_back(0);
	}
	
	// Retreive spanning sequences
	sequenceEnds[0] = IntegerVec(breakSequences.size(), 0);
	sequenceEnds[1] = IntegerVec(breakSequences.size(), 0);
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		for (SequenceMapConstIter sequenceIter = mSpanningSequences[clusterEnd].begin(); sequenceIter != mSpanningSequences[clusterEnd].end(); sequenceIter++)
		{
			int readID = sequenceIter->first;
			const string& readSeq = sequenceIter->second;
			
			breakSequences.push_back(readSeq);
			sequenceEnds[clusterEnd].push_back(1);
			sequenceEnds[1-clusterEnd].push_back(0);
		}
	}
	
	cout << "\tFinding Paths" << endl;
	
	// Find all paths between proximal sequences
	PathFinder pf(k,maxLength,false,2);
	pf.FindPath(breakSequences, sequenceEnds[0], sequenceEnds[1]);
	
	// No paths found
	if (pf.GetSequences().size() == 0)
	{
		return false;
	}
	
	// Use the first sequence path
	const string& deNovoSeq = pf.GetSequences()[0];
	int deNovoMinCount = pf.GetMinCounts()[0];
	
	// Output reference sequences to split alignments file
	breakOut << mDeNovoAlignSeq[0] << string(mNumBreakPadding, 'N') << mDeNovoAlignSeq[1] << endl;
	
	cout << "\tDoing Split Alignment" << endl;
	
	// Align the first sequence path
	SplitReadAligner splitReadAligner(matchScore, mismatchScore, gapScore, true, minAnchor * matchScore, mDeNovoAlignSeq[0], mDeNovoAlignSeq[1]);
	splitReadAligner.Align(deNovoSeq);
	
	// Retreive the split alignments
	SplitReadAlignVec splitAlignments;
	splitReadAligner.GetAlignments(splitAlignments, 0, true, true);
	
	// Check if a split alignment exists
	if (splitAlignments.size() == 0)
	{
		return false;
	}
	
	// Use the first split alignment
	const SplitReadAlignment& splitAlignment = splitAlignments[0];
	
	breakOut << "Denovo Prediction" << endl;
	
	// Output the alignment and calculate number bases on each side of the breakpoint
	int prevMatch = -1;
	for (int matchIndex = 0; matchIndex < splitAlignment.matches1.size(); matchIndex++)
	{
		const IntegerPair& match = splitAlignment.matches1[matchIndex];
		
		int numRefGap = match.first - prevMatch - 1;
		breakOut << string(numRefGap,(prevMatch == -1) ? ' ' : '-');
		
		breakOut << deNovoSeq[match.second];
		
		prevMatch = match.first;
	}
	
	int numRemaining = mDeNovoAlignSeq[0].length() - prevMatch - 1 + mNumBreakPadding;
	breakOut << string(numRemaining,'-');
	
	prevMatch = -1;
	for (int matchIndex = 0; matchIndex < splitAlignment.matches2.size(); matchIndex++)
	{
		const IntegerPair& match = splitAlignment.matches2[matchIndex];
		
		int numRefGap = match.first - prevMatch - 1;
		breakOut << string(numRefGap,'-');
		
		breakOut << deNovoSeq[match.second];
		
		prevMatch = match.first;
	}
	
	breakOut << endl;
	
	// Double check that there are enough matched nucleotides
	if (splitAlignment.matches1.size() <= reconstructOffset || splitAlignment.matches2.size() <= reconstructOffset)
	{
		return false;
	}

	// Reconstruct sequence half way into minimum anchor
	int reconstructPos1 = splitAlignment.matches1[reconstructOffset].first;
	int reconstructPos2 = splitAlignment.matches2[splitAlignment.matches2.size() - reconstructOffset - 1].first;
	string reconstructDenovoSeq = deNovoSeq.substr(reconstructOffset, deNovoSeq.size() - 2*reconstructOffset);
	
	// Construct break sequence
	breakSequence = mDeNovoAlignSeq[0].substr(0, reconstructPos1);
	breakSequence += reconstructDenovoSeq;
	breakSequence += mDeNovoAlignSeq[1].substr(reconstructPos2 + 1);
	
	if (mAlignStrand[0] == PlusStrand)
	{
		breakPos1 = mDeNovoAlignSeqPos[0].start + reconstructPos1 - 1;
	}
	else
	{
		breakPos1 = mDeNovoAlignSeqPos[0].end - reconstructPos1 + 1;
	}
	
	if (mAlignStrand[1] == PlusStrand)
	{
		breakPos2 = mDeNovoAlignSeqPos[1].start + reconstructPos2 + 1;
	}
	else
	{
		breakPos2 = mDeNovoAlignSeqPos[1].end - reconstructPos2 - 1;
	}
	
	interLen = reconstructDenovoSeq.size();
	
	// Return minimum kmer count
	minCount = deNovoMinCount;
		
	return true;
}

bool FusionSequence::FindBreakSplitReads(string& sequence, int& breakPos1, int& breakPos2, int& splitReadCount, double& splitPosAvg, double& splitMinAvg) const
{
	typedef map<IntegerPair,SplitReadAlignVec> SplitAlignmentsMap;
	typedef map<IntegerPair,IntegerVec> SplitReadIDsMap;
	typedef map<IntegerPair,StringVec> SplitSequencesMap;
	typedef map<IntegerPair,Counter> SplitScoreMap;
	typedef map<IntegerPair,Counter> SplitCountMap;
	typedef set<IntegerPair> SplitSet;
	
	const int matchScore = 2;
	const int mismatchScore = -1;
	const int gapScore = -2;
	const int minAnchor = 4;
	
	const int kmerHash = 7;
	
	if (mCandidateSequences.size() == 0)
	{
		return false;
	}
	
	cout << "\tSplit Prediction" << endl;
	
	SplitReadAligner splitReadAligner(matchScore, mismatchScore, gapScore, false, minAnchor * matchScore, mSplitAlignSeq[0], mSplitAlignSeq[1]);
	
	breakOut << mSplitAlignSeq[0] << string(mNumBreakPadding, 'N') << mSplitAlignSeq[1] << endl;
	
	SplitAlignmentsMap splitAlignmentsMap;
	SplitReadIDsMap splitReadIDsMap;
	SplitSequencesMap splitSequencesMap;
	SplitScoreMap splitScoreMap;
	SplitCountMap splitCountMap;

	cout << "\tSplit Read Alignment" << endl;
	
	for (SequenceMapConstIter sequenceIter = mCandidateSequences.begin(); sequenceIter != mCandidateSequences.end(); sequenceIter++)
	{
		int readID = sequenceIter->first;
		const string& readSeq = sequenceIter->second;
				
		splitReadAligner.Align(readSeq);
		
		SplitReadAlignVec splitAlignments;
		splitReadAligner.GetAlignments(splitAlignments, (int)((float)readSeq.length() * (float)matchScore * 0.90), true, false);
		
		SplitSet readSplits;
		
		for (int splitAlignIndex = 0; splitAlignIndex < splitAlignments.size(); splitAlignIndex++)
		{
			const SplitReadAlignment& splitAlignment = splitAlignments[splitAlignIndex];
			
			if (readSplits.find(splitAlignment.split) != readSplits.end())
			{
				continue;
			}
			
			readSplits.insert(splitAlignment.split);
			splitAlignmentsMap[splitAlignment.split].push_back(splitAlignment);
			splitReadIDsMap[splitAlignment.split].push_back(readID);
			splitSequencesMap[splitAlignment.split].push_back(readSeq);
			splitScoreMap[splitAlignment.split] += min(splitAlignment.score1,splitAlignment.score2);
			splitCountMap[splitAlignment.split] += 1;
		}
	}
	
	cout << "\tFinding splits" << endl;
	
	sequence = "";
	splitPosAvg = -1.0;
	splitMinAvg = -1.0;
	while (splitScoreMap.size() > 0)
	{
		int maxScore = -1;
		SplitScoreMap::iterator maxSplitScoreIter = splitScoreMap.end();
		for (SplitScoreMap::iterator splitScoreIter = splitScoreMap.begin(); splitScoreIter != splitScoreMap.end(); splitScoreIter++)
		{
			if (splitScoreIter->second > maxScore)
			{
				maxSplitScoreIter = splitScoreIter;
				maxScore = splitScoreIter->second;
			}
		}
		
		if (maxSplitScoreIter == splitScoreMap.end())
		{
			cerr << "Error: Unable to find max score split" << endl;
			break;
		}
		
		if (sequence == "")
		{
			IntegerPair bestSplit = maxSplitScoreIter->first;
			
			DebugCheck(bestSplit.first <= mSplitAlignSeq[0].length());
			DebugCheck(bestSplit.second + 1 < mSplitAlignSeq[1].length());
						
			string alignBreak1 = mSplitRemainderSeq[0] + mSplitAlignSeq[0].substr(0, bestSplit.first);
			string alignBreak2 = mSplitAlignSeq[1].substr(bestSplit.second + 1) + mSplitRemainderSeq[1];

			sequence = alignBreak1 + "|" + alignBreak2;

			if (mAlignStrand[0] == PlusStrand)
			{
				breakPos1 = mSplitAlignSeqPos[0].start + bestSplit.first - 1;
			}
			else
			{
				breakPos1 = mSplitAlignSeqPos[0].end - bestSplit.first + 1;				
			}
			
			if (mAlignStrand[1] == PlusStrand)
			{
				breakPos2 = mSplitAlignSeqPos[1].start + bestSplit.second + 1;
			}
			else
			{
				breakPos2 = mSplitAlignSeqPos[1].end - bestSplit.second - 1;
			}
			
			splitReadCount = splitCountMap[bestSplit];
		}
		
		const SplitReadAlignVec& splitAlignments = splitAlignmentsMap[maxSplitScoreIter->first];
		const IntegerVec& readIDs = splitReadIDsMap[maxSplitScoreIter->first];
		const StringVec& sequences = splitSequencesMap[maxSplitScoreIter->first];
		
		IntegerVec leftBaseCount;
		IntegerVec rightBaseCount;
		
		for (int splitAlignIndex = 0; splitAlignIndex < splitAlignments.size(); splitAlignIndex++)
		{
			const SplitReadAlignment& splitAlignment = splitAlignments[splitAlignIndex];
			ReadID readID;
			readID.id = readIDs[splitAlignIndex];
			const string& readSeq = sequences[splitAlignIndex];
			
			breakOut << readID.fragmentIndex << ((readID.readEnd == 0) ? "/1" : "/2") << endl;
			
			leftBaseCount.push_back(splitAlignment.matches1.size());
			rightBaseCount.push_back(splitAlignment.matches2.size());
			
			int prevMatch = -1;
			for (int matchIndex = 0; matchIndex < splitAlignment.matches1.size(); matchIndex++)
			{
				const IntegerPair& match = splitAlignment.matches1[matchIndex];
				
				int numRefGap = match.first - prevMatch - 1;
				breakOut << string(numRefGap,(prevMatch == -1) ? ' ' : '-');
				breakOut << readSeq[match.second];
				
				prevMatch = match.first;
			}
			
			int numRemaining = mSplitAlignSeq[0].length() - prevMatch - 1 + mNumBreakPadding;
			breakOut << string(numRemaining,'-');
			
			prevMatch = -1;
			for (int matchIndex = 0; matchIndex < splitAlignment.matches2.size(); matchIndex++)
			{
				const IntegerPair& match = splitAlignment.matches2[matchIndex];
				
				int numRefGap = match.first - prevMatch - 1;
				breakOut << string(numRefGap,'-');
				breakOut << readSeq[match.second];
				
				prevMatch = match.first;
			}
			
			breakOut << endl;
		}
		
		if (splitPosAvg == -1.0 || splitMinAvg == -1.0)
		{
			double posSum = 0.0;
			double minSum = 0.0;
			
			for (int splitIndex = 0; splitIndex < leftBaseCount.size(); splitIndex++)
			{
				double posRange = (double)(leftBaseCount[splitIndex] + rightBaseCount[splitIndex] - 2*minAnchor);
				double posValue = max(0, leftBaseCount[splitIndex] - minAnchor);

				double minRange = 0.5 * (double)(leftBaseCount[splitIndex] + rightBaseCount[splitIndex] - 2*minAnchor);
				double minValue = max(0, min(leftBaseCount[splitIndex] - minAnchor, rightBaseCount[splitIndex] - minAnchor));
				
				posSum += posValue / posRange;
				minSum += minValue / minRange;
			}
			
			splitPosAvg = posSum / leftBaseCount.size();
			splitMinAvg = minSum / leftBaseCount.size();
		}
		
		splitScoreMap.erase(maxSplitScoreIter);
	}
	
	return true;
}

void FusionSequence::CalculateBreakRegion(int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength)
{
	int alignRegionLength = alignEnd - alignStart + 1;
	
	// Push back the break region maximum read length nucleotides from 3' end of reads
	// or half way into the reads, whichever is minimum
	int pushBreakRegion = min(mMaxReadLength, (int)(0.5 * alignRegionLength));
	
	breakLength = mMaxFragmentLength - alignRegionLength - mMinReadLength + 2 * pushBreakRegion;
	
	if (strand == PlusStrand)
	{
		breakStart = alignEnd - pushBreakRegion + 1;
	}
	else
	{
		breakStart = alignStart + pushBreakRegion - 1;
	}
}

void FusionSequence::CalculateSplitMateRegion(int breakStart, int breakLength, int strand, int& mateRegionStart, int& mateRegionEnd)
{
	if (strand == PlusStrand)
	{
		int breakEnd = breakStart + breakLength - 1;

		mateRegionStart = breakStart - mMaxFragmentLength + mMinReadLength;
		mateRegionEnd = breakEnd - mMinFragmentLength + mMaxReadLength;
	}
	else
	{
		int breakEnd = breakStart - breakLength + 1;

		mateRegionStart = breakEnd + mMinFragmentLength - mMaxReadLength;
		mateRegionEnd = breakStart + mMaxFragmentLength - mMinReadLength;
	}
}

