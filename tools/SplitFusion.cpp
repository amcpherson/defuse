/*
 *  SplitFusion.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "SplitFusion.h"
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

const int SplitFusion::mNumBreakPadding = 10;

SplitFusion::SplitFusion(const AlignmentIndex& discordant, const AlignmentIndex& anchored, const ReadIndex& reads, const FastaIndex& reference,
						 const ExonRegions& exonRegions, double fragmentLengthMean, double fragmentLengthStdDev, int minReadLength, int maxReadLength,
						 ostream& logFile)
: mDiscordant(discordant), mAnchored(anchored), mReads(reads), mReference(reference), mFragmentLengthMean(fragmentLengthMean),
  mExonRegions(exonRegions), mFragmentLengthStdDev(fragmentLengthStdDev), mMinReadLength(minReadLength), mMaxReadLength(maxReadLength),
  mLogFile(logFile)
{
	mMinFragmentLength = mFragmentLengthMean - 3 * mFragmentLengthStdDev;
	mMaxFragmentLength = mFragmentLengthMean + 3 * mFragmentLengthStdDev;
}

bool SplitFusion::Calculate(const LocationVec& alignPair, IntegerVec& breakPos, string& sequence, int& splitReadCount, double& splitPosAvg, double& splitMinAvg)
{
	if (alignPair.size() != 2)
	{
		cerr << "Error: Incorrect input for SplitFusion::Calculate()" << endl;
		return false;
	}
	
	breakPos = IntegerVec(2);
	
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
		
		mSplitSeqStrand[clusterEnd] = refSeqStrand;
		
		// Assumption: the break region we get from the paired end analysis will not involve a read that pushes more than 50%
		// into the breakpoint.  This is aligner dependent.
		if (alignStrand == PlusStrand)
		{
			mSplitAlignSeqStart[clusterEnd] = breakRegionStart - mMaxReadLength;
			mSplitAlignSeqLength[clusterEnd] = breakRegionLength + mMaxReadLength;
		}
		else
		{
			mSplitAlignSeqStart[clusterEnd] = breakRegionStart - breakRegionLength + 1;
			mSplitAlignSeqLength[clusterEnd] = breakRegionLength + mMaxReadLength;
		}
		
		// Retrieve split reference alignment sequence
		mReference.Get(alignRefName, refSeqStrand, mSplitAlignSeqStart[clusterEnd], mSplitAlignSeqLength[clusterEnd], mSplitAlignSeq[clusterEnd]);
		
		// Remainder sequence for reconstruction
		mSplitRemainderSeq[clusterEnd] = "";
		if (alignStrand == PlusStrand)
		{
			if (alignStart < mSplitAlignSeqStart[clusterEnd])
			{
				int remainderStart = alignStart;
				int remainderEnd = mSplitAlignSeqStart[clusterEnd] - 1;
				int remainderLength = remainderEnd - remainderStart + 1;
				mReference.Get(alignRefName, refSeqStrand, remainderStart, remainderLength, mSplitRemainderSeq[clusterEnd]);
			}
		}
		else
		{
			if (alignEnd > mSplitAlignSeqStart[clusterEnd] + mSplitAlignSeqLength[clusterEnd] - 1)
			{
				int remainderStart = mSplitAlignSeqStart[clusterEnd] + mSplitAlignSeqLength[clusterEnd];
				int remainderEnd = alignEnd;
				int remainderLength = remainderEnd - remainderStart + 1;
				
				mReference.Get(alignRefName, refSeqStrand, remainderStart, remainderLength, mSplitRemainderSeq[clusterEnd]);
			}
		}
		
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
	
	SplitReadAligner splitReadAligner(matchScore, mismatchScore, gapScore, false, minAnchor * matchScore, mSplitAlignSeq[0], mSplitAlignSeq[1]);
	
	mLogFile << mSplitAlignSeq[0] << string(mNumBreakPadding, 'N') << mSplitAlignSeq[1] << endl;
	
	SplitAlignmentsMap splitAlignmentsMap;
	SplitReadIDsMap splitReadIDsMap;
	SplitSequencesMap splitSequencesMap;
	SplitScoreMap splitScoreMap;
	SplitCountMap splitCountMap;

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

			if (mSplitSeqStrand[0] == PlusStrand)
			{
				breakPos[0] = mSplitAlignSeqStart[0] + bestSplit.first - 1;
			}
			else
			{
				breakPos[0] = mSplitAlignSeqStart[0] + mSplitAlignSeqLength[0] - bestSplit.first;				
			}
			
			if (mSplitSeqStrand[1] == PlusStrand)
			{
				breakPos[1] = mSplitAlignSeqStart[1] + bestSplit.second + 1;
			}
			else
			{
				breakPos[1] = mSplitAlignSeqStart[1] + mSplitAlignSeqLength[1] - bestSplit.second - 2;
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
			
			mLogFile << readID.fragmentIndex << ((readID.readEnd == 0) ? "/1" : "/2") << endl;
			
			leftBaseCount.push_back(splitAlignment.matches1.size());
			rightBaseCount.push_back(splitAlignment.matches2.size());
			
			int prevMatch = -1;
			for (int matchIndex = 0; matchIndex < splitAlignment.matches1.size(); matchIndex++)
			{
				const IntegerPair& match = splitAlignment.matches1[matchIndex];
				
				int numRefGap = match.first - prevMatch - 1;
				mLogFile << string(numRefGap,(prevMatch == -1) ? ' ' : '-');
				mLogFile << readSeq[match.second];
				
				prevMatch = match.first;
			}
			
			int numRemaining = mSplitAlignSeq[0].length() - prevMatch - 1 + mNumBreakPadding;
			mLogFile << string(numRemaining,'-');
			
			prevMatch = -1;
			for (int matchIndex = 0; matchIndex < splitAlignment.matches2.size(); matchIndex++)
			{
				const IntegerPair& match = splitAlignment.matches2[matchIndex];
				
				int numRefGap = match.first - prevMatch - 1;
				mLogFile << string(numRefGap,'-');
				mLogFile << readSeq[match.second];
				
				prevMatch = match.first;
			}
			
			mLogFile << endl;
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

void SplitFusion::AddAnchoredReads(const string& transcript, int strand, bool revComp, int start, int end)
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

void SplitFusion::CalculateBreakRegion(int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength)
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

void SplitFusion::CalculateSplitMateRegion(int breakStart, int breakLength, int strand, int& mateRegionStart, int& mateRegionEnd)
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

