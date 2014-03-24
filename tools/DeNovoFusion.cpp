/*
 *  DeNovoFusion.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "DeNovoFusion.h"
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

const int DeNovoFusion::mNumBreakPadding = 10;

DeNovoFusion::DeNovoFusion(const AlignmentIndex& discordant, const AlignmentIndex& anchored, const ReadIndex& reads, const FastaIndex& reference,
                           const ExonRegions& exonRegions, double fragmentLengthMean, double fragmentLengthStdDev, int minReadLength, int maxReadLength,
						   ostream& logFile)
: mDiscordant(discordant), mAnchored(anchored), mReads(reads), mReference(reference), mFragmentLengthMean(fragmentLengthMean),
  mExonRegions(exonRegions), mFragmentLengthStdDev(fragmentLengthStdDev), mMinReadLength(minReadLength), mMaxReadLength(maxReadLength),
  mLogFile(logFile)
{
	mMinFragmentLength = mFragmentLengthMean - 3 * mFragmentLengthStdDev;
	mMaxFragmentLength = mFragmentLengthMean + 3 * mFragmentLengthStdDev;
}

bool DeNovoFusion::Calculate(const LocationVec& alignPair, IntegerVec& breakPos, string& sequence, int& interLen, int& minCount)
{
	if (alignPair.size() != 2)
	{
		cerr << "Error: Incorrect input for DeNovoFusion::Calculate()" << endl;
		return false;
	}
	
	breakPos = IntegerVec(2);
	
	AddSpanningReads(alignPair);
	
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		const string& alignRefName = alignPair[clusterEnd].refName;
		int alignStrand = alignPair[clusterEnd].strand;
		int alignStart = alignPair[clusterEnd].start;
		int alignEnd = alignPair[clusterEnd].end;
		
		int refSeqStrand = (clusterEnd == 0) ? alignStrand : OtherStrand(alignStrand);
		bool revCompReads = (clusterEnd == 0) ? true : false;
		
		mDeNovoSeqStrand[clusterEnd] = refSeqStrand;
		
		int breakRegionFirst;
		int breakRegionLength;
		CalculateBreakRegion(alignStart, alignEnd, alignStrand, breakRegionFirst, breakRegionLength);

		int breakRegionStart;
		int breakRegionEnd;
		if (alignStrand == PlusStrand)
		{
			breakRegionStart = breakRegionFirst;
			breakRegionEnd = breakRegionFirst + breakRegionLength - 1;
		}
		else
		{
			breakRegionStart = breakRegionFirst - breakRegionLength + 1;
			breakRegionEnd = breakRegionFirst;
		}
		
		int alignBreakRegionStart = min(alignStart, breakRegionStart);
		int alignBreakRegionEnd = max(alignEnd, breakRegionEnd);
		
		mDeNovoAlignSeqStart[clusterEnd] = alignBreakRegionStart;
		mDeNovoAlignSeqLength[clusterEnd] = alignBreakRegionEnd - alignBreakRegionStart + 1;
		
		// Retrieve denovo reference alignment sequence
		mReference.Get(alignRefName, refSeqStrand, mDeNovoAlignSeqStart[clusterEnd], mDeNovoAlignSeqLength[clusterEnd], mDeNovoAlignSeq[clusterEnd]);
		
		// Find gene for this transcript
		const string& gene = mExonRegions.GetTranscriptGeneName(alignRefName);
		
		// Find all transcripts of the same gene
		const StringVec& geneTranscripts = mExonRegions.GetGeneTranscriptNames(gene);
		
		// Iterate through all transcripts of this gene
		for (StringVecConstIter geneTranscriptIter = geneTranscripts.begin(); geneTranscriptIter != geneTranscripts.end(); geneTranscriptIter++)
		{
			const string& remapTranscript = *geneTranscriptIter;
			
			// Remap break region start to this transcript if possible
			int remappedBreakRegionFirst;
			if (!mExonRegions.RemapTranscriptToTranscript(alignRefName, breakRegionFirst, remapTranscript, remappedBreakRegionFirst))
			{
				continue;
			}
			
			DebugCheck(remapTranscript != alignRefName || breakRegionFirst == remappedBreakRegionFirst);
			
			// Calculate the anchored mate alignment region in this transcript
			int mateRegionStart;
			int mateRegionEnd;
			CalculateSplitMateRegion(remappedBreakRegionFirst, breakRegionLength, alignStrand, mateRegionStart, mateRegionEnd);

			AddAnchoredReads(remapTranscript, alignStrand, revCompReads, mateRegionStart, mateRegionEnd);
		}
	}
	
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
	mLogFile << mDeNovoAlignSeq[0] << string(mNumBreakPadding, 'N') << mDeNovoAlignSeq[1] << endl;
	
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
	
	// Output the alignment and calculate number bases on each side of the breakpoint
	int prevMatch = -1;
	for (int matchIndex = 0; matchIndex < splitAlignment.matches1.size(); matchIndex++)
	{
		const IntegerPair& match = splitAlignment.matches1[matchIndex];
		
		int numRefGap = match.first - prevMatch - 1;
		mLogFile << string(numRefGap,(prevMatch == -1) ? ' ' : '-');
		
		mLogFile << deNovoSeq[match.second];
		
		prevMatch = match.first;
	}
	
	int numRemaining = mDeNovoAlignSeq[0].length() - prevMatch - 1 + mNumBreakPadding;
	mLogFile << string(numRemaining,'-');
	
	prevMatch = -1;
	for (int matchIndex = 0; matchIndex < splitAlignment.matches2.size(); matchIndex++)
	{
		const IntegerPair& match = splitAlignment.matches2[matchIndex];
		
		int numRefGap = match.first - prevMatch - 1;
		mLogFile << string(numRefGap,'-');
		
		mLogFile << deNovoSeq[match.second];
		
		prevMatch = match.first;
	}
	
	mLogFile << endl;
	
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
	sequence = mDeNovoAlignSeq[0].substr(0, reconstructPos1);
	sequence += reconstructDenovoSeq;
	sequence += mDeNovoAlignSeq[1].substr(reconstructPos2 + 1);
	
	if (mDeNovoSeqStrand[0] == PlusStrand)
	{
		breakPos[0] = mDeNovoAlignSeqStart[0] + reconstructPos1 - 1;
	}
	else
	{
		breakPos[0] = mDeNovoAlignSeqStart[0] + mDeNovoAlignSeqLength[0] - reconstructPos1;
	}
	
	if (mDeNovoSeqStrand[1] == PlusStrand)
	{
		breakPos[1] = mDeNovoAlignSeqStart[1] + reconstructPos2 + 1;
	}
	else
	{
		breakPos[1] = mDeNovoAlignSeqStart[1] + mDeNovoAlignSeqLength[1] - reconstructPos2 - 2;
	}
	
	// Return length of intervening sequence
	interLen = reconstructDenovoSeq.size();
	
	// Return minimum kmer count across path
	minCount = deNovoMinCount;
		
	return true;
}

void DeNovoFusion::AddAnchoredReads(const string& transcript, int strand, bool revComp, int start, int end)
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

void DeNovoFusion::AddSpanningReads(const LocationVec& alignPair)
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
			
			string spanningSeq = readSeqs[spanningReadID.readEnd];
			if (clusterEnd == 1)
			{
				ReverseComplement(spanningSeq);
			}
		
			mSpanningSequences[clusterEnd][spanningReadID.id] = spanningSeq;
		}
	}
}

void DeNovoFusion::CalculateBreakRegion(int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength)
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

void DeNovoFusion::CalculateSplitMateRegion(int breakStart, int breakLength, int strand, int& mateRegionStart, int& mateRegionEnd)
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

