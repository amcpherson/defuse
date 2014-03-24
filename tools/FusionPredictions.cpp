/*
 *  FusionPredictions.cpp
 *  findbreaks
 *
 *  Created by Andrew McPherson on 28/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Alignments.h"
#include "AlignmentStream.h"
#include "FusionPredictions.h"
#include "GenomeIndex.h"
#include "PathFinder.h"
#include "SplitReadAligner.h"
#include "DebugCheck.h"

#include <map>
#include <list>
#include <fstream>
#include <sstream>



double CalculateCoverageStatistic(const vector<int>& points)
{
	unordered_set<int> uniquepoints(points.begin(), points.end());
	
	return uniquepoints.size();
}

inline bool ClusterOrder(const IntegerPairVec& vec1, const IntegerPairVec& vec2)
{
	return vec1.size() > vec2.size();
}

const int FusionPredictions::mNumBreakPadding = 10;

FusionPredictions::FusionPredictions(const Alignments& alignments, const Parameters& parameters, const IntegerPairTable& clusters, const ExonRegions& exonRegions)
: mAlignments(alignments), mParameters(parameters), mClusters(clusters), mExonRegions(exonRegions)
{
	sort(mClusters.begin(), mClusters.end(), ClusterOrder);
	InitializeClusterBreaks();
}

void FusionPredictions::CalculateBreakRegion(int strand, const Region& alignRegion, const Parameters& parameters, Region& breakRegion) const
{
	if (strand == PlusStrand)
	{
		breakRegion.start = alignRegion.start + parameters.GetMinReadAlignment();
		breakRegion.end = breakRegion.start + parameters.GetMaxBreakRegionSize() - 1;		
	}
	else if (strand == MinusStrand)
	{
		breakRegion.end = alignRegion.end - parameters.GetMinReadAlignment();		
		breakRegion.start = breakRegion.end - parameters.GetMaxBreakRegionSize() + 1;
	}
	else
	{
		DebugCheck(strand == PlusStrand || strand == MinusStrand);
	}
}

void FusionPredictions::FilterLowCoverage()
{
	const double threshold = 1.2;

	double maxReadLength = mParameters.GetMaxReadLength();
	double maxFragmentLength = mParameters.GetMaxFragmentSize();	
	double thresholdCoverageStat = (double)mParameters.GetMaxReadLength() * threshold / (double)(maxFragmentLength - maxReadLength);

	IntegerPairTable newClusters;
	
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		float coverageStat1;
		float coverageStat2;
		CalculateClusterCoverageStat(clusterIndex, coverageStat1, coverageStat2);

		if (coverageStat1 >= thresholdCoverageStat && coverageStat2 >= thresholdCoverageStat)
		{
			newClusters.push_back(mClusters[clusterIndex]);
		}
	}

	swap(mClusters,newClusters);
}

struct ClusterBreak
{
	int clusterIndex;
	int end;
	RefStrand refStrand;
	Region region;
	
	typedef unsigned int ChromosomeType;
	
	unsigned int GetChromosome()
	{
		return refStrand.id;
	}
	
	const Region& GetRegion()
	{
		return region;
	}
};


void FusionPredictions::ClearClusterReads()
{
	mClusterReads.resize(mClusters.size());
	mClusterReadRevComp.resize(mClusters.size());
	
	mReadSequences.clear();
}

void FusionPredictions::FindSpanningClusterReads()
{
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{		
		IntegerVec readIDs1;
		IntegerVec readIDs2;
		CalculateClusterReads(clusterIndex, readIDs1, readIDs2);
		
		for (int readIDPairIndex = 0; readIDPairIndex < readIDs1.size(); readIDPairIndex++)
		{
			int readID1 = readIDs1[readIDPairIndex];
			int readID2 = readIDs2[readIDPairIndex];
			
			mClusterReads[clusterIndex].push_back(readID1);
			mClusterReads[clusterIndex].push_back(readID2);
			
			mClusterReadRevComp[clusterIndex].push_back(0);
			mClusterReadRevComp[clusterIndex].push_back(1);
			
			mReadSequences[readID1] = "";
			mReadSequences[readID2] = "";
		}
	}
}

void FusionPredictions::InitializeClusterBreaks()
{
	mClusterBreaks.clear();

	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{		
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			RefStrand refStrand;
			GetClusterRefStrand(clusterIndex, clusterEnd, refStrand);
			
			int breakRegionStart;
			int breakRegionLength;
			CalculateClusterBreakRegion(clusterIndex, clusterEnd, breakRegionStart, breakRegionLength);
						
			const string& transcript = Alignments::ReferenceNames.Get(refStrand.referenceIndex);
			
			// Find gene for this transcript
			const string& gene = mExonRegions.GetTranscriptGeneName(transcript);
			
			// Find all transcripts of the same gene
			const StringVec& geneTranscripts = mExonRegions.GetGeneTranscriptNames(gene);
			
			// Iterate through all transcripts of this gene
			for (StringVecConstIter geneTranscriptIter = geneTranscripts.begin(); geneTranscriptIter != geneTranscripts.end(); geneTranscriptIter++)
			{
				const string& remapTranscript = *geneTranscriptIter;
				
				// Remap break region start to this transcript if possible
				int remappedBreakRegionStart;
				if (!mExonRegions.RemapTranscriptToTranscript(transcript, breakRegionStart, remapTranscript, remappedBreakRegionStart))
				{
					continue;
				}
				
				DebugCheck(remapTranscript != transcript || breakRegionStart == remappedBreakRegionStart);					
				
				// Calculate break region in this transcript
				// Also calculate union of align region and break region in this transcript
				Region breakRegion;
				if (refStrand.strand == PlusStrand)
				{
					breakRegion.start = remappedBreakRegionStart;
					breakRegion.end = remappedBreakRegionStart + breakRegionLength - 1;
				}
				else
				{
					breakRegion.start = remappedBreakRegionStart - breakRegionLength + 1;
					breakRegion.end = remappedBreakRegionStart;
				}
				
				// Add cluster break for this transcript
				ClusterBreak clusterBreak;
				clusterBreak.clusterIndex = clusterIndex;
				clusterBreak.end = clusterEnd;
				clusterBreak.refStrand.referenceIndex = Alignments::ReferenceNames.Index(remapTranscript);
				clusterBreak.refStrand.strand = refStrand.strand;
				clusterBreak.region = breakRegion;
				mClusterBreaks.push_back(clusterBreak);
			}
		}
	}
}


void FusionPredictions::FindHangingClusterReads(IAlignmentStream* alignments)
{
	GenomeIndex<ClusterBreak> clusterBreakGenomeIndex;
	for (int clusterBreakIndex = 0; clusterBreakIndex < mClusterBreaks.size(); clusterBreakIndex++)
	{
		clusterBreakGenomeIndex.Insert(&mClusterBreaks[clusterBreakIndex]);
	}
		
	RawAlignment rawAlignment;
	while (alignments->GetNextAlignment(rawAlignment))
	{
		ReadID readID;
		readID.fragmentIndex = Alignments::FragmentNames.Index(rawAlignment.fragment);
		readID.readEnd = rawAlignment.readEnd;
		
		if (mReadMisMatches.find(readID.id) == mReadMisMatches.end())
		{
			mReadMisMatches[readID.id] = rawAlignment.numMis;
		}
		else
		{
			mReadMisMatches[readID.id] = min(mReadMisMatches[readID.id], rawAlignment.numMis);
		}
		
		RefStrand refStrand;
		refStrand.referenceIndex = Alignments::ReferenceNames.Index(rawAlignment.reference);
		refStrand.strand = rawAlignment.strand;
		
		// Assuming reads are same size for each end of the fragment
		int readLength = rawAlignment.sequence.length();
		
		// Calculate region of possible alignment for other read
		Region otherReadAlignment;
		if (refStrand.strand == PlusStrand)
		{
			otherReadAlignment.start = rawAlignment.region.start + mParameters.GetMinFragmentSize() - readLength;
			otherReadAlignment.end = rawAlignment.region.start + mParameters.GetMaxFragmentSize() - 1;
		}
		else
		{
			otherReadAlignment.start = rawAlignment.region.end - mParameters.GetMaxFragmentSize() + 1;
			otherReadAlignment.end = rawAlignment.region.end - mParameters.GetMinFragmentSize() + readLength;
		}		
		
		set<ClusterBreak*> overlapClusterBreaks;
		clusterBreakGenomeIndex.GetOverlapping(refStrand.id, otherReadAlignment, overlapClusterBreaks);
		
		for (set<ClusterBreak*>::const_iterator overlapIter = overlapClusterBreaks.begin(); overlapIter != overlapClusterBreaks.end(); overlapIter++)
		{
			const ClusterBreak* clusterBreak = *overlapIter;
			
			ReadID otherReadID;
			otherReadID.fragmentIndex = readID.fragmentIndex;
			otherReadID.readEnd = (readID.readEnd == 0) ? 1 : 0;
			
			if (clusterBreak->clusterIndex >= mClusters.size())
			{
				cerr << "Error: corruption with break region index" << endl;
				exit(1);
			}
			
			int revComp = (clusterBreak->end == 0) ? 1 : 0;
			
			mClusterReads[clusterBreak->clusterIndex].push_back(otherReadID.id);
			mClusterReadRevComp[clusterBreak->clusterIndex].push_back(revComp);
			
			mReadSequences[readID.id] = "";
			mReadSequences[otherReadID.id] = "";
		}
	}
}

void FusionPredictions::RetrieveReadSequences(IReadStream* reads)
{
	RawRead rawRead;
	while (reads->GetNextRead(rawRead))
	{
		ReadID readID;
		readID.fragmentIndex = Alignments::FragmentNames.Index(rawRead.fragment);
		readID.readEnd = rawRead.readEnd;		
		
		SequenceMapIter seqIter = mReadSequences.find(readID.id);
		if (seqIter != mReadSequences.end())
		{
			seqIter->second = rawRead.sequence;
		}
	}
}

void FusionPredictions::RetrieveClusterSequences(Reference& reference)
{
	LocationVec locations;
	
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		Location location1;
		Location location2;
		
		CalculateClusterBreakRegions(clusterIndex, location1.refStrand, location1.region, location2.refStrand, location2.region);
		
		// Assumption: the break region we get from the paired end analysis will not involve a read that pushes more than 50%
		// into the breakpoint.  This is aligner dependent.
		location1.region.start -= (int)(1.0 * mParameters.GetMaxReadLength());
		location1.region.end += (int)(1.0 * mParameters.GetMaxReadLength());
		
		location2.region.start -= (int)(1.0 * mParameters.GetMaxReadLength());
		location2.region.end += (int)(1.0 * mParameters.GetMaxReadLength());
		
		locations.push_back(location1);
		locations.push_back(location2);
		
		mClusterSequencePos[0].push_back(location1);
		mClusterSequencePos[1].push_back(location2);
	}
	
	SequenceVec sequences;
	reference.GetSequences(locations, sequences);
	
	mClusterSequences[0].resize(mClusters.size());
	mClusterSequences[1].resize(mClusters.size());
	
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		RefStrand refStrand1;
		RefStrand refStrand2;
		
		Region alignRegion1;
		Region alignRegion2;
		
		CalculateClusterAlignRegions(clusterIndex, refStrand1, alignRegion1, refStrand2, alignRegion2);
		
		int seq1Index = clusterIndex * 2;
		int seq2Index = clusterIndex * 2 + 1;
		
		if (seq1Index >= sequences.size() && seq2Index >= sequences.size())
		{
			cerr << "Error: problem finding reference sequences" << endl;
			exit(1);
		}
		
		string seq1 = sequences[seq1Index];
		string seq2 = sequences[seq2Index];
		
		if (seq1.size() == 0 || seq2.size() == 0)
		{
			continue;
		}
		
		DebugCheck(refStrand1.id == mClusterSequencePos[0][clusterIndex].refStrand.id);
		DebugCheck(refStrand2.id == mClusterSequencePos[1][clusterIndex].refStrand.id);

		DebugCheck(seq1.size() == mClusterSequencePos[0][clusterIndex].region.end - mClusterSequencePos[0][clusterIndex].region.start + 1);		
		DebugCheck(seq2.size() == mClusterSequencePos[1][clusterIndex].region.end - mClusterSequencePos[1][clusterIndex].region.start + 1);		

		if (refStrand1.strand == MinusStrand)
		{
			ReverseComplement(seq1);
		}
		
		if (refStrand2.strand == PlusStrand)
		{
			ReverseComplement(seq2);
		}
		
		mClusterSequences[0][clusterIndex] = seq1;
		mClusterSequences[1][clusterIndex] = seq2;
	}
}

void FusionPredictions::RetrievePartnerSequences(Reference& reference)
{
	LocationVec locations;
	
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		Location location1;
		Location location2;
		
		CalculateClusterAlignRegions(clusterIndex, location1.refStrand, location1.region, location2.refStrand, location2.region);
				
		locations.push_back(location1);
		locations.push_back(location2);

		mPartnerSequencePos[0].push_back(location1);
		mPartnerSequencePos[1].push_back(location2);
	}
	
	SequenceVec sequences;
	reference.GetSequences(locations, sequences);
	
	mPartnerSequences[0].resize(mClusters.size());
	mPartnerSequences[1].resize(mClusters.size());
	
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		RefStrand refStrand1;
		RefStrand refStrand2;
		
		Region alignRegion1;
		Region alignRegion2;
		
		CalculateClusterAlignRegions(clusterIndex, refStrand1, alignRegion1, refStrand2, alignRegion2);
		
		int seq1Index = clusterIndex * 2;
		int seq2Index = clusterIndex * 2 + 1;
		
		if (seq1Index >= sequences.size() && seq2Index >= sequences.size())
		{
			cerr << "Error: problem finding reference sequences" << endl;
			exit(1);
		}
		
		string seq1 = sequences[seq1Index];
		string seq2 = sequences[seq2Index];
		
		if (refStrand1.strand == MinusStrand)
		{
			ReverseComplement(seq1);
		}
		
		if (refStrand2.strand == PlusStrand)
		{
			ReverseComplement(seq2);
		}
		
		mPartnerSequences[0][clusterIndex] = seq1;
		mPartnerSequences[1][clusterIndex] = seq2;
	}
}

void FusionPredictions::FindBreaks()
{
	mClusterBreakSequences.resize(mClusters.size());
	mSplitReadCount.resize(mClusters.size());
	mClusterFragmentMean.resize(mClusters.size());
	mSplitPosAvg.resize(mClusters.size());
	mSplitMinAvg.resize(mClusters.size());

	// Generate output
	string splitAlignOutputFilename = mParameters.GetOutputDirectory() + "/split_alignments.txt";
	ofstream splitAlignOutputFile(splitAlignOutputFilename.c_str());

	if (!splitAlignOutputFile.good())
	{
		cerr << "Error: Unable to open " << splitAlignOutputFilename << endl;
		exit(1);
	}

	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		if (mClusterSequences[0][clusterIndex].size() == 0 || mClusterSequences[1][clusterIndex].size() == 0)
		{
			continue;
		}
		
		cout << "Cluster Break: " << clusterIndex << " " << mClusters[clusterIndex].size() << " " << mClusterReads[clusterIndex].size() << endl;
		
		string deNovoBreakSequence;
		int deNovoMinCount;
		double deNovoFragmentMean;

		bool deNovoResult = FindBreakDeNovo(clusterIndex, splitAlignOutputFile, deNovoBreakSequence, deNovoMinCount, deNovoFragmentMean);

		string splitBreakSequence;
		int splitReadCount;
		double splitPosAvg;
		double splitMinAvg;
		double splitFragmentMean;
		
		bool splitResult = false;
		if (mClusterReads[clusterIndex].size() < 20000)
		{
			splitResult = FindBreakSplitReads(clusterIndex, splitAlignOutputFile, splitBreakSequence, splitReadCount, splitPosAvg, splitMinAvg, splitFragmentMean);
		}

		// Select the prediction with the sample fragment length mean closest to the overall fragment length mean
		double deNovoFragmentMeanDiff = abs(mParameters.GetFragmentSizeMean() - deNovoFragmentMean);
		double splitFragmentMeanDiff = abs(mParameters.GetFragmentSizeMean() - splitFragmentMean);
		
		if (deNovoResult && (!splitResult || deNovoFragmentMeanDiff < splitFragmentMeanDiff))
		{
			mClusterBreakSequences[clusterIndex] = deNovoBreakSequence;
			mSplitReadCount[clusterIndex] = deNovoMinCount;
			mClusterFragmentMean[clusterIndex] = deNovoFragmentMean;
			mSplitPosAvg[clusterIndex] = 1.0;
			mSplitMinAvg[clusterIndex] = 1.0;
		}
		else if (splitResult)
		{
			mClusterBreakSequences[clusterIndex] = splitBreakSequence;
			mSplitReadCount[clusterIndex] = splitReadCount;
			mClusterFragmentMean[clusterIndex] = splitFragmentMean;
			mSplitPosAvg[clusterIndex] = splitPosAvg;
			mSplitMinAvg[clusterIndex] = splitMinAvg;
		}
		else
		{
			mClusterBreakSequences[clusterIndex] = "";
			mSplitReadCount[clusterIndex] = 0;
			mClusterFragmentMean[clusterIndex] = 0.0;
			mSplitPosAvg[clusterIndex] = 0.0;
			mSplitMinAvg[clusterIndex] = 0.0;
			
			cout << "\tFailed to find break" << endl;
		}
	}
}

bool FusionPredictions::FindBreakDeNovo(int clusterIndex, ostream& breakOut, string& breakSequence, int& minCount, double& fragmentMean) const
{
	int k = 21;
	int maxLength = 2 * mParameters.GetMaxFragmentSize();
	int reconstructOffset = k / 2;
	
	const int matchScore = 2;
	const int mismatchScore = -1;
	const int gapScore = -2;
	const int minAnchor = k;
	
	if (mClusterSequences[0][clusterIndex].size() == 0 || mClusterSequences[1][clusterIndex].size() == 0)
	{
		return false;
	}
	
	IntegerVec readIDs1;
	IntegerVec readIDs2;
	CalculateClusterReads(clusterIndex, readIDs1, readIDs2);
	
	unordered_set<int> leftReadIDs(readIDs1.begin(), readIDs1.end());
	unordered_set<int> rightReadIDs(readIDs2.begin(), readIDs2.end());
	
	// Retreive all sequences
	SequenceVec breakSequences;
	IntegerVec leftSequence;
	IntegerVec rightSequence;
	unordered_set<int> uniqueReadIDs;
	for (int clusterReadIndex = 0; clusterReadIndex < mClusterReads[clusterIndex].size(); clusterReadIndex++)
	{
		ReadID readID;
		readID.id = mClusterReads[clusterIndex][clusterReadIndex];
		
		if (uniqueReadIDs.find(readID.id) != uniqueReadIDs.end())
		{
			continue;
		}
		uniqueReadIDs.insert(readID.id);
		
		int revComp = mClusterReadRevComp[clusterIndex][clusterReadIndex];
		
		SequenceMapConstIter seqIter = mReadSequences.find(readID.id);
		if (seqIter == mReadSequences.end() || seqIter->second == "")
		{
			cerr << "Error: Unable to find sequence for " << Alignments::FragmentNames.Get(readID.fragmentIndex) << endl;
			continue;
		}
		
		string readSeq = seqIter->second;
		
		if (revComp)
		{
			ReverseComplement(readSeq);
		}
		
		breakSequences.push_back(readSeq);
		
		int left = (leftReadIDs.find(readID.id) != leftReadIDs.end()) ? 1 : 0;
		int right = (rightReadIDs.find(readID.id) != rightReadIDs.end()) ? 1 : 0;
		
		leftSequence.push_back(left);
		rightSequence.push_back(right);
	}
	
	cout << "\tFinding Paths" << endl;
	
	// Find all paths between proximal sequences
	PathFinder pf(k,maxLength,false,2);
	pf.FindPath(breakSequences, leftSequence, rightSequence);
	
	// No paths found
	if (pf.GetSequences().size() == 0)
	{
		return false;
	}
	
	// Use the first sequence path
	const string& deNovoSeq = pf.GetSequences()[0];
	int deNovoMinCount = pf.GetMinCounts()[0];
	
	// Create combined partner breakpoint sequences
	string partnerBreakpointSequence[2];	
	CalculateAlignBreakSequence(clusterIndex, partnerBreakpointSequence[0], partnerBreakpointSequence[1]);
	
	// Output reference sequences to split alignments file
	breakOut << "cluster: " << clusterIndex << endl;
	breakOut << partnerBreakpointSequence[0] << string(mNumBreakPadding, 'N') << partnerBreakpointSequence[1] << endl;
	
	cout << "\tDoing Split Alignment" << endl;
	
	// Align the first sequence path
	SplitReadAligner splitReadAligner(matchScore, mismatchScore, gapScore, true, minAnchor * matchScore, partnerBreakpointSequence[0], partnerBreakpointSequence[1]);
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
	
	int numRemaining = partnerBreakpointSequence[0].length() - prevMatch - 1 + mNumBreakPadding;
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
	
	// Calculate fragment size mean
	const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
	
	if (alignmentPairs.size() < 1)
	{
		return false;
	}

	const Location& partner1Pos = mPartnerSequencePos[0][clusterIndex];
	const Location& partner2Pos = mPartnerSequencePos[1][clusterIndex];
		
	double fragmentSizeSum = 0;
	double fragmentSizeSumSq = 0;
	
	// Double check that there are enough matched nucleotides
	if (splitAlignment.matches1.size() <= reconstructOffset || splitAlignment.matches2.size() <= reconstructOffset)
	{
		return false;
	}

	// Reconstruct sequence half way into minimum anchor
	int reconstructPos1 = splitAlignment.matches1[reconstructOffset].first;
	int reconstructPos2 = splitAlignment.matches2[splitAlignment.matches2.size() - reconstructOffset - 1].first;
	string reconstructDenovoSeq = deNovoSeq.substr(reconstructOffset, deNovoSeq.size() - 2*reconstructOffset);
	
	for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
	{
		int end1AlignIndex = alignPairIter->first;
		int end2AlignIndex = alignPairIter->second;
		
		const Alignment& end1Alignment = mAlignments.GetAlignment(end1AlignIndex);
		const Alignment& end2Alignment = mAlignments.GetAlignment(end2AlignIndex);
		
		DebugCheck(end1Alignment.refStrand.id == partner1Pos.refStrand.id);
		DebugCheck(end2Alignment.refStrand.id == partner2Pos.refStrand.id);
		
		int localStart;
		if (partner1Pos.refStrand.strand == PlusStrand)
		{
			localStart = end1Alignment.region.start - partner1Pos.region.start + 1;
		}
		else
		{
			localStart = partner1Pos.region.end - end1Alignment.region.end + 1;
		}
		
		int partnerAlign1Size = reconstructPos1 + 1 - localStart;
		
		int localEnd;
		if (partner2Pos.refStrand.strand == MinusStrand)
		{
			localEnd = partner2Pos.region.end - end2Alignment.region.end + 1;
		}
		else
		{
			localEnd = end2Alignment.region.start - partner2Pos.region.start + 1;
		}
		
		int partnerAlign2Size = partnerBreakpointSequence[1].size() - localEnd - reconstructPos2;
		
		int fragmentSize = reconstructDenovoSeq.size() + partnerAlign1Size + partnerAlign2Size;
		
		fragmentSizeSum += fragmentSize;
		fragmentSizeSumSq += fragmentSize * fragmentSize;
	}
	
	// Construct break sequence
	breakSequence = partnerBreakpointSequence[0].substr(0, reconstructPos1);
	breakSequence += reconstructDenovoSeq;
	breakSequence += partnerBreakpointSequence[1].substr(reconstructPos2 + 1);			
	
	// Return fragment length mean
	fragmentMean = fragmentSizeSum / (double)alignmentPairs.size();
	
	// Return minimum kmer count
	minCount = deNovoMinCount;
		
	return true;
}

bool FusionPredictions::FindBreakSplitReads(int clusterIndex, ostream& breakOut, string& breakSequence, int& splitReadCount, double& splitPosAvg, double& splitMinAvg, double& fragmentMean) const
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
	
	if (mClusterSequences[0][clusterIndex].size() == 0 || mClusterSequences[1][clusterIndex].size() == 0)
	{
		return false;
	}
	
	cout << "\tSplit Prediction" << endl;
	cout << "\tHashing reference" << endl;
	
	// Create a hash of kmers in each partner sequence
	unordered_set<string> clusterSequenceHash1;
	for (int start = 0; start <= (int)mClusterSequences[0][clusterIndex].size() - kmerHash; start++)
	{
		clusterSequenceHash1.insert(mClusterSequences[0][clusterIndex].substr(start,kmerHash));
	}
	
	unordered_set<string> clusterSequenceHash2;
	for (int start = 0; start <= (int)mClusterSequences[1][clusterIndex].size() - kmerHash; start++)
	{
		clusterSequenceHash2.insert(mClusterSequences[1][clusterIndex].substr(start,kmerHash));
	}
	
	SplitReadAligner splitReadAligner(matchScore, mismatchScore, gapScore, false, minAnchor * matchScore, mClusterSequences[0][clusterIndex], mClusterSequences[1][clusterIndex]);
	
	breakOut << "cluster: " << clusterIndex << endl;
	breakOut << mClusterSequences[0][clusterIndex] << string(mNumBreakPadding, 'N') << mClusterSequences[1][clusterIndex] << endl;
	
	const IntegerVec& clusterReads = mClusterReads[clusterIndex];
	const IntegerVec& clusterReadRevComp = mClusterReadRevComp[clusterIndex];
	
	StringVec clusterReadSequences;
	IntegerVec clusterReadIDs;
	
	cout << "\tCollecting reads" << endl;
	
	unordered_set<int> uniqueReadIDs;
	for (int clusterReadIndex = 0; clusterReadIndex < clusterReads.size(); clusterReadIndex++)
	{
		ReadID readID;
		readID.id = clusterReads[clusterReadIndex];
		
		if (uniqueReadIDs.find(readID.id) != uniqueReadIDs.end())
		{
			continue;
		}
		uniqueReadIDs.insert(readID.id);
		
		IntegerMapConstIter mismatchesIter = mReadMisMatches.find(readID.id);
		if (mismatchesIter != mReadMisMatches.end())
		{
			int bestNumMismatches = mismatchesIter->second;
			if (bestNumMismatches <= 1)
			{
				continue;
			}
		}
		
		int revComp = clusterReadRevComp[clusterReadIndex];
		
		SequenceMapConstIter seqIter = mReadSequences.find(readID.id);
		if (seqIter == mReadSequences.end())
		{
			cerr << "Error: Unable to find sequence for " << Alignments::FragmentNames.Get(readID.fragmentIndex) << endl;
			continue;
		}
		
		string readSeq = seqIter->second;
		
		if (revComp)
		{
			ReverseComplement(readSeq);
		}
		
		clusterReadSequences.push_back(readSeq);
		clusterReadIDs.push_back(readID.id);
	}
	
	DebugCheck(clusterReads.size() == clusterReadRevComp.size());
	
	SplitAlignmentsMap splitAlignmentsMap;
	SplitReadIDsMap splitReadIDsMap;
	SplitSequencesMap splitSequencesMap;
	SplitScoreMap splitScoreMap;
	SplitCountMap splitCountMap;

	cout << "\tSplit Read Alignment" << endl;
	
	for (int clusterReadSequenceIndex = 0; clusterReadSequenceIndex < clusterReadSequences.size(); clusterReadSequenceIndex++)
	{
		const string& readSeq = clusterReadSequences[clusterReadSequenceIndex];
		int readID = clusterReadIDs[clusterReadSequenceIndex];
		
		bool found1 = false;
		bool found2 = false;
		for (int start = 0; start <= readSeq.size() - kmerHash; start++)
		{
			if (clusterSequenceHash1.find(readSeq.substr(start,kmerHash)) != clusterSequenceHash1.end())
			{
				found1 = true;
			}
			
			if (clusterSequenceHash2.find(readSeq.substr(start,kmerHash)) != clusterSequenceHash2.end())
			{
				found2 = true;
			}
		}
		
		if (!found1 || !found2)
		{
			continue;
		}
		
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
	
	breakSequence = "";
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
		
		if (breakSequence == "")
		{
			IntegerPair bestSplit = maxSplitScoreIter->first;
			
			DebugCheck(maxSplitScoreIter->first.first <= mClusterSequences[0][clusterIndex].length());
			DebugCheck(maxSplitScoreIter->first.second + 1 < mClusterSequences[1][clusterIndex].length());
			
			const Location& cluster1Pos = mClusterSequencePos[0][clusterIndex];
			const Location& cluster2Pos = mClusterSequencePos[1][clusterIndex];
			
			string alignBreak1;
			string alignBreak2;
			CalculateAlignBreakSequence(clusterIndex, bestSplit.first, bestSplit.second, alignBreak1, alignBreak2);
			breakSequence = alignBreak1 + "|" + alignBreak2;

			splitReadCount = splitCountMap[bestSplit];
			
			const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
			
			int breakPos1;
			if (cluster1Pos.refStrand.strand == PlusStrand)
			{
				breakPos1 = cluster1Pos.region.start + bestSplit.first;
			}
			else
			{
				breakPos1 = cluster1Pos.region.end - bestSplit.first;
			}
			
			int breakPos2;
			if (cluster2Pos.refStrand.strand == MinusStrand)
			{
				breakPos2 = cluster2Pos.region.start + bestSplit.second;
			}
			else
			{
				breakPos2 = cluster2Pos.region.end - bestSplit.second;
			}
			
			double fragmentSizeSum = 0;
			double fragmentSizeSumSq = 0;
			
			for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
			{
				int end1AlignIndex = alignPairIter->first;
				int end2AlignIndex = alignPairIter->second;
				
				const Alignment& end1Alignment = mAlignments.GetAlignment(end1AlignIndex);
				const Alignment& end2Alignment = mAlignments.GetAlignment(end2AlignIndex);
				
				DebugCheck(end1Alignment.refStrand.id == cluster1Pos.refStrand.id);
				DebugCheck(end2Alignment.refStrand.id == cluster2Pos.refStrand.id);
				
				int end1Size;
				if (cluster1Pos.refStrand.strand == PlusStrand)
				{
					end1Size = breakPos1 - end1Alignment.region.start;
				}
				else
				{
					end1Size = end1Alignment.region.end - breakPos1;
				}
				
				int end2Size;
				if (cluster2Pos.refStrand.strand == MinusStrand)
				{
					end2Size = end2Alignment.region.end - breakPos2;
				}
				else
				{
					end2Size = breakPos2 - end2Alignment.region.start;
				}
				
				int fragmentSize = end1Size + end2Size;
				
				fragmentSizeSum += fragmentSize;
				fragmentSizeSumSq += fragmentSize * fragmentSize;
			}
			
			fragmentMean = fragmentSizeSum / (double)alignmentPairs.size();
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
			
			breakOut << Alignments::FragmentNames.Get(readID.fragmentIndex) << ((readID.readEnd == 0) ? "/1" : "/2") << endl;
			
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
			
			int numRemaining = mClusterSequences[0][clusterIndex].length() - prevMatch - 1 + mNumBreakPadding;
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

bool FusionPredictions::GetClusterRefStrand(int clusterIndex, int clusterEnd, RefStrand& refStrand)
{
	const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
	
	if (alignmentPairs.size() == 0)
	{
		return false;
	}
	
	for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
	{
		int endAlignIndices[2];
		endAlignIndices[0] = alignPairIter->first;
		endAlignIndices[1] = alignPairIter->second;
		
		int endAlignIndex = endAlignIndices[clusterEnd];
		
		const Alignment& alignment = mAlignments.GetAlignment(endAlignIndex);
		
		DebugCheck(alignPairIter == alignmentPairs.begin() || refStrand.id == alignment.refStrand.id);
		
		refStrand = alignment.refStrand;
	}
	
	return true;
}

bool FusionPredictions::CalculateClusterAlignRegion(int clusterIndex, int clusterEnd, int& alignRegionStart, int& alignRegionLength)
{
	const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
	
	if (alignmentPairs.size() == 0)
	{
		return false;
	}
	
	Region clusterAlignRegion;
	
	for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
	{
		int endAlignIndices[2];
		endAlignIndices[0] = alignPairIter->first;
		endAlignIndices[1] = alignPairIter->second;
		
		int endAlignIndex = endAlignIndices[clusterEnd];
		
		const Alignment& alignment = mAlignments.GetAlignment(endAlignIndex);
				
		Region breakRegion;
		CalculateBreakRegion(alignment.refStrand.strand, alignment.region, mParameters, breakRegion);
		
		if (alignPairIter == alignmentPairs.begin())
		{
			clusterAlignRegion = alignment.region;
		}
		else
		{
			clusterAlignRegion.start = min(clusterAlignRegion.start, alignment.region.start);
			clusterAlignRegion.end = max(clusterAlignRegion.end, alignment.region.end);
		}
	}
	
	return true;
}

bool FusionPredictions::CalculateClusterBreakRegion(int clusterIndex, int end, int& start, int& length)
{
	const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
	
	if (alignmentPairs.size() == 0)
	{
		return false;
	}
	
	RefStrand refStrand;
	Region clusterAlignRegion;
	Region clusterBreakRegion;
	
	for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
	{
		int endAlignIndices[2];
		endAlignIndices[0] = alignPairIter->first;
		endAlignIndices[1] = alignPairIter->second;
		
		int endAlignIndex = endAlignIndices[end];
		
		const Alignment& alignment = mAlignments.GetAlignment(endAlignIndex);
		
		refStrand = alignment.refStrand;
		
		Region breakRegion;
		CalculateBreakRegion(alignment.refStrand.strand, alignment.region, mParameters, breakRegion);
		
		if (alignPairIter == alignmentPairs.begin())
		{
			clusterAlignRegion = alignment.region;
			clusterBreakRegion = breakRegion;
		}
		else
		{
			clusterAlignRegion.start = min(clusterAlignRegion.start, alignment.region.start);
			clusterAlignRegion.end = max(clusterAlignRegion.end, alignment.region.end);

			clusterBreakRegion.start = max(clusterBreakRegion.start, breakRegion.start);
			clusterBreakRegion.end = min(clusterBreakRegion.end, breakRegion.end);
		}
	}
	
	// Push back the break region maximum read length nucleotides from 3' end of reads
	// or half way into the reads, whichever is minimum
	int pushBreakRegion = min(mParameters.GetMaxReadLength(), (int)(0.5 * (clusterAlignRegion.end - clusterAlignRegion.start + 1)));

	if (refStrand.strand == PlusStrand)
	{
		clusterBreakRegion.start -= pushBreakRegion;
		
		start = clusterBreakRegion.start;
		length = clusterBreakRegion.end - clusterBreakRegion.start + 1;
	}
	else
	{
		clusterBreakRegion.end += pushBreakRegion;
		
		start = clusterBreakRegion.end;
		length = clusterBreakRegion.end - clusterBreakRegion.start + 1;
	}
		
	return true;
}

bool FusionPredictions::CalculateClusterBreakRegions(int clusterIndex, RefStrand& refStrand1, Region& breakRegion1, RefStrand& refStrand2, Region& breakRegion2) const
{
	const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
	
	if (alignmentPairs.size() == 0)
	{
		return false;
	}
	
	for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
	{
		int end1AlignIndex = alignPairIter->first;
		int end2AlignIndex = alignPairIter->second;
		
		const Alignment& end1Alignment = mAlignments.GetAlignment(end1AlignIndex);
		const Alignment& end2Alignment = mAlignments.GetAlignment(end2AlignIndex);
		
		refStrand1 = end1Alignment.refStrand;
		refStrand2 = end2Alignment.refStrand;
		
		Region end1BreakRegion;
		Region end2BreakRegion;
		
		CalculateBreakRegion(end1Alignment.refStrand.strand, end1Alignment.region, mParameters, end1BreakRegion);
		CalculateBreakRegion(end2Alignment.refStrand.strand, end2Alignment.region, mParameters, end2BreakRegion);
		
		if (alignPairIter == alignmentPairs.begin())
		{
			breakRegion1 = end1BreakRegion;
			breakRegion2 = end2BreakRegion;
		}
		else
		{
			breakRegion1.start = max(breakRegion1.start, end1BreakRegion.start);
			breakRegion2.start = max(breakRegion2.start, end2BreakRegion.start);
			
			breakRegion1.end = min(breakRegion1.end, end1BreakRegion.end);
			breakRegion2.end = min(breakRegion2.end, end2BreakRegion.end);
		}
	}
	
	return true;
}

bool FusionPredictions::CalculateClusterAlignRegions(int clusterIndex, RefStrand& refStrand1, Region& alignRegion1, RefStrand& refStrand2, Region& alignRegion2) const
{
	const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
	
	if (alignmentPairs.size() == 0)
	{
		return false;
	}
	
	for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
	{
		int end1AlignIndex = alignPairIter->first;
		int end2AlignIndex = alignPairIter->second;
		
		const Alignment& end1Alignment = mAlignments.GetAlignment(end1AlignIndex);
		const Alignment& end2Alignment = mAlignments.GetAlignment(end2AlignIndex);
		
		refStrand1 = end1Alignment.refStrand;
		refStrand2 = end2Alignment.refStrand;
		
		if (alignPairIter == alignmentPairs.begin())
		{
			alignRegion1 = end1Alignment.region;
			alignRegion2 = end2Alignment.region;
		}
		else
		{
			alignRegion1.start = min(alignRegion1.start, end1Alignment.region.start);
			alignRegion2.start = min(alignRegion2.start, end2Alignment.region.start);
			
			alignRegion1.end = max(alignRegion1.end, end1Alignment.region.end);
			alignRegion2.end = max(alignRegion2.end, end2Alignment.region.end);
		}
	}
	
	return true;
}

void FusionPredictions::CalculateAlignBreakSequence(int clusterIndex, int breakPos1, int breakPos2, string& alignBreak1, string& alignBreak2) const
{	
	DebugCheck(breakPos1 >= 0 && breakPos1 <= (int)mClusterSequences[0][clusterIndex].size());
	DebugCheck(breakPos2 >= -1 && breakPos2 <= (int)mClusterSequences[1][clusterIndex].size() - 1);
	
	const Location& cluster1Pos = mClusterSequencePos[0][clusterIndex];
	const Location& cluster2Pos = mClusterSequencePos[1][clusterIndex];
	
	const Location& partner1Pos = mPartnerSequencePos[0][clusterIndex];
	const Location& partner2Pos = mPartnerSequencePos[1][clusterIndex];
	
	int partner1Size;
	if (cluster1Pos.refStrand.strand == PlusStrand)
	{
		partner1Size = cluster1Pos.region.start - partner1Pos.region.start;
	}
	else
	{
		partner1Size = partner1Pos.region.end - cluster1Pos.region.end;
	}
	
	int cluster1Start = 0;
	int cluster1Size = breakPos1;
	if (partner1Size < 0)
	{
		cluster1Start = -partner1Size;
		cluster1Size = max(0, breakPos1 + partner1Size);
		partner1Size = 0;
	}
	
	DebugCheck(partner1Size >= 0 && partner1Size <= mPartnerSequences[0][clusterIndex].size());
	
	int partner2Start;
	if (cluster2Pos.refStrand.strand == PlusStrand)
	{
		partner2Start = partner2Pos.region.end - cluster2Pos.region.start + 1;
	}
	else
	{
		partner2Start = cluster2Pos.region.end - partner2Pos.region.start + 1;
	}
	
	int cluster2Size = mClusterSequences[1][clusterIndex].size() - breakPos2 - 1;
	if (partner2Start > mPartnerSequences[1][clusterIndex].size())
	{
		cluster2Size = max(0, cluster2Size - (partner2Start - (int)mPartnerSequences[1][clusterIndex].size()));
		partner2Start = mPartnerSequences[1][clusterIndex].size();
	}
	
	DebugCheck(partner2Start >= 0 && partner2Start <= mPartnerSequences[1][clusterIndex].size());
	
	alignBreak1 = mPartnerSequences[0][clusterIndex].substr(0, partner1Size);
	alignBreak1 += mClusterSequences[0][clusterIndex].substr(cluster1Start, cluster1Size);
	
	alignBreak2 = mClusterSequences[1][clusterIndex].substr(breakPos2 + 1, cluster2Size);
	alignBreak2 += mPartnerSequences[1][clusterIndex].substr(partner2Start);
}

void FusionPredictions::CalculateAlignBreakSequence(int clusterIndex, string& alignBreak1, string& alignBreak2) const
{	
	CalculateAlignBreakSequence(clusterIndex, mClusterSequences[0][clusterIndex].size(), -1, alignBreak1, alignBreak2);
}

bool FusionPredictions::CalculateClusterReads(int clusterIndex, IntegerVec& readIDs1, IntegerVec& readIDs2) const
{
	const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
	
	if (alignmentPairs.size() == 0)
	{
		return false;
	}
	
	for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
	{
		int end1AlignIndex = alignPairIter->first;
		int end2AlignIndex = alignPairIter->second;
		
		const Alignment& end1Alignment = mAlignments.GetAlignment(end1AlignIndex);
		const Alignment& end2Alignment = mAlignments.GetAlignment(end2AlignIndex);
		
		readIDs1.push_back(end1Alignment.readID.id);
		readIDs2.push_back(end2Alignment.readID.id);
	}
	
	return true;
}

bool FusionPredictions::CalculateClusterCoverageStat(int clusterIndex, float& coverageStat1, float& coverageStat2) const
{
	vector<int> pointsEnd1;
	vector<int> pointsEnd2;
	for (IntegerPairVecConstIter alignPairIter = mClusters[clusterIndex].begin(); alignPairIter != mClusters[clusterIndex].end(); alignPairIter++)
	{
		int end1AlignIndex = alignPairIter->first;
		int end2AlignIndex = alignPairIter->second;
		
		const Alignment& end1Alignment = mAlignments.GetAlignment(end1AlignIndex);
		const Alignment& end2Alignment = mAlignments.GetAlignment(end2AlignIndex);
		
		for (int pos = end1Alignment.region.start; pos <= end1Alignment.region.end; pos++)
		{
			pointsEnd1.push_back(pos);
		}
		for (int pos = end2Alignment.region.start; pos <= end2Alignment.region.end; pos++)
		{
			pointsEnd2.push_back(pos);
		}
	}
	
	double maxReadLength = mParameters.GetMaxReadLength();
	double maxFragmentLength = mParameters.GetMaxFragmentSize();
	
	coverageStat1 = CalculateCoverageStatistic(pointsEnd1) / (double)(maxFragmentLength - maxReadLength);
	coverageStat2 = CalculateCoverageStatistic(pointsEnd2) / (double)(maxFragmentLength - maxReadLength);
	
	return true;
}

int FusionPredictions::GetNumClusters()
{
	return mClusters.size();
}

void FusionPredictions::PrintClusters()
{
	// Generate output
	string clusterOutputFilename = mParameters.GetOutputDirectory() + "/fusions.txt";
	ofstream clusterOutputFile(clusterOutputFilename.c_str());
	
	if (!clusterOutputFile.good())
	{
		cerr << "Error: Unable to open " << clusterOutputFilename << endl;
	}

	// Generate output
	string clusterFragmentsOutputFilename = mParameters.GetOutputDirectory() + "/fusion_fragments.txt";
	ofstream clusterFragmentsOutputFile(clusterFragmentsOutputFilename.c_str());
	
	if (!clusterFragmentsOutputFile.good())
	{
		cerr << "Error: Unable to open " << clusterFragmentsOutputFilename << endl;
		exit(1);
	}

	for (int clusterIndex = 0; clusterIndex != mClusters.size(); clusterIndex++)
	{
		const IntegerPairVec& alignmentPairs = mClusters[clusterIndex];
		
		if (alignmentPairs.size() < 1)
		{
			continue;
		}
				
		RefStrand refStrand1;
		RefStrand refStrand2;
		
		Region alignRegion1;
		Region alignRegion2;
		
		Region breakRegion1;
		Region breakRegion2;
		
		float coverageStat1;
		float coverageStat2;
		
		CalculateClusterAlignRegions(clusterIndex, refStrand1, alignRegion1, refStrand2, alignRegion2);
		CalculateClusterBreakRegions(clusterIndex, refStrand1, breakRegion1, refStrand2, breakRegion2);
		CalculateClusterCoverageStat(clusterIndex, coverageStat1, coverageStat2);

		const string& referenceName1 = Alignments::ReferenceNames.Get(refStrand1.referenceIndex);
		const string& referenceName2 = Alignments::ReferenceNames.Get(refStrand2.referenceIndex);
		
		clusterOutputFile << clusterIndex << "\t";
		clusterOutputFile << alignmentPairs.size() << "\t";
		clusterOutputFile << referenceName1 << "\t";
		clusterOutputFile << ((refStrand1.strand == PlusStrand) ? '+' : '-') << "\t";
		clusterOutputFile << alignRegion1.start << "\t" << alignRegion1.end << "\t";
		clusterOutputFile << breakRegion1.start << "\t" << breakRegion1.end << "\t";
		clusterOutputFile << coverageStat1 << "\t";
		clusterOutputFile << mPartnerSequences[0][clusterIndex] << "\t";
		clusterOutputFile << referenceName2 << "\t";
		clusterOutputFile << ((refStrand2.strand == PlusStrand) ? '+' : '-') << "\t";
		clusterOutputFile << alignRegion2.start << "\t" << alignRegion2.end << "\t";
		clusterOutputFile << breakRegion2.start << "\t" << breakRegion2.end << "\t";
		clusterOutputFile << coverageStat2 << "\t";
		clusterOutputFile << mPartnerSequences[1][clusterIndex] << "\t";
		clusterOutputFile << mClusterBreakSequences[clusterIndex] << "\t";
		clusterOutputFile << mSplitReadCount[clusterIndex] << "\t";
		clusterOutputFile << mSplitPosAvg[clusterIndex] << "\t";
		clusterOutputFile << mSplitMinAvg[clusterIndex] << "\t";
		clusterOutputFile << mClusterFragmentMean[clusterIndex] << "\t";		
		clusterOutputFile << endl;
		
		for (IntegerPairVecConstIter alignPairIter = alignmentPairs.begin(); alignPairIter != alignmentPairs.end(); alignPairIter++)
		{
			int end1AlignIndex = alignPairIter->first;
			int end2AlignIndex = alignPairIter->second;
			
			const Alignment& end1Alignment = mAlignments.GetAlignment(end1AlignIndex);
			const Alignment& end2Alignment = mAlignments.GetAlignment(end2AlignIndex);
			
			clusterFragmentsOutputFile << clusterIndex << "\t";
			clusterFragmentsOutputFile << Alignments::FragmentNames.Get(end1Alignment.readID.fragmentIndex) << "\t";
			clusterFragmentsOutputFile << mReadSequences[end1Alignment.readID.id] << "\t";
			clusterFragmentsOutputFile << referenceName1 << "\t";
			clusterFragmentsOutputFile << ((end1Alignment.refStrand.strand == PlusStrand) ? '+' : '-') << "\t";
			clusterFragmentsOutputFile << end1Alignment.region.start << "\t";
			clusterFragmentsOutputFile << end1Alignment.region.end << "\t";
			clusterFragmentsOutputFile << mReadSequences[end2Alignment.readID.id] << "\t";
			clusterFragmentsOutputFile << referenceName2 << "\t";
			clusterFragmentsOutputFile << ((end2Alignment.refStrand.strand == PlusStrand) ? '+' : '-') << "\t";
			clusterFragmentsOutputFile << end2Alignment.region.start << "\t";
			clusterFragmentsOutputFile << end2Alignment.region.end << "\t";
			clusterFragmentsOutputFile << endl;
		}		
	}	
}


