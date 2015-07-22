/*
 *  SplitAlignment.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "SplitAlignment.h"
#include "DebugCheck.h"
#include "SplitReadAligner.h"
#include "Parsers.h"

#include <map>
#include <list>
#include <fstream>
#include <sstream>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

const int numBreakPadding = 10;
const int matchScore = 2;
const int mismatchScore = -1;
const int gapScore = -2;
const int minAnchor = 4;

bool SplitAlignmentTask::Initialize(int id, const LocationVec& alignPair, const FastaIndex& reference, 
		const ExonRegions& exonRegions, double fragmentLengthMean, 
		double fragmentLengthStdDev, int minReadLength, int maxReadLength)
{
	mFusionID = id;

	int minFragmentLength = (int)(fragmentLengthMean - 3 * fragmentLengthStdDev);
	int maxFragmentLength = (int)(fragmentLengthMean + 3 * fragmentLengthStdDev);

	if (alignPair.size() != 2)
	{
		cerr << "Error: Incorrect input for SplitAlignment::Calculate()" << endl;
		return false;
	}
	
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		const string& alignRefName = alignPair[clusterEnd].refName;
		int alignStrand = alignPair[clusterEnd].strand;
		int alignStart = alignPair[clusterEnd].start;
		int alignEnd = alignPair[clusterEnd].end;
		
		mAlignRefName[clusterEnd] = alignRefName;
		mAlignStrand[clusterEnd] = alignStrand;
		
		int refSeqStrand = (clusterEnd == 0) ? alignStrand : OtherStrand(alignStrand);
		bool revCompReads = (clusterEnd == 0) ? 1 : 0;
		
		int breakRegionStart;
		int breakRegionLength;
		CalculateBreakRegion(minReadLength, maxReadLength, maxFragmentLength, alignStart, alignEnd, alignStrand, breakRegionStart, breakRegionLength);
		
		mSplitSeqStrand[clusterEnd] = refSeqStrand;
		
		// Assumption: the break region we get from the paired end analysis will not involve a read that pushes more than 50%
		// into the breakpoint.  This is aligner dependent.
		if (alignStrand == PlusStrand)
		{
			mSplitAlignSeqStart[clusterEnd] = breakRegionStart - maxReadLength;
			mSplitAlignSeqLength[clusterEnd] = breakRegionLength + maxReadLength;
		}
		else
		{
			mSplitAlignSeqStart[clusterEnd] = breakRegionStart - breakRegionLength + 1;
			mSplitAlignSeqLength[clusterEnd] = breakRegionLength + maxReadLength;
		}
		
		// Retrieve split reference alignment sequence
		reference.Get(alignRefName, refSeqStrand, mSplitAlignSeqStart[clusterEnd], mSplitAlignSeqLength[clusterEnd], mSplitAlignSeq[clusterEnd]);
		
		// Remainder sequence for reconstruction
		mSplitRemainderSeq[clusterEnd] = "";
		if (alignStrand == PlusStrand)
		{
			if (alignStart < mSplitAlignSeqStart[clusterEnd])
			{
				int remainderStart = alignStart;
				int remainderEnd = mSplitAlignSeqStart[clusterEnd] - 1;
				int remainderLength = remainderEnd - remainderStart + 1;
				
				reference.Get(alignRefName, refSeqStrand, remainderStart, remainderLength, mSplitRemainderSeq[clusterEnd]);
			}
		}
		else
		{
			if (alignEnd > mSplitAlignSeqStart[clusterEnd] + mSplitAlignSeqLength[clusterEnd] - 1)
			{
				int remainderStart = mSplitAlignSeqStart[clusterEnd] + mSplitAlignSeqLength[clusterEnd];
				int remainderEnd = alignEnd;
				int remainderLength = remainderEnd - remainderStart + 1;
				
				reference.Get(alignRefName, refSeqStrand, remainderStart, remainderLength, mSplitRemainderSeq[clusterEnd]);
			}
		}
		
		string chromosome;
		int genomeAlignStrand;
		int genomeBreakRegionStart;
		string alignGene;
		string alignTranscript;
		if (ParseTranscriptID(alignRefName, alignGene, alignTranscript) && exonRegions.IsTranscript(alignTranscript))
		{
			// Remap break region start to the genome
			bool remapResult = exonRegions.RemapTranscriptToGenome(alignTranscript, alignStrand, breakRegionStart, chromosome, genomeAlignStrand, genomeBreakRegionStart);
			DebugCheck(remapResult);
		}
		else
		{
			chromosome = alignRefName;
			genomeAlignStrand = alignStrand;
			genomeBreakRegionStart = breakRegionStart;
		}
		
		// Find min and max range of mates
		int mateMin = minFragmentLength - breakRegionLength - maxReadLength + 1;
		int mateMax = maxFragmentLength - minReadLength;
		
		// Find start and end of mate region in the genome
		Region genomeMateRegion;
		if (genomeAlignStrand == PlusStrand)
		{
			genomeMateRegion.start = genomeBreakRegionStart - mateMax;
			genomeMateRegion.end = genomeBreakRegionStart - mateMin;
		}
		else
		{
			genomeMateRegion.start = genomeBreakRegionStart + mateMin;
			genomeMateRegion.end = genomeBreakRegionStart + mateMax;
		}
		
		mMateRegions[clusterEnd].push_back(Location());
		mMateRegions[clusterEnd].back().refName = chromosome;
		mMateRegions[clusterEnd].back().strand = genomeAlignStrand;
		mMateRegions[clusterEnd].back().start = genomeMateRegion.start;
		mMateRegions[clusterEnd].back().end = genomeMateRegion.end;
		
		// Find overlapping genes
		StringVec transcripts;
		exonRegions.GetRegionTranscripts(chromosome, genomeMateRegion, transcripts);
			
		for (StringVecConstIter transcriptIter = transcripts.begin(); transcriptIter != transcripts.end(); transcriptIter++)
		{
			const string& remapGene = exonRegions.GetTranscriptGene(*transcriptIter);
			const string& remapTranscript = *transcriptIter;
			
			string remapTranscriptID = remapGene + "|" + remapTranscript;
			
			int mateRegionStart;
			int mateRegionEnd;
			int remapAlignStrand;
			if (exonRegions.RemapThroughTranscript(remapTranscript, genomeBreakRegionStart, 1-genomeAlignStrand, mateMin, mateMax, remapAlignStrand, mateRegionStart, mateRegionEnd))
			{
				remapAlignStrand = 1-remapAlignStrand;
				
				mMateRegions[clusterEnd].push_back(Location());
				mMateRegions[clusterEnd].back().refName = remapTranscriptID;
				mMateRegions[clusterEnd].back().strand = remapAlignStrand;
				mMateRegions[clusterEnd].back().start = mateRegionStart;
				mMateRegions[clusterEnd].back().end = mateRegionEnd;
			}
		}
	}
	
	return true;
}

BinnedLocations::BinnedLocations(int binSpacing)
	: mBinSpacing(binSpacing)
{}

void BinnedLocations::Add(int id, const Location& location)
{
	int idx = (int)mIDs.size();

	Region region;
	region.start = location.start;
	region.end = location.end;

	mIDs.push_back(id);
	mRegions.push_back(region);

	int startBin = location.start / mBinSpacing;
	int endBin = location.end / mBinSpacing;
	
	for (int bin = startBin; bin <= endBin; bin++)
	{
		mBinned[location.strand][location.refName][bin].push_back(idx);
	}
}

void BinnedLocations::Overlapping(const RawAlignment& alignment, unordered_set<int>& ids) const
{
	unordered_map<string,unordered_map<int,IntegerVec> >::const_iterator findRefIter = mBinned[alignment.strand].find(alignment.reference);
	if (findRefIter != mBinned[alignment.strand].end())
	{
		int startBin = alignment.region.start / mBinSpacing;
		int endBin = alignment.region.end / mBinSpacing;
		
		for (int bin = startBin; bin <= endBin; bin++)
		{
			unordered_map<int,IntegerVec>::const_iterator findBinIter = findRefIter->second.find(bin);
			if (findBinIter != findRefIter->second.end())
			{
				for (IntegerVecConstIter iter = findBinIter->second.begin(); iter != findBinIter->second.end(); iter++)
				{
					int idx = *iter;

					int id = mIDs[idx];
					const Region& region = mRegions[idx];
					
					if (region.start <= alignment.region.end && region.end >= alignment.region.start)
					{
						ids.insert(id);
					}
				}
			}
		}
	}
}

SplitReadRealigner::SplitReadRealigner()
	: mBinnedMateRegions(2000),
	  mAligner(matchScore, mismatchScore, gapScore, false, minAnchor * matchScore)
{}

void SplitReadRealigner::AddTask(const SplitAlignmentTask& alignTask)
{
	mAlignTasks[alignTask.mFusionID] = alignTask;

	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		ClusterID clusterID;
		clusterID.clusterID = alignTask.mFusionID;
		clusterID.clusterEnd = clusterEnd;

		for (int mateRegionIndex = 0; mateRegionIndex < alignTask.mMateRegions[clusterEnd].size(); mateRegionIndex++)
		{
			mBinnedMateRegions.Add(clusterID.id, alignTask.mMateRegions[clusterEnd][mateRegionIndex]);
		}
	}
}

void SplitReadRealigner::AddReads(IReadStream* readStream)
{
	RawRead rawRead;
	while (readStream->GetNextRead(rawRead))
	{
		ReadID readID;
		readID.fragmentIndex = lexical_cast<int>(rawRead.fragment);
		readID.readEnd = rawRead.readEnd;

		mReads[readID.id] = rawRead.sequence;
	}
}

void SplitReadRealigner::DoAlignment(AlignmentStream* mateAlignments, ostream& alignmentsOut)
{
	unordered_map<int,unordered_set<IntegerPair> > candidateUnique;
	
	RawAlignment mateAlignment;
	while (mateAlignments->GetNextAlignment(mateAlignment))
	{
		unordered_set<int> overlapping;
		mBinnedMateRegions.Overlapping(mateAlignment, overlapping);
		
		for (unordered_set<int>::const_iterator olapIter = overlapping.begin(); olapIter != overlapping.end(); olapIter++)
		{
			ClusterID clusterID;
			clusterID.id = *olapIter;

			SplitReadInfo readInfo;
			readInfo.readID.fragmentIndex = lexical_cast<int>(mateAlignment.fragment);
			readInfo.readID.readEnd = (mateAlignment.readEnd == 0) ? 1 : 0;
			readInfo.revComp = (clusterID.clusterEnd == 0);
			
			string readSeq = mReads[readInfo.readID.id];
			if (readInfo.revComp)
			{
				ReverseComplement(readSeq);
			}
			
			if (candidateUnique[clusterID.clusterID].insert(IntegerPair(readInfo.readID.id,readInfo.revComp)).second)
			{
				vector<SplitAlignment> alignments = mAlignTasks[clusterID.clusterID].Align(mAligner, readInfo, readSeq);
				
				for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
				{
					alignments[alignmentIndex].WriteAlignment(alignmentsOut);
				}
			}
		}
	}
}

void SplitAlignment::WriteAlignment(ostream& out)
{
	out << fusionID << "\t";
	out << readInfo.readID.fragmentIndex << "\t";
	out << readInfo.readID.readEnd << "\t";
	out << readInfo.revComp << "\t";
	out << refSplit.first << "\t";
	out << refSplit.second << "\t";
	out << readSplit.first << "\t";
	out << readSplit.second << "\t";
	out << score << "\t";
	out << endl;
}

vector<SplitAlignment> SplitAlignment::ReadSortedAlignments(istream& in)
{
	vector<SplitAlignment> alignments;

	streampos prevPos = in.tellg();

	int firstFusionID = -1;
	bool firstEntry = true;

	string line;
	while (getline(in, line))
	{
		vector<string> fields;
		split(fields, line, is_any_of("\t"));
		
		if (fields.size() < 7)
		{
			cerr << "Error: Format error for candidate reads line:" << endl << line << endl;
			exit(1);
		}
		
		int fusionID = lexical_cast<int>(fields[0]);

		if (firstEntry)
		{
			firstFusionID = fusionID;
			firstEntry = false;
		}
		else if (fusionID != firstFusionID)
		{
			in.seekg(prevPos);
			break;
		}

		SplitAlignment alignment;

		alignment.fusionID = fusionID;
		alignment.readInfo.readID.fragmentIndex = lexical_cast<int>(fields[1]);
		alignment.readInfo.readID.readEnd = lexical_cast<int>(fields[2]);
		alignment.readInfo.revComp = lexical_cast<bool>(fields[3]);
		alignment.refSplit = make_pair(lexical_cast<int>(fields[4]), lexical_cast<int>(fields[5]));
		alignment.readSplit = make_pair(lexical_cast<int>(fields[6]), lexical_cast<int>(fields[7]));
		alignment.score = lexical_cast<int>(fields[8]);

		alignments.push_back(alignment);

		prevPos = in.tellg();
	}

	return alignments;
}

vector<SplitAlignment> SplitAlignmentTask::Align(SplitReadAligner& aligner, const SplitReadInfo& readInfo,
		const string& readSeq)
{
	vector<SplitAlignment> alignments;

	aligner.Align(readSeq, mSplitAlignSeq[0], mSplitAlignSeq[1]);
	
	SplitReadAlignVec splitAlignments;
	aligner.GetAlignments(splitAlignments, (int)((float)readSeq.length() * (float)matchScore * 0.90), true, false, GenerateAlignmentText);
	
	unordered_set<IntegerPair> readSplits;
	
	for (int splitAlignIndex = 0; splitAlignIndex < splitAlignments.size(); splitAlignIndex++)
	{
		const SplitReadAlignment& splitAlignment = splitAlignments[splitAlignIndex];
		
		if (readSplits.find(splitAlignment.refSplit) != readSplits.end())
		{
			continue;
		}
		
		readSplits.insert(splitAlignment.refSplit);

		SplitAlignment alignment;

		alignment.fusionID = mFusionID;
		alignment.readInfo = readInfo;
		alignment.refSplit = splitAlignment.refSplit;
		alignment.readSplit = splitAlignment.readSplit;
		alignment.score = min(splitAlignment.score1,splitAlignment.score2);
		
		if (GenerateAlignmentText)
		{
			stringstream alignmentText;
			
			alignmentText << readInfo.readID.fragmentIndex << ((readInfo.readID.readEnd == 0) ? "/1" : "/2") << endl;
			
			int prevMatch = -1;
			for (int matchIndex = 0; matchIndex < splitAlignment.matches1.size(); matchIndex++)
			{
				const IntegerPair& match = splitAlignment.matches1[matchIndex];
				
				int numRefGap = match.first - prevMatch - 1;
				alignmentText << string(numRefGap,(prevMatch == -1) ? ' ' : '-');
				alignmentText << readSeq[match.second];
				
				prevMatch = match.first;
			}
			
			int numRemaining = mSplitAlignSeq[0].length() - prevMatch - 1 + 1;
			alignmentText << string(numRemaining,'-');
			
			prevMatch = -1;
			for (int matchIndex = 0; matchIndex < splitAlignment.matches2.size(); matchIndex++)
			{
				const IntegerPair& match = splitAlignment.matches2[matchIndex];
				
				int numRefGap = match.first - prevMatch - 1;
				alignmentText << string(numRefGap,'-');
				alignmentText << readSeq[match.second];
				
				prevMatch = match.first;
			}
			
			alignmentText << endl;
			
			alignment.text = alignmentText.str();
		}

		alignments.push_back(alignment);
	}

	return alignments;
}

SplitAlignment SplitAlignmentTask::ReAlign(SplitReadAligner& aligner, const ReadIndex& readIndex,
		const SplitAlignment& alignment)
{
	string readSeq;
	readIndex.Find(alignment.readInfo.readID.fragmentIndex, alignment.readInfo.readID.readEnd, readSeq);

	if (alignment.readInfo.revComp)
	{
		ReverseComplement(readSeq);
	}

	vector<SplitAlignment> realignments = Align(aligner, alignment.readInfo, readSeq);

	for (int alignmentIndex = 0; alignmentIndex < realignments.size(); alignmentIndex++)
	{
		if (realignments[alignmentIndex].refSplit == alignment.refSplit)
		{
			return realignments[alignmentIndex];
		}
	}

	DebugCheck(false);
	return SplitAlignment();
}

vector<SplitAlignment> SplitAlignmentTask::ReAlign(SplitReadAligner& aligner,
		const ReadIndex& readIndex, const vector<SplitAlignment>& alignments)
{
	vector<SplitAlignment> realignments;

	for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
	{
		realignments.push_back(ReAlign(aligner, readIndex, alignments[alignmentIndex]));
	}

	return realignments;
}

BreakPrediction SplitAlignmentTask::Evaluate(const vector<SplitAlignment>& alignments) const
{
	BreakPrediction prediction;

	prediction.fusionID = mFusionID;
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		prediction.alignRefName[clusterEnd] = mAlignRefName[clusterEnd];
		prediction.alignStrand[clusterEnd] = mAlignStrand[clusterEnd];
	}

	prediction.sequence = "N";
	prediction.splitReadCount = 0;
	prediction.splitPosAvg = -1.0;
	prediction.splitMinAvg = -1.0;
		
	if (alignments.size() == 0)
	{
		return prediction;
	}

	unordered_map<IntegerPair,int> splitScore;
	for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
	{
		const IntegerPair& split = alignments[alignmentIndex].refSplit;
		int score = alignments[alignmentIndex].score;

		splitScore.insert(make_pair(split, 0)).first->second += score;
	}
	
	int maxScore = -1;
	IntegerPair bestSplit;
	for (unordered_map<IntegerPair,int>::iterator splitScoreIter = splitScore.begin();
			splitScoreIter != splitScore.end(); splitScoreIter++)
	{
		const IntegerPair& split = splitScoreIter->first;
		int score = splitScoreIter->second;

		if (score > maxScore)
		{
			bestSplit = split;
			maxScore = score;
		}
	}

	if (maxScore == -1)
	{
		cerr << "Error: Unable to find max score split" << endl;
		return prediction;
	}
	
	for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
	{
		const IntegerPair& split = alignments[alignmentIndex].refSplit;

		if (split == bestSplit)
		{
			prediction.alignments.push_back(alignments[alignmentIndex]);
		}
	}

	DebugCheck(bestSplit.first <= mSplitAlignSeq[0].length());
	DebugCheck(bestSplit.second + 1 < mSplitAlignSeq[1].length());
	
	string alignBreak1 = mSplitRemainderSeq[0] + mSplitAlignSeq[0].substr(0, bestSplit.first);
	string alignBreak2 = mSplitAlignSeq[1].substr(bestSplit.second + 1) + mSplitRemainderSeq[1];
	
	prediction.sequence = alignBreak1 + "|" + alignBreak2;
	
	if (mSplitSeqStrand[0] == PlusStrand)
	{
		prediction.breakPos[0] = mSplitAlignSeqStart[0] + bestSplit.first - 1;
	}
	else
	{
		prediction.breakPos[0] = mSplitAlignSeqStart[0] + mSplitAlignSeqLength[0] - bestSplit.first;				
	}
	
	if (mSplitSeqStrand[1] == PlusStrand)
	{
		prediction.breakPos[1] = mSplitAlignSeqStart[1] + bestSplit.second + 1;
	}
	else
	{
		prediction.breakPos[1] = mSplitAlignSeqStart[1] + mSplitAlignSeqLength[1] - bestSplit.second - 2;
	}
	
	double posSum = 0.0;
	double minSum = 0.0;
	
	for (int alignmentIndex = 0; alignmentIndex < prediction.alignments.size(); alignmentIndex++)
	{
		int leftBaseCount = prediction.alignments[alignmentIndex].readSplit.first;
		int rightBaseCount = prediction.alignments[alignmentIndex].readSplit.second;

		double posRange = (double)(leftBaseCount + rightBaseCount - 2*minAnchor);
		double posValue = max(0, leftBaseCount - minAnchor);
		
		double minRange = floor(0.5 * (double)(leftBaseCount + rightBaseCount - 2*minAnchor));
		double minValue = max(0, min(leftBaseCount - minAnchor, rightBaseCount - minAnchor));
		
		posSum += posValue / posRange;
		minSum += minValue / minRange;
	}
	
	prediction.splitReadCount = prediction.alignments.size();
	prediction.splitPosAvg = posSum / (double)prediction.alignments.size();
	prediction.splitMinAvg = minSum / prediction.alignments.size();
		
	return prediction;
}

void BreakPrediction::WriteSequence(ostream& out)
{
	out << fusionID << "\t";
	out << sequence << "\t";
	out << "0" << "\t";
	out << splitReadCount << "\t";
	out << splitPosAvg << "\t";
	out << splitMinAvg << endl;
}

void BreakPrediction::WriteBreak(ostream& out)
{
	for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
	{
		out << fusionID << "\t";
		out << clusterEnd << "\t";
		out << alignRefName[clusterEnd] << "\t";
		out << (alignStrand[clusterEnd] == PlusStrand ? "+" : "-") << "\t";
		out << breakPos[clusterEnd] << endl;
	}
}

void BreakPrediction::WriteAlignments(ostream& out)
{
	for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
	{
		alignments[alignmentIndex].WriteAlignment(out);
	}
}

void SplitAlignmentTask::WriteAlignText(ostream& out, const vector<SplitAlignment>& alignments)
{
	out << mFusionID << endl;
	out << mSplitAlignSeq[0] << "|" << mSplitAlignSeq[1] << endl;
	for (int alignmentIndex = 0; alignmentIndex < alignments.size(); alignmentIndex++)
	{
		DebugCheck(!alignments[alignmentIndex].text.empty());
		out << alignments[alignmentIndex].text;
	}
}

void SplitAlignmentTask::CalculateBreakRegion(int minReadLength, int maxReadLength, int maxFragmentLength, int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength)
{
	int alignRegionLength = alignEnd - alignStart + 1;
	
	// Push back the break region maximum read length nucleotides from 3' end of reads
	// or half way into the reads, whichever is minimum
	int pushBreakRegion = min(maxReadLength, (int)(0.5 * alignRegionLength));
	
	breakLength = maxFragmentLength - alignRegionLength - minReadLength + 2 * pushBreakRegion;
	
	if (strand == PlusStrand)
	{
		breakStart = alignEnd - pushBreakRegion + 1;
	}
	else
	{
		breakStart = alignStart + pushBreakRegion - 1;
	}
}

unordered_map<int,SplitAlignmentTask> CreateTasks(const string& referenceFasta,
		const string& exonRegionsFilename, double fragmentLengthMean, double fragmentLengthStdDev,
		int minReadLength, int maxReadLength, const LocationVecMap& fusionRegions)
{
	FastaIndex reference;
	ExonRegions exonRegions;
	
	reference.Open(referenceFasta);
	
	ifstream exonRegionsFile(exonRegionsFilename.c_str());
	if (!exonRegionsFile.good() || !exonRegions.Read(exonRegionsFile))
	{
		cerr << "Error: Unable to read exon regions file " << exonRegionsFilename << endl;
		exit(1);
	}
	exonRegionsFile.close();

	unordered_map<int,SplitAlignmentTask> alignTasks;

	for (LocationVecMapConstIter pairIter = fusionRegions.begin(); pairIter != fusionRegions.end(); pairIter++)
	{
		int id = pairIter->first;
		const LocationVec& alignRegionPair = pairIter->second;

		alignTasks[id].Initialize(id, alignRegionPair, reference, exonRegions,
			fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength);
	}

	return alignTasks;
}


