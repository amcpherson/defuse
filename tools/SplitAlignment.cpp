/*
 *  SplitAlignment.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "SplitAlignment.h"
#include "DebugCheck.h"
#include "SplitReadAligner.h"

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

bool SplitAlignment::FindCandidates(const LocationVec& alignPair, const AlignmentIndex& discordant, const AlignmentIndex& anchored, const FastaIndex& reference, 
							        const ExonRegions& exonRegions, double fragmentLengthMean, double fragmentLengthStdDev, int minReadLength, int maxReadLength)
{
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
		
		AddAnchoredReads(discordant, anchored, chromosome, genomeAlignStrand, revCompReads, genomeMateRegion.start, genomeMateRegion.end);
		
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
				AddAnchoredReads(discordant, anchored, remapTranscriptID, remapAlignStrand, revCompReads, mateRegionStart, mateRegionEnd);
			}
		}
	}
	
	return true;
}

void SplitAlignment::WriteCandidateReads(ostream& out, SplitAlignmentMap& splitAlignments)
{
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;

		for (int candidateIndex = 0; candidateIndex < splitAlignment.mCandidateReadID.size(); candidateIndex++)
		{
			ReadID readID;
			readID.id = splitAlignment.mCandidateReadID[candidateIndex];

			out << id << "\t";
			out << readID.fragmentIndex << "\t";
			out << readID.readEnd << "\t";
			out << splitAlignment.mCandidateRevComp[candidateIndex] << "\t";
			out << endl;
		}
	}
}

void SplitAlignment::WriteCandidateRegions(ostream& out, SplitAlignmentMap& splitAlignments)
{
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;
		
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			out << id << "\t";
			out << clusterEnd << "\t";
			out << splitAlignment.mAlignRefName[clusterEnd] << "\t";
			out << splitAlignment.mAlignStrand[clusterEnd] << "\t";
			out << splitAlignment.mSplitAlignSeqStart[clusterEnd] << "\t";
			out << splitAlignment.mSplitAlignSeqLength[clusterEnd] << "\t";
			out << splitAlignment.mSplitSeqStrand[clusterEnd] << "\t";
			out << splitAlignment.mSplitAlignSeq[clusterEnd] << "\t";
			out << splitAlignment.mSplitRemainderSeq[clusterEnd] << "\t";
			out << endl;
		}
	}
}
		
void SplitAlignment::ReadCandidateReads(istream& in, SplitAlignmentMap& splitAlignments)
{
	string line;
	while (getline(in, line))
	{
		vector<string> fields;
		split(fields, line, is_any_of("\t"));
		
		if (fields.size() < 4)
		{
			cerr << "Error: Format error for candidate reads line:" << endl << line << endl;
			exit(1);
		}
		
		int id = lexical_cast<int>(fields[0]);
		int fragmentIndex = lexical_cast<int>(fields[1]);
		int readEnd = lexical_cast<int>(fields[2]);
		int revComp = lexical_cast<int>(fields[3]);
		
		ReadID readID;
		readID.fragmentIndex = fragmentIndex;
		readID.readEnd = readEnd;

		splitAlignments[id].mCandidateReadID.push_back(readID.id);
		splitAlignments[id].mCandidateRevComp.push_back(revComp);
	}
}

void SplitAlignment::ReadCandidateRegions(istream& in, SplitAlignmentMap& splitAlignments)
{
	string line;
	while (getline(in, line))
	{
		vector<string> fields;
		split(fields, line, is_any_of("\t"));
		
		if (fields.size() < 9)
		{
			cerr << "Error: Format error for candidate reads line:" << endl << line << endl;
			exit(1);
		}
		
		int id = lexical_cast<int>(fields[0]);
		int clusterEnd = lexical_cast<int>(fields[1]);
		string alignRefName = fields[2];
		int alignStrand = lexical_cast<int>(fields[3]);
		int alignSeqStart = lexical_cast<int>(fields[4]);
		int alignSeqLength = lexical_cast<int>(fields[5]);
		int seqStrand = lexical_cast<int>(fields[6]);
		string alignSeq = fields[7];
		string remainderSeq = fields[8];
		
		splitAlignments[id].mAlignRefName[clusterEnd] = alignRefName;
		splitAlignments[id].mAlignStrand[clusterEnd] = alignStrand;
		splitAlignments[id].mSplitAlignSeqStart[clusterEnd] = alignSeqStart;
		splitAlignments[id].mSplitAlignSeqLength[clusterEnd] = alignSeqLength;
		splitAlignments[id].mSplitSeqStrand[clusterEnd] = seqStrand;
		splitAlignments[id].mSplitAlignSeq[clusterEnd] = alignSeq;
		splitAlignments[id].mSplitRemainderSeq[clusterEnd] = remainderSeq;
	}
}

void SplitAlignment::ReadCandidateSequences(IReadStream* readStream, SplitAlignmentMap& splitAlignments)
{
	unordered_map<int,IntegerPairVec> candidateReadMap;
	
	for (SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		SplitAlignment& splitAlignment = splitAlignIter->second;
		
		splitAlignment.mCandidateSequence.resize(splitAlignment.mCandidateReadID.size());

		for (int candidateIndex = 0; candidateIndex < splitAlignment.mCandidateReadID.size(); candidateIndex++)
		{
			int readID = splitAlignment.mCandidateReadID[candidateIndex];
			candidateReadMap[readID].push_back(IntegerPair(id,candidateIndex));
		}
	}
	
	RawRead rawRead;
	while (readStream->GetNextRead(rawRead))
	{
		ReadID readID;
		readID.fragmentIndex = lexical_cast<int>(rawRead.fragment);
		readID.readEnd = rawRead.readEnd;
		
		if (candidateReadMap.find(readID.id) == candidateReadMap.end())
		{
			continue;
		}
		
		for (IntegerPairVecConstIter candidateIter = candidateReadMap[readID.id].begin(); candidateIter != candidateReadMap[readID.id].end(); candidateIter++)
		{
			int id = candidateIter->first;
			int candidateIndex = candidateIter->second;
			
			splitAlignments.find(id)->second.mCandidateSequence[candidateIndex] = rawRead.sequence;
		}
	}
}

void SplitAlignment::ReadCandidateSequences(const ReadIndex& readIndex, SplitAlignmentMap& splitAlignments)
{
	unordered_map<int,IntegerPairVec> candidateReadMap;
	
	for (SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		SplitAlignment& splitAlignment = splitAlignIter->second;
		
		splitAlignment.mCandidateSequence.resize(splitAlignment.mCandidateReadID.size());
		
		for (int candidateIndex = 0; candidateIndex < splitAlignment.mCandidateReadID.size(); candidateIndex++)
		{
			ReadID readID;
			readID.id = splitAlignment.mCandidateReadID[candidateIndex];
			
			string readSeq;
			readIndex.Find(readID.fragmentIndex, readID.readEnd, readSeq);
			
			splitAlignment.mCandidateSequence[candidateIndex] = readSeq;
		}
	}
}

bool SplitAlignment::Align(bool generateAlignmentText)
{
	if (mCandidateSequence.size() == 0)
	{
		return false;
	}
	
	SplitReadAligner splitReadAligner(matchScore, mismatchScore, gapScore, false, minAnchor * matchScore, mSplitAlignSeq[0], mSplitAlignSeq[1]);
	
	for (int candidateIndex = 0; candidateIndex < mCandidateReadID.size(); candidateIndex++)
	{
		ReadID readID;
		readID.id = mCandidateReadID[candidateIndex];
		int revComp = mCandidateRevComp[candidateIndex];
		string readSeq = mCandidateSequence[candidateIndex];
		
		if (revComp)
		{
			ReverseComplement(readSeq);
		}
		
		splitReadAligner.Align(readSeq);
		
		SplitReadAlignVec splitAlignments;
		splitReadAligner.GetAlignments(splitAlignments, (int)((float)readSeq.length() * (float)matchScore * 0.90), true, false);
		
		unordered_set<IntegerPair> readSplits;
		
		for (int splitAlignIndex = 0; splitAlignIndex < splitAlignments.size(); splitAlignIndex++)
		{
			const SplitReadAlignment& splitAlignment = splitAlignments[splitAlignIndex];
			
			if (readSplits.find(splitAlignment.split) != readSplits.end())
			{
				continue;
			}
			
			readSplits.insert(splitAlignment.split);
			
			mAlignmentReadID.push_back(readID.id);
			mAlignmentSplit.push_back(splitAlignment.split);
			mAlignmentMatches.push_back(IntegerPair(splitAlignment.matches1.size(),splitAlignment.matches2.size()));
			mAlignmentScore.push_back(min(splitAlignment.score1,splitAlignment.score2));
			
			if (generateAlignmentText)
			{
				stringstream alignmentText;
				
				alignmentText << readID.fragmentIndex << ((readID.readEnd == 0) ? "/1" : "/2") << endl;
				
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
				
				mAlignmentText.push_back(alignmentText.str());
			}
		}
	}
	
	return true;
}

void SplitAlignment::WriteAlignments(ostream& out, SplitAlignmentMap& splitAlignments)
{
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;
		
		for (int alignmentIndex = 0; alignmentIndex < splitAlignment.mAlignmentReadID.size(); alignmentIndex++)
		{
			out << id << "\t";
			out << splitAlignment.mAlignmentReadID[alignmentIndex] << "\t";
			out << splitAlignment.mAlignmentSplit[alignmentIndex].first << "\t" << splitAlignment.mAlignmentSplit[alignmentIndex].second << "\t";
			out << splitAlignment.mAlignmentMatches[alignmentIndex].first << "\t" << splitAlignment.mAlignmentMatches[alignmentIndex].second << "\t";
			out << splitAlignment.mAlignmentScore[alignmentIndex] << "\t";
			out << endl;
		}
	}
}

void SplitAlignment::ReadAlignments(istream& in, SplitAlignmentMap& splitAlignments)
{
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
		
		int id = lexical_cast<int>(fields[0]);
		int readID = lexical_cast<int>(fields[1]);
		int splitFirst = lexical_cast<int>(fields[2]);
		int splitSecond = lexical_cast<int>(fields[3]);
		int matchFirst = lexical_cast<int>(fields[4]);
		int matchSecond = lexical_cast<int>(fields[5]);
		int score = lexical_cast<int>(fields[6]);
		
		splitAlignments[id].mAlignmentReadID.push_back(readID);
		splitAlignments[id].mAlignmentSplit.push_back(IntegerPair(splitFirst,splitSecond));
		splitAlignments[id].mAlignmentMatches.push_back(IntegerPair(matchFirst,matchSecond));
		splitAlignments[id].mAlignmentScore.push_back(score);
	}
}

bool SplitAlignment::Evaluate()
{
	mBreakPos = IntegerVec(2);
	mSequence = "N";
	mSplitReadCount = 0;
	mSplitPosAvg = -1.0;
	mSplitMinAvg = -1.0;
		
	if (mAlignmentReadID.size() == 0)
	{
		return false;
	}
	
	typedef map<IntegerPair,IntegerVec> IntPairIntVecMap;
	typedef map<IntegerPair,Counter> SplitScoreMap;
	typedef map<IntegerPair,Counter> SplitCountMap;
	
	IntPairIntVecMap splitReadIDsMap;
	IntPairIntVecMap leftBaseCounts;
	IntPairIntVecMap rightBaseCounts;
	SplitScoreMap splitScoreMap;
	SplitCountMap splitCountMap;

	for (int alignmentIndex = 0; alignmentIndex < mAlignmentReadID.size(); alignmentIndex++)
	{
		const IntegerPair& split = mAlignmentSplit[alignmentIndex];
		
		splitReadIDsMap[split].push_back(mAlignmentReadID[alignmentIndex]);
		leftBaseCounts[split].push_back(mAlignmentMatches[alignmentIndex].first);
		rightBaseCounts[split].push_back(mAlignmentMatches[alignmentIndex].second);
		splitScoreMap[split] += mAlignmentScore[alignmentIndex];
		splitCountMap[split]++;
	}
	
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
		return false;
	}
	
	mBestSplit = maxSplitScoreIter->first;
	
	DebugCheck(mBestSplit.first <= mSplitAlignSeq[0].length());
	DebugCheck(mBestSplit.second + 1 < mSplitAlignSeq[1].length());
	
	string alignBreak1 = mSplitRemainderSeq[0] + mSplitAlignSeq[0].substr(0, mBestSplit.first);
	string alignBreak2 = mSplitAlignSeq[1].substr(mBestSplit.second + 1) + mSplitRemainderSeq[1];
	
	mSequence = alignBreak1 + "|" + alignBreak2;
	
	if (mSplitSeqStrand[0] == PlusStrand)
	{
		mBreakPos[0] = mSplitAlignSeqStart[0] + mBestSplit.first - 1;
	}
	else
	{
		mBreakPos[0] = mSplitAlignSeqStart[0] + mSplitAlignSeqLength[0] - mBestSplit.first;				
	}
	
	if (mSplitSeqStrand[1] == PlusStrand)
	{
		mBreakPos[1] = mSplitAlignSeqStart[1] + mBestSplit.second + 1;
	}
	else
	{
		mBreakPos[1] = mSplitAlignSeqStart[1] + mSplitAlignSeqLength[1] - mBestSplit.second - 2;
	}
	
	mSplitReadCount = splitCountMap[mBestSplit];
	
	const IntegerVec& leftBaseCount = leftBaseCounts[maxSplitScoreIter->first];
	const IntegerVec& rightBaseCount = rightBaseCounts[maxSplitScoreIter->first];
	
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
	
	mSplitPosAvg = posSum / leftBaseCount.size();
	mSplitMinAvg = minSum / leftBaseCount.size();
		
	return true;
}

void SplitAlignment::WriteSequences(ostream& out, SplitAlignmentMap& splitAlignments)
{
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;
		
		out << id << "\t" << splitAlignment.mSequence << "\t" << "0" << "\t" << splitAlignment.mSplitReadCount << "\t" << splitAlignment.mSplitPosAvg << "\t" << splitAlignment.mSplitMinAvg << endl;
	}
}

void SplitAlignment::WriteBreaks(ostream& out, SplitAlignmentMap& splitAlignments)
{
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;
		
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			out << id << "\t" << clusterEnd << "\t" << splitAlignment.mAlignRefName[clusterEnd] << "\t" << (splitAlignment.mAlignStrand[clusterEnd] == PlusStrand ? "+" : "-") << "\t" << splitAlignment.mBreakPos[clusterEnd] << endl;
		}
	}
}

void SplitAlignment::WriteAlignText(ostream& out, SplitAlignmentMap& splitAlignments)
{
	for (SplitAlignmentMapConstIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		int id = splitAlignIter->first;
		const SplitAlignment& splitAlignment = splitAlignIter->second;
		
		DebugCheck(!splitAlignment.mAlignmentText.empty());
		
		out << id << endl;
		out << splitAlignment.mSplitAlignSeq[0] << "|" << splitAlignment.mSplitAlignSeq[1] << endl;
		for (int alignmentIndex = 0; alignmentIndex < splitAlignment.mAlignmentReadID.size(); alignmentIndex++)
		{
			const IntegerPair& split = splitAlignment.mAlignmentSplit[alignmentIndex];
			
			if (split == splitAlignment.mBestSplit)
			{
				out << splitAlignment.mAlignmentText[alignmentIndex];
			}
		}
	}
}

void SplitAlignment::AddAnchoredReads(const AlignmentIndex& discordant, const AlignmentIndex& anchored, const string& transcript, int strand, int revComp, int start, int end)
{
	CompAlignVec candidateMates;
	
	NameIndex referenceNames;
	discordant.Find(transcript, strand, start, end, referenceNames, candidateMates);
	anchored.Find(transcript, strand, start, end, referenceNames, candidateMates);
	
	for (int candidateIndex = 0; candidateIndex < candidateMates.size(); candidateIndex++)
	{
		ReadID candidateReadID;
		candidateReadID.fragmentIndex = candidateMates[candidateIndex].readID.fragmentIndex;
		candidateReadID.readEnd = (candidateMates[candidateIndex].readID.readEnd == 0) ? 1 : 0;
		
		if (mCandidateUnique.insert(IntegerPair(candidateReadID.id,revComp)).second)
		{
			mCandidateReadID.push_back(candidateReadID.id);
			mCandidateRevComp.push_back(revComp);
		}
	}
}

void SplitAlignment::CalculateBreakRegion(int minReadLength, int maxReadLength, int maxFragmentLength, int alignStart, int alignEnd, int strand, int& breakStart, int& breakLength)
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

void SplitAlignment::CalculateSplitMateRegion(int minReadLength, int maxReadLength, int minFragmentLength, int maxFragmentLength, int breakStart, int breakLength, int strand, int& mateRegionStart, int& mateRegionEnd)
{
	if (strand == PlusStrand)
	{
		int breakEnd = breakStart + breakLength - 1;

		mateRegionStart = breakStart - maxFragmentLength + minReadLength;
		mateRegionEnd = breakEnd - minFragmentLength + maxReadLength;
	}
	else
	{
		int breakEnd = breakStart - breakLength + 1;

		mateRegionStart = breakEnd + minFragmentLength - maxReadLength;
		mateRegionEnd = breakStart + maxFragmentLength - minReadLength;
	}
}


