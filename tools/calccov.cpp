/*
 *  splitseq.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "AlignmentIndex.h"
#include "Common.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
#include "FastaIndex.h"
#include "SplitFusion.h"
#include "Parsers.h"
#include "PathFinder.h"
#include "ReadIndex.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/unordered_set.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;

double CalculateLength(CompAlignVec& pair);
double CalculateSplitPos(const Region& alignment, int position, int splitMinAnchor);
double CalculateSplitMin(const Region& alignment, int position, int splitMinAnchor);
void FindSpanning(const AlignmentIndex& concordant, const string& refName, int position, int maxFragmentLength, CompAlignTable& spanning);
void FindSplit(const AlignmentIndex& concordant, const string& refName, int position, int splitMinAnchor, CompAlignVec& split);

int main(int argc, char* argv[])
{
	string concordantBamFilename;
	string covarianceStatsFilename;
	string geneTranscriptsFilename;
	double covSampleDensity;
	int splitMinAnchor;
	int maxFragmentLength;
	
	try
	{
		TCLAP::CmdLine cmd("Calculate covariance stats from concordant alignments");
		TCLAP::ValueArg<string> concordantBamFilenameArg("c","conc","Concordant Alignments Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> covarianceStatsFilenameArg("v","vstat","Covariance Stats Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> geneTranscriptsFilenameArg("g","genetran","Gene Transcripts Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> covSampleDensityArg("d","density","Covariance Sampling Density",true,0.0,"float",cmd);
		TCLAP::ValueArg<int> splitMinAnchorArg("a","anchor","Gene Transcripts Filename",true,0,"integer",cmd);
		TCLAP::ValueArg<int> maxFragmentLengthArg("m","maxfrag","Maximum Fragment Length",true,0,"integer",cmd);
		cmd.parse(argc,argv);

		concordantBamFilename = concordantBamFilenameArg.getValue();
		covarianceStatsFilename = covarianceStatsFilenameArg.getValue();
		geneTranscriptsFilename = geneTranscriptsFilenameArg.getValue();
		covSampleDensity = covSampleDensityArg.getValue();
		splitMinAnchor = splitMinAnchorArg.getValue();
		maxFragmentLength = maxFragmentLengthArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	AlignmentIndex concordant;
	concordant.Open(concordantBamFilename);
	
	ExonRegions geneTranscripts;
	ifstream geneTranscriptsFile(geneTranscriptsFilename.c_str());
	if (!geneTranscriptsFile.good() || !geneTranscripts.ReadGeneTranscripts(geneTranscriptsFile))
	{
		cerr << "Error: Unable to gene transcripts file " << geneTranscriptsFilename << endl;
		exit(1);
	}
	geneTranscriptsFile.close();
	
	// Initialize RNG for predictability
	srand(1);
	
	// Covariance statistics
	double spanLenSum1 = 0.0;
	double spanLenSum2 = 0.0;
	double spanLenCrossSum = 0.0;
	double spanLenCount = 0.0;
	
	double splitPosSum1 = 0.0;
	double splitPosSum2 = 0.0;
	double splitPosCrossSum = 0.0;
	double splitPosCount = 0.0;
	
	double splitMinSum1 = 0.0;
	double splitMinSum2 = 0.0;
	double splitMinCrossSum = 0.0;
	double splitMinCount = 0.0;
	
	// Iterate through a list of acceptable unspliced transcripts
	const StringVec& geneNames = geneTranscripts.GetGeneNames();
	for (StringVecConstIter geneIter = geneNames.begin(); geneIter != geneNames.end(); geneIter++)
	{
		if (geneTranscripts.GetGeneTranscriptNames(*geneIter).size() == 1)
		{
			const string& transcript = geneTranscripts.GetGeneTranscriptNames(*geneIter)[0];

			int length = geneTranscripts.GetTranscriptLength(transcript);
			int numMarkers = length * covSampleDensity;
			for (int positionIndex = 0; positionIndex < numMarkers; positionIndex++)
			{
				int position = rand() % length + 1;
				
				CompAlignTable spanning;
				FindSpanning(concordant, transcript, position, maxFragmentLength, spanning);
				
				// Require at least 2 spanning
				if (spanning.size() >= 2)
				{
					// Randomly select two spanning fragments and calculate lengths
					double fragmentLength1 = CalculateLength(spanning[rand() % spanning.size()]);
					double fragmentLength2 = CalculateLength(spanning[rand() % spanning.size()]);
					
					spanLenSum1 += fragmentLength1;
					spanLenSum2 += fragmentLength2;
					spanLenCrossSum += fragmentLength1 * fragmentLength2;
					spanLenCount++;
				}
								
				CompAlignVec split;
				FindSplit(concordant, transcript, position, splitMinAnchor, split);
				
				// Require at least 2 split
				if (split.size() >= 2)
				{
					// Randomly select the alignment regions of two split reads
					Region splitAlignment1 = split[rand() % split.size()].region;
					Region splitAlignment2 = split[rand() % split.size()].region;
					
					double splitPos1 = CalculateSplitPos(splitAlignment1, position, splitMinAnchor);
					double splitPos2 = CalculateSplitPos(splitAlignment2, position, splitMinAnchor);
					
					double splitMin1 = CalculateSplitMin(splitAlignment1, position, splitMinAnchor);
					double splitMin2 = CalculateSplitMin(splitAlignment2, position, splitMinAnchor);
					
					splitPosSum1 += splitPos1;
					splitPosSum2 += splitPos2;
					splitPosCrossSum += splitPos1 * splitPos2; 
					splitPosCount++;
					
					splitMinSum1 += splitMin1;
					splitMinSum2 += splitMin2;
					splitMinCrossSum += splitMin1 * splitMin2; 
					splitMinCount++;
				}
			}
		}
	}
	
	double spanLenMean = (spanLenSum1 + spanLenSum2) / (2.0 * spanLenCount);
	double spanLenCov = (spanLenCrossSum - spanLenSum1 * spanLenSum2 / spanLenCount) / spanLenCount;
	
	double splitPosMean = (splitPosSum1 + splitPosSum2) / (2.0 * splitPosCount);
	double splitPosCov = (splitPosCrossSum - splitPosSum1 * splitPosSum2 / splitPosCount) / splitPosCount;
	
	double splitMinMean = (splitMinSum1 + splitMinSum2) / (2.0 * splitMinCount);
	double splitMinCov = (splitMinCrossSum - splitMinSum1 * splitMinSum2 / splitMinCount) / splitMinCount;
	
	spanLenCov = max(0.0, spanLenCov);
	splitPosCov = max(0.0, splitPosCov);
	splitMinCov = max(0.0, splitMinCov);
	
	ofstream covarianceStatsFile(covarianceStatsFilename.c_str());
	CheckFile(covarianceStatsFile, covarianceStatsFilename);
	
	covarianceStatsFile << "span_mean\tspan_covariance\tsplit_pos_mean\tsplit_pos_covariance\tsplit_min_mean\tsplit_min_covariance\n";
	covarianceStatsFile << spanLenMean << "\t" << spanLenCov << "\t" << splitPosMean << "\t" << splitPosCov << "\t" << splitMinMean << "\t" << splitMinCov << endl;
}

double CalculateLength(CompAlignVec& pair)
{
	if (pair.size() != 2)
	{
		cerr << "Error: Invalid call to CalculateLength" << endl;
		exit(1);
	}

	return pair[1].region.end - pair[0].region.start + 1;
}

double CalculateSplitPos(const Region& alignment, int position, int splitMinAnchor)
{
	double splitPosValue = max(0.0, (double)(position - alignment.start - splitMinAnchor));
	double splitPosRange = alignment.end - alignment.start + 1.0 - 2.0 * splitMinAnchor;
	
	return splitPosValue / splitPosRange;
}

double CalculateSplitMin(const Region& alignment, int position, int splitMinAnchor)
{
	double splitMinValue = max(0.0, (double)min(position - alignment.start - splitMinAnchor, alignment.end - position - splitMinAnchor));
	double splitMinRange = 0.5 * (alignment.end - alignment.start + 1.0 - 2.0 * splitMinAnchor);
	
	return splitMinValue / splitMinRange;
}

void FindSpanning(const AlignmentIndex& concordant, const string& refName, int position, int maxFragmentLength, CompAlignTable& spanning)
{
	NameIndex referenceNames;

	unordered_set<int> fragments;
	unordered_map<int,CompactAlignment> paired[2];

	CompAlignVec upstream;
	concordant.Find(refName, PlusStrand, position - maxFragmentLength, position, referenceNames, upstream);
	
	for (CompAlignVecConstIter alignIter = upstream.begin(); alignIter != upstream.end(); alignIter++)
	{
		if (alignIter->region.end < position)
		{
			paired[0][alignIter->readID.fragmentIndex] = *alignIter;
			fragments.insert(alignIter->readID.fragmentIndex);
		}
	}
	
	CompAlignVec downstream;
	concordant.Find(refName, MinusStrand, position, position + maxFragmentLength, referenceNames, downstream);
	
	for (CompAlignVecConstIter alignIter = downstream.begin(); alignIter != downstream.end(); alignIter++)
	{
		if (alignIter->region.start > position)
		{
			paired[1][alignIter->readID.fragmentIndex] = *alignIter;
			fragments.insert(alignIter->readID.fragmentIndex);
		}
	}
	
	for (unordered_set<int>::const_iterator fragmentIter = fragments.begin(); fragmentIter != fragments.end(); fragmentIter++)
	{
		int fragmentIndex = *fragmentIter;
		
		if (paired[0].find(fragmentIndex) != paired[0].end() && paired[1].find(fragmentIndex) != paired[1].end() && 
		    paired[0].find(fragmentIndex)->second.readID.readEnd != paired[1].find(fragmentIndex)->second.readID.readEnd)
		{
			spanning.push_back(CompAlignVec());
			spanning.back().push_back(paired[0].find(fragmentIndex)->second);
			spanning.back().push_back(paired[1].find(fragmentIndex)->second);
		}
	}
}

void FindSplit(const AlignmentIndex& concordant, const string& refName, int position, int splitMinAnchor, CompAlignVec& split)
{
	NameIndex referenceNames;
	
	CompAlignVec overlapping;
	concordant.Find(refName, PlusStrand, position, position, referenceNames, overlapping);
	concordant.Find(refName, MinusStrand, position, position, referenceNames, overlapping);
	
	for (CompAlignVecConstIter alignIter = overlapping.begin(); alignIter != overlapping.end(); alignIter++)
	{
		if (position - alignIter->region.start >= splitMinAnchor && alignIter->region.end - position >= splitMinAnchor)
		{
			split.push_back(*alignIter);
		}
	}
}

