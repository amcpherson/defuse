/*
 *  splitseq.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "AlignmentStream.h"
#include "Common.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
#include "Indexer.h"

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

double CalculateSplitPos(int start, int end, int position, int splitMinAnchor);
double CalculateSplitMin(int start, int end, int position, int splitMinAnchor);

class BinnedSamplePositions
{
public:
	BinnedSamplePositions(int binSpacing) : mBinSpacing(binSpacing) {}
	
	void Load( IntegerPairVec& samplePositions)
	{
		for (int samplePosIndex = 0; samplePosIndex < samplePositions.size(); samplePosIndex++)
		{
			int bin = samplePositions[samplePosIndex].second / mBinSpacing;
			mBinned[IntegerPair(samplePositions[samplePosIndex].first,bin)].push_back(samplePosIndex);
		}
	}
	
	void ApproxContained(int ref, int start, int end, unordered_set<int>& indices) const
	{
		int startBin = start / mBinSpacing;
		int endBin = end / mBinSpacing;
		
		for (int bin = startBin; bin <= endBin; bin++)
		{
			unordered_map<IntegerPair,IntegerVec>::const_iterator findIter = mBinned.find(IntegerPair(ref,bin));
			if (findIter != mBinned.end())
			{
				indices.insert(findIter->second.begin(),findIter->second.end());
			}
		}
	}
	
private:	
	int mBinSpacing;
	unordered_map<IntegerPair,IntegerVec> mBinned;
};

int main(int argc, char* argv[])
{
	string concordantSamFilename;
	string geneTranscriptsFilename;
	string lengthSamplesFilename;
	string splitPosSamplesFilename;
	string splitMinSamplesFilename;
	double covSampleDensity;
	int splitMinAnchor;
	int trimLength;
	bool multiExonTransripts;

	try
	{
		TCLAP::CmdLine cmd("Calculate covariance stats from concordant alignments");
		TCLAP::ValueArg<string> concordantSamFilenameArg("c","conc","Concordant Sam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> geneTranscriptsFilenameArg("g","genetran","Gene Transcripts Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> lengthSamplesFilenameArg("l","len","Spanning Length Samples Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> splitPosSamplesFilenameArg("p","pos","Split Position Samples Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> splitMinSamplesFilenameArg("m","min","Split Minimum Samples Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> covSampleDensityArg("d","density","Covariance Sampling Density",true,0.0,"float",cmd);
		TCLAP::ValueArg<int> splitMinAnchorArg("a","anchor","Gene Transcripts Filename",true,0,"integer",cmd);
		TCLAP::ValueArg<int> trimLengthArg("t","trim","Trim Length for Spanning Alignments",true,0,"integer",cmd);
		TCLAP::SwitchArg multiExonTransriptsArg("","multiexon","Use Multi-Exon Transcripts",cmd);
		cmd.parse(argc,argv);
		
		concordantSamFilename = concordantSamFilenameArg.getValue();
		geneTranscriptsFilename = geneTranscriptsFilenameArg.getValue();
		lengthSamplesFilename = lengthSamplesFilenameArg.getValue();
		splitPosSamplesFilename = splitPosSamplesFilenameArg.getValue();
		splitMinSamplesFilename = splitMinSamplesFilenameArg.getValue();
		covSampleDensity = covSampleDensityArg.getValue();
		splitMinAnchor = splitMinAnchorArg.getValue();
		trimLength = trimLengthArg.getValue();
		multiExonTransripts = multiExonTransriptsArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ExonRegions geneTranscripts;
	ifstream geneTranscriptsFile(geneTranscriptsFilename.c_str());
	if (!geneTranscriptsFile.good() || !geneTranscripts.Read(geneTranscriptsFile))
	{
		cerr << "Error: Unable to gene transcripts file " << geneTranscriptsFilename << endl;
		exit(1);
	}
	geneTranscriptsFile.close();
	
	// Initialize RNG for predictability
	srand(11);

	IntegerPairVec samplePositions;
	NameIndex refNames;
	
	// Iterate through a list of acceptable unspliced transcripts
	const StringVec& genes = geneTranscripts.GetGenes();
	for (StringVecConstIter geneIter = genes.begin(); geneIter != genes.end(); geneIter++)
	{
		if (geneTranscripts.GetGeneTranscripts(*geneIter).size() == 1 || multiExonTransripts)
		{
			const string& transcript = geneTranscripts.GetGeneTranscripts(*geneIter)[0];
			
			string transcriptID = *geneIter + "|" + transcript;
			
			int refNameIndex = refNames.Index(transcriptID);
			
			int length = geneTranscripts.GetTranscriptLength(transcript);
			int numMarkers = (int)(length * covSampleDensity);
			for (int positionIndex = 0; positionIndex < numMarkers; positionIndex++)
			{
				int position = rand() % length + 1;
				samplePositions.push_back(IntegerPair(refNameIndex,position));
			}
		}
	}
	
	BinnedSamplePositions binnedSamplePositions(10000);
	binnedSamplePositions.Load(samplePositions);
	
	vector<pair<int,int> > spanLengths;
	vector<pair<int,double> > splitPositions;
	vector<pair<int,double> > splitMinimums;
	
	AlignmentStream* alignmentStream = new SamAlignmentStream(concordantSamFilename);
	FragmentAlignmentStream fragmentAlignments(alignmentStream);
	
	RawAlignmentVec alignments;
	while (fragmentAlignments.GetNextAlignments(alignments))
	{
		DebugCheck(alignments.size() > 0);
		
		if (alignments.size() != 2)
		{
			cerr << "Error: expected 2 alignments per fragment" << endl;
			cerr << "retrieved " << alignments.size() << " alignments for " << alignments[0].fragment << endl;		
			exit(1);
		}
		
		if (!refNames.Contains(alignments[0].reference))
		{
			continue;
		}
		
		int refNameIndex = refNames.Index(alignments[0].reference);
		
		int unseqStart = min(alignments[0].region.start + trimLength, alignments[1].region.start + trimLength);
		int unseqEnd = max(alignments[0].region.end - trimLength, alignments[1].region.end - trimLength);
		
		int fragmentLength = max(alignments[0].region.end, alignments[1].region.end) - min(alignments[0].region.start, alignments[1].region.start);
		
		unordered_set<int> unseqContained;
		binnedSamplePositions.ApproxContained(refNameIndex, unseqStart, unseqEnd, unseqContained);
		
		for (unordered_set<int>::const_iterator sampleIter = unseqContained.begin(); sampleIter != unseqContained.end(); sampleIter++)
		{
			int samplePos = samplePositions[*sampleIter].second;
			if (samplePos >= unseqStart && samplePos <= unseqEnd)
			{
				spanLengths.push_back(pair<int,int>(*sampleIter,fragmentLength));
			}
		}
		
		for (int readEnd = 0; readEnd <= 1; readEnd++)
		{
			int anchoredStart = alignments[readEnd].region.start + splitMinAnchor;
			int anchoredEnd = alignments[readEnd].region.end - splitMinAnchor + 1;
			
			unordered_set<int> readContained;
			binnedSamplePositions.ApproxContained(refNameIndex, anchoredStart, anchoredEnd, readContained);
			
			for (unordered_set<int>::const_iterator sampleIter = readContained.begin(); sampleIter != readContained.end(); sampleIter++)
			{
				int samplePos = samplePositions[*sampleIter].second;
				if (samplePos >= anchoredStart && samplePos <= anchoredEnd)
				{
					splitPositions.push_back(pair<int,double>(*sampleIter,CalculateSplitPos(alignments[readEnd].region.start, alignments[readEnd].region.end, samplePos, splitMinAnchor)));
					splitMinimums.push_back(pair<int,double>(*sampleIter,CalculateSplitMin(alignments[readEnd].region.start, alignments[readEnd].region.end, samplePos, splitMinAnchor)));
				}
			}
		}
	}
	
	ofstream lengthSamplesFile(lengthSamplesFilename.c_str());
	CheckFile(lengthSamplesFile, lengthSamplesFilename);
	for (vector<pair<int,int> >::const_iterator sampleIter = spanLengths.begin(); sampleIter != spanLengths.end(); sampleIter++)
	{
		lengthSamplesFile << sampleIter->first << "\t" << sampleIter->second << endl;
	}
	lengthSamplesFile.close();
	
	ofstream splitPosSamplesFile(splitPosSamplesFilename.c_str());
	CheckFile(splitPosSamplesFile, splitPosSamplesFilename);
	for (vector<pair<int,double> >::const_iterator sampleIter = splitPositions.begin(); sampleIter != splitPositions.end(); sampleIter++)
	{
		splitPosSamplesFile << sampleIter->first << "\t" << sampleIter->second << endl;
	}
	splitPosSamplesFile.close();
	
	ofstream splitMinSamplesFile(splitMinSamplesFilename.c_str());
	CheckFile(splitMinSamplesFile, splitMinSamplesFilename);
	for (vector<pair<int,double> >::const_iterator sampleIter = splitMinimums.begin(); sampleIter != splitMinimums.end(); sampleIter++)
	{
		splitMinSamplesFile << sampleIter->first << "\t" << sampleIter->second << endl;
	}
	splitMinSamplesFile.close();
}

double CalculateSplitPos(int start, int end, int position, int splitMinAnchor)
{
	double splitPosValue = max(0.0, (double)(position - start - splitMinAnchor));
	double splitPosRange = end - start + 1.0 - 2.0 * splitMinAnchor;
	
	return splitPosValue / splitPosRange;
}

double CalculateSplitMin(int start, int end, int position, int splitMinAnchor)
{
	double splitMinValue = max(0.0, (double)min(position - start - splitMinAnchor, end + 1 - position - splitMinAnchor));
	double splitMinRange = floor(0.5 * (end - start + 1.0 - 2.0 * splitMinAnchor));
	
	return splitMinValue / splitMinRange;
}

