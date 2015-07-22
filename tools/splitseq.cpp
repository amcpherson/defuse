/*
 *  splitseq.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "SplitAlignment.h"
#include "Parsers.h"
#include "ReadIndex.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


int main(int argc, char* argv[])
{
	string referenceFasta;
	string exonRegionsFilename;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	int minReadLength;
	int maxReadLength;
	string fusionRegionsFilename;

	string readsPrefix;
	string splitAlignmentsFilename;
	int queryFusionID;
	
	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction by split reads");
		TCLAP::ValueArg<string> referenceFastaArg("f","fasta","Reference Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> exonRegionsFilenameArg("e","exons","Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","ufrag","Fragment Length Mean",true,0.0,"float",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","sfrag","Fragment Length Standard Deviation",true,0.0,"float",cmd);
		TCLAP::ValueArg<int> minReadLengthArg("n","minread","Minimum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> maxReadLengthArg("x","maxread","Maximum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> fusionRegionsFilenameArg("r","regions","Fusion Regions Filename",true,"","string",cmd);

		TCLAP::ValueArg<string> readsPrefixArg("p","prefix","Reads Filename Prefix",true,"","string",cmd);
		TCLAP::ValueArg<string> splitAlignmentsFilenameArg("a","align","Split Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> queryFusionIDArg("i","id","Query Fusion ID",false,-1,"integer",cmd);

		cmd.parse(argc,argv);

		referenceFasta = referenceFastaArg.getValue();
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		minReadLength = minReadLengthArg.getValue();
		maxReadLength = maxReadLengthArg.getValue();
		fusionRegionsFilename = fusionRegionsFilenameArg.getValue();

		readsPrefix = readsPrefixArg.getValue();
		splitAlignmentsFilename = splitAlignmentsFilenameArg.getValue();
		queryFusionID = queryFusionIDArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	LocationVecMap fusionRegions;
	ReadAlignRegionPairs(fusionRegionsFilename, fusionRegions);

	if (queryFusionID >= 0)
	{
		LocationVecMap queryFusionRegions;

		if (fusionRegions.find(queryFusionID) == fusionRegions.end())
		{
			cerr << "Error: Unable to find fusion " << queryFusionID << endl;
			exit(1);
		}

		queryFusionRegions[queryFusionID] = fusionRegions[queryFusionID];
		fusionRegions = queryFusionRegions;
	}
	
	unordered_map<int,SplitAlignmentTask> alignTasks = CreateTasks(referenceFasta, exonRegionsFilename,
		fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength, fusionRegions);

	SplitReadRealigner aligner;

	ReadIndex reads;
	reads.Open(readsPrefix);

	ifstream splitAlignmentsFile(splitAlignmentsFilename.c_str());
	CheckFile(splitAlignmentsFile, splitAlignmentsFilename);

	while (splitAlignmentsFile)
	{
		vector<SplitAlignment> alignments = SplitAlignment::ReadSortedAlignments(splitAlignmentsFile);

		if (alignments.empty())
		{
			break;
		}

		int fusionID = alignments.front().fusionID;

		if (alignTasks.find(fusionID) == alignTasks.end())
		{
			continue;
		}

		SplitAlignmentTask& alignTask = alignTasks[fusionID];
		alignTask.GenerateAlignmentText = true;

		vector<SplitAlignment> realignments = alignTask.ReAlign(aligner.GetAligner(), reads, alignments);
		alignTask.WriteAlignText(cout, realignments);
	}
}

