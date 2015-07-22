/*
 *  evalsplitalign.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "SplitAlignment.h"
#include "Parsers.h"

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

	string splitAlignmentsFilename;

	string seqFilename;
	string breakFilename;
	string predAlignFilename;

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

		TCLAP::ValueArg<string> splitAlignmentsFilenameArg("a","align","Split Alignments Filename",true,"","string",cmd);

		TCLAP::ValueArg<string> seqFilenameArg("q","seq","Sequences Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakFilenameArg("b","break","Break Positions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> predAlignFilenameArg("p","predalign","Prediction Split Alignments Filename",true,"","string",cmd);
		
		cmd.parse(argc,argv);

		referenceFasta = referenceFastaArg.getValue();
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		minReadLength = minReadLengthArg.getValue();
		maxReadLength = maxReadLengthArg.getValue();
		fusionRegionsFilename = fusionRegionsFilenameArg.getValue();

		splitAlignmentsFilename = splitAlignmentsFilenameArg.getValue();

		seqFilename = seqFilenameArg.getValue();
		breakFilename = breakFilenameArg.getValue();
		predAlignFilename = predAlignFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	LocationVecMap fusionRegions;
	ReadAlignRegionPairs(fusionRegionsFilename, fusionRegions);
	
	unordered_map<int,SplitAlignmentTask> alignTasks = CreateTasks(referenceFasta, exonRegionsFilename,
		fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength, fusionRegions);

	ifstream splitAlignmentsFile(splitAlignmentsFilename.c_str());
	ofstream seqFile(seqFilename.c_str());
	ofstream breakFile(breakFilename.c_str());
	ofstream predAlignFile(predAlignFilename.c_str());
	
	CheckFile(splitAlignmentsFile, splitAlignmentsFilename);
	CheckFile(seqFile, seqFilename);
	CheckFile(breakFile, breakFilename);
	CheckFile(predAlignFile, predAlignFilename);

	while (splitAlignmentsFile)
	{
		vector<SplitAlignment> alignments = SplitAlignment::ReadSortedAlignments(splitAlignmentsFile);

		if (alignments.empty())
		{
			break;
		}

		int fusionID = alignments.front().fusionID;

		const SplitAlignmentTask& alignTask = alignTasks[fusionID];

		BreakPrediction prediction = alignTask.Evaluate(alignments);

		prediction.WriteSequence(seqFile);
		prediction.WriteBreak(breakFile);
		prediction.WriteAlignments(predAlignFile);
	}
}

