/*
 *  dosplitalign.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "SplitAlignment.h"
#include "Common.h"
#include "DebugCheck.h"
#include "ReadStream.h"
#include "Parsers.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>

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

	string alignmentsFilename;
	string readSeqs1Filename;
	string readSeqs2Filename;

	string splitAlignmentsFilename;

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

		TCLAP::ValueArg<string> alignmentsFilenameArg("i","improper","Improper Alignments Sam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs1FilenameArg("1","seq1","End 1 Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs2FilenameArg("2","seq2","End 2 Sequences",true,"","string",cmd);
		
		TCLAP::ValueArg<string> splitAlignmentsFilenameArg("a","align","Split Alignments Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		referenceFasta = referenceFastaArg.getValue();
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		minReadLength = minReadLengthArg.getValue();
		maxReadLength = maxReadLengthArg.getValue();
		fusionRegionsFilename = fusionRegionsFilenameArg.getValue();

		alignmentsFilename = alignmentsFilenameArg.getValue();
		readSeqs1Filename = readSeqs1FilenameArg.getValue();
		readSeqs2Filename = readSeqs2FilenameArg.getValue();
		
		splitAlignmentsFilename = splitAlignmentsFilenameArg.getValue();
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

	SplitReadRealigner aligner;

	for (unordered_map<int,SplitAlignmentTask>::const_iterator alignTaskIter = alignTasks.begin();
		alignTaskIter != alignTasks.end(); alignTaskIter++)
	{
		aligner.AddTask(alignTaskIter->second);
	}

	IReadStream* readSeqs1Stream = IReadStream::Create(readSeqs1Filename);
	IReadStream* readSeqs2Stream = IReadStream::Create(readSeqs2Filename);
	
	if (!readSeqs1Stream || !readSeqs2Stream)
	{
		cout << "Error: unable to read sequences" << endl;
		exit(1);
	}
	
	aligner.AddReads(readSeqs1Stream);
	aligner.AddReads(readSeqs2Stream);

	AlignmentStream* alignmentStream = new SamAlignmentStream(alignmentsFilename);

	ofstream splitAlignmentsFile(splitAlignmentsFilename.c_str());
	CheckFile(splitAlignmentsFile, splitAlignmentsFilename);

	aligner.DoAlignment(alignmentStream, splitAlignmentsFile);
}



