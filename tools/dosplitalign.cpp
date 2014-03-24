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
	string candidateReadsFilename;
	string candidateRegionsFilename;
	string splitAlignmentsFilename;
	string readSeqs1Filename;
	string readSeqs2Filename;

	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction by split reads");
		TCLAP::ValueArg<string> candidateReadsFilenameArg("r","reads","Candidate Reads Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> candidateRegionsFilenameArg("c","cseqs","Candidate Reference Sequences Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> splitAlignmentsFilenameArg("a","align","Split Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs1FilenameArg("1","seq1","End 1 Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs2FilenameArg("2","seq2","End 2 Sequences",true,"","string",cmd);
		cmd.parse(argc,argv);

		candidateReadsFilename = candidateReadsFilenameArg.getValue();
		candidateRegionsFilename = candidateRegionsFilenameArg.getValue();
		splitAlignmentsFilename = splitAlignmentsFilenameArg.getValue();
		readSeqs1Filename = readSeqs1FilenameArg.getValue();
		readSeqs2Filename = readSeqs2FilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ifstream candidateReadsFile(candidateReadsFilename.c_str());
	ifstream candidateRegionsFile(candidateRegionsFilename.c_str());
	ofstream splitAlignmentsFile(splitAlignmentsFilename.c_str());
	
	CheckFile(candidateReadsFile, candidateReadsFilename);
	CheckFile(candidateRegionsFile, candidateRegionsFilename);
	CheckFile(splitAlignmentsFile, splitAlignmentsFilename);
	
	SplitAlignment::SplitAlignmentMap splitAlignments;
	
	SplitAlignment::ReadCandidateReads(candidateReadsFile, splitAlignments);
	SplitAlignment::ReadCandidateRegions(candidateRegionsFile, splitAlignments);
	
	IReadStream* readSeqs1Stream = IReadStream::Create(readSeqs1Filename);
	IReadStream* readSeqs2Stream = IReadStream::Create(readSeqs2Filename);
	
	if (!readSeqs1Stream || !readSeqs2Stream)
	{
		cout << "Error: unable to read sequences" << endl;
		exit(1);
	}
	
	SplitAlignment::ReadCandidateSequences(readSeqs1Stream, splitAlignments);
	SplitAlignment::ReadCandidateSequences(readSeqs2Stream, splitAlignments);
	
	for (SplitAlignment::SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		SplitAlignment& splitAlignment = splitAlignIter->second;		
		splitAlignment.Align(false);
	}
	
	SplitAlignment::WriteAlignments(splitAlignmentsFile, splitAlignments);

	candidateReadsFile.close();
	candidateRegionsFile.close();
	splitAlignmentsFile.close();
}



