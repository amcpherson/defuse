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
	string alignmentsFilename;
	string candidateRefSeqsFilename;
	string candidateMateRegionsFilename;
	string splitAlignmentsFilename;
	string readSeqs1Filename;
	string readSeqs2Filename;

	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction by split reads");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Alignments Sam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> candidateRefSeqsFilenameArg("r","refseqs","Candidate Reference Sequences Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> candidateMateRegionsFilenameArg("m","materegions","Candidate Mate Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> splitAlignmentsFilenameArg("s","splits","Split Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs1FilenameArg("1","seq1","End 1 Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs2FilenameArg("2","seq2","End 2 Sequences",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		alignmentsFilename = alignmentsFilenameArg.getValue();
		candidateRefSeqsFilename = candidateRefSeqsFilenameArg.getValue();
		candidateMateRegionsFilename = candidateMateRegionsFilenameArg.getValue();
		splitAlignmentsFilename = splitAlignmentsFilenameArg.getValue();
		readSeqs1Filename = readSeqs1FilenameArg.getValue();
		readSeqs2Filename = readSeqs2FilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ifstream candidateRefSeqsFile(candidateRefSeqsFilename.c_str());
	ifstream candidateMateRegionsFile(candidateMateRegionsFilename.c_str());
	ofstream splitAlignmentsFile(splitAlignmentsFilename.c_str());
	
	CheckFile(candidateRefSeqsFile, candidateRefSeqsFilename);
	CheckFile(candidateMateRegionsFile, candidateMateRegionsFilename);
	CheckFile(splitAlignmentsFile, splitAlignmentsFilename);
	
	SplitAlignment::SplitAlignmentMap splitAlignments;
	
	SplitAlignment::ReadCandidateRefSeqs(candidateRefSeqsFile, splitAlignments);
	SplitAlignment::ReadCandidateMateRegions(candidateMateRegionsFile, splitAlignments);
	
	AlignmentStream* alignmentStream = new SamAlignmentStream(alignmentsFilename);
	
	SplitAlignment::FindCandidates(alignmentStream, splitAlignments);
	
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

	candidateRefSeqsFile.close();
	candidateMateRegionsFile.close();
	splitAlignmentsFile.close();
}



