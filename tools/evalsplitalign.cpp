/*
 *  evalsplitalign.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "SplitAlignment.h"

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
	string candidateRefSeqsFilename;
	string splitAlignmentsFilename;
	string seqFilename;
	string breakFilename;

	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction by split reads");
		TCLAP::ValueArg<string> candidateRefSeqsFilenameArg("r","refseqs","Candidate Reference Sequences Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> splitAlignmentsFilenameArg("a","align","Split Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> seqFilenameArg("q","seq","Sequences Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakFilenameArg("b","break","Break Positions Filename",true,"","string",cmd);
		cmd.parse(argc,argv);

		candidateRefSeqsFilename = candidateRefSeqsFilenameArg.getValue();
		splitAlignmentsFilename = splitAlignmentsFilenameArg.getValue();
		seqFilename = seqFilenameArg.getValue();
		breakFilename = breakFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ifstream candidateRefSeqsFile(candidateRefSeqsFilename.c_str());
	ifstream splitAlignmentsFile(splitAlignmentsFilename.c_str());
	ofstream seqFile(seqFilename.c_str());
	ofstream breakFile(breakFilename.c_str());
	
	CheckFile(candidateRefSeqsFile, candidateRefSeqsFilename);
	CheckFile(splitAlignmentsFile, splitAlignmentsFilename);
	CheckFile(seqFile, seqFilename);
	CheckFile(breakFile, breakFilename);
	
	SplitAlignment::SplitAlignmentMap splitAlignments;
	
	SplitAlignment::ReadAlignments(splitAlignmentsFile, splitAlignments);
	SplitAlignment::ReadCandidateRefSeqs(candidateRefSeqsFile, splitAlignments);

	for (SplitAlignment::SplitAlignmentMapIter splitAlignIter = splitAlignments.begin(); splitAlignIter != splitAlignments.end(); splitAlignIter++)
	{
		SplitAlignment& splitAlignment = splitAlignIter->second;		

		splitAlignment.Evaluate();
	}
	
	SplitAlignment::WriteSequences(seqFile, splitAlignments);
	SplitAlignment::WriteBreaks(breakFile, splitAlignments);	
	
	candidateRefSeqsFile.close();
	splitAlignmentsFile.close();
	seqFile.close();
	breakFile.close();
}

