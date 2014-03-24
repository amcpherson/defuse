/*
 *  initsplitalign.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
#include "FastaIndex.h"
#include "SplitAlignment.h"
#include "Parsers.h"
#include "PathFinder.h"

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
	string alignFilename;
	string candidateRefSeqsFilename;
	string candidateMateRegionsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Initialize Fusion sequence prediction by split reads");
		TCLAP::ValueArg<string> referenceFastaArg("f","fasta","Reference Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> exonRegionsFilenameArg("e","exons","Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","ufrag","Fragment Length Mean",true,0.0,"float",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","sfrag","Fragment Length Standard Deviation",true,0.0,"float",cmd);
		TCLAP::ValueArg<int> minReadLengthArg("n","minread","Minimum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> maxReadLengthArg("x","maxread","Maximum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> alignFilenameArg("i","input","Input Alignment Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> candidateRefSeqsFilenameArg("r","refseqs","Candidate Reference Sequences Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> candidateMateRegionsFilenameArg("m","materegions","Candidate Mate Regions Filename",true,"","string",cmd);
		cmd.parse(argc,argv);

		referenceFasta = referenceFastaArg.getValue();
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		minReadLength = minReadLengthArg.getValue();
		maxReadLength = maxReadLengthArg.getValue();
		alignFilename = alignFilenameArg.getValue();
		candidateRefSeqsFilename = candidateRefSeqsFilenameArg.getValue();
		candidateMateRegionsFilename = candidateMateRegionsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	LocationVecMap alignRegionPairs;
	ReadAlignRegionPairs(alignFilename, alignRegionPairs);
	
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
	
	ofstream candidateRefSeqsFile(candidateRefSeqsFilename.c_str());
	ofstream candidateMateRegionsFile(candidateMateRegionsFilename.c_str());
	
	CheckFile(candidateRefSeqsFile, candidateRefSeqsFilename);
	CheckFile(candidateMateRegionsFile, candidateMateRegionsFilename);

	for (LocationVecMapConstIter pairIter = alignRegionPairs.begin(); pairIter != alignRegionPairs.end(); pairIter++)
	{
		int id = pairIter->first;
		const LocationVec& alignRegionPair = pairIter->second;
		
		SplitAlignment::SplitAlignmentMap splitAlignments;
		
		splitAlignments[id].Initialize(alignRegionPair, reference, exonRegions, fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength);
		
		SplitAlignment::WriteCandidateRefSeqs(candidateRefSeqsFile, splitAlignments);
		SplitAlignment::WriteCandidateMateRegions(candidateMateRegionsFile, splitAlignments);
	}
	
	candidateRefSeqsFile.close();
	candidateMateRegionsFile.close();
}

