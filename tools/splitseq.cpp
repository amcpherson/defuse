/*
 *  splitseq.cpp
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
	string alignmentsFilename;
	string readsPrefix;
	string referenceFasta;
	string exonRegionsFilename;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	int minReadLength;
	int maxReadLength;
	string align1String;
	string align2String;
	string alignFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction by split reads");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Alignments Sam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("f","fasta","Reference Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readsPrefixArg("r","reads","Reads Filename Prefix",true,"","string",cmd);
		TCLAP::ValueArg<string> exonRegionsFilenameArg("e","exons","Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","ufrag","Fragment Length Mean",true,0.0,"float",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","sfrag","Fragment Length Standard Deviation",true,0.0,"float",cmd);
		TCLAP::ValueArg<int> minReadLengthArg("m","minread","Minimum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> maxReadLengthArg("x","maxread","Maximum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> align1Arg("1","align1","Alignment Region 1 (format: chromosomestrand:start-end)",false,"","string",cmd);
		TCLAP::ValueArg<string> align2Arg("2","align2","Alignment Region 2 (format: chromosomestrand:start-end)",false,"","string",cmd);
		TCLAP::ValueArg<string> alignFilenameArg("i","input","Input Alignment Regions Filename",false,"","string",cmd);
		cmd.parse(argc,argv);

		alignmentsFilename = alignmentsFilenameArg.getValue();
		readsPrefix = readsPrefixArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		minReadLength = minReadLengthArg.getValue();
		maxReadLength = maxReadLengthArg.getValue();
		align1String = align1Arg.getValue();
		align2String = align2Arg.getValue();
		alignFilename = alignFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "Error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	LocationVecMap alignRegionPairs;
	if (!align1String.empty() && !align2String.empty())
	{
		alignRegionPairs[0] = LocationVec(2);
		
		IntepretAlignString(align1String, alignRegionPairs[0][0]);
		IntepretAlignString(align2String, alignRegionPairs[0][1]);
	}
	else if (!alignFilename.empty())
	{
		ReadAlignRegionPairs(alignFilename, alignRegionPairs);
	}
	else
	{
		cerr << "Error: You must specify either a filename or 2 alignment strings" << endl;
		exit(1);
	}
	
	ReadIndex reads;
	FastaIndex reference;
	ExonRegions exonRegions;
		
	reads.Open(readsPrefix);
	reference.Open(referenceFasta);
	
	ifstream exonRegionsFile(exonRegionsFilename.c_str());
	if (!exonRegionsFile.good() || !exonRegions.Read(exonRegionsFile))
	{
		cerr << "Error: Unable to read exon regions file " << exonRegionsFilename << endl;
		exit(1);
	}
	exonRegionsFile.close();
	
	SplitAlignment::SplitAlignmentMap splitAlignments;
	for (LocationVecMapConstIter pairIter = alignRegionPairs.begin(); pairIter != alignRegionPairs.end(); pairIter++)
	{
		int id = pairIter->first;
		const LocationVec& alignRegionPair = pairIter->second;
		
		splitAlignments[id].Initialize(alignRegionPair, reference, exonRegions, fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength);
	}
	
	AlignmentStream* alignmentStream = new SamAlignmentStream(alignmentsFilename);
	
	SplitAlignment::FindCandidates(alignmentStream, splitAlignments);
	SplitAlignment::ReadCandidateSequences(reads, splitAlignments);
	
	for (LocationVecMapConstIter pairIter = alignRegionPairs.begin(); pairIter != alignRegionPairs.end(); pairIter++)
	{
		int id = pairIter->first;
		
		splitAlignments[id].Align(true);
		splitAlignments[id].Evaluate();
	}

	SplitAlignment::WriteAlignText(cout, splitAlignments);
}

