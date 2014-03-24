/*
 *  denovoseq.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "AlignmentIndex.h"
#include "Common.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
#include "FastaIndex.h"
#include "DeNovoFusion.h"
#include "Parsers.h"
#include "PathFinder.h"
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
	string discordantBamFilename;
	string anchoredBamFilename;
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
	string seqFilename;
	string breakFilename;
	string logFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction by denovo assembly");
		TCLAP::ValueArg<string> discordantBamFilenameArg("d","discordant","Discordant Alignments Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> anchoredBamFilenameArg("a","anchored","Anchored Alignments Bam Filename",true,"","string",cmd);
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
		TCLAP::ValueArg<string> seqFilenameArg("q","seq","Sequences Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakFilenameArg("b","break","Break Positions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> logFilenameArg("l","log","Log Filename",true,"","string",cmd);
		cmd.parse(argc,argv);

		discordantBamFilename = discordantBamFilenameArg.getValue();
		anchoredBamFilename = anchoredBamFilenameArg.getValue();
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
		seqFilename = seqFilenameArg.getValue();
		breakFilename = breakFilenameArg.getValue();
		logFilename = logFilenameArg.getValue();
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
		
	AlignmentIndex discordant;
	AlignmentIndex anchored;
	ReadIndex reads;
	FastaIndex reference;
	ExonRegions exonRegions;
		
	discordant.Open(discordantBamFilename);
	anchored.Open(anchoredBamFilename);
	reads.Open(readsPrefix);
	reference.Open(referenceFasta);
	
	ifstream exonRegionsFile(exonRegionsFilename.c_str());
	if (!exonRegionsFile.good() || !exonRegions.Read(exonRegionsFile))
	{
		cerr << "Error: Unable to read exon regions file " << exonRegionsFilename << endl;
		exit(1);
	}
	exonRegionsFile.close();
	
	ofstream seqFile(seqFilename.c_str());
	ofstream breakFile(breakFilename.c_str());
	ofstream logFile(logFilename.c_str());
	
	CheckFile(seqFile, seqFilename);
	CheckFile(breakFile, breakFilename);
	CheckFile(logFile, logFilename);

	for (LocationVecMapConstIter pairIter = alignRegionPairs.begin(); pairIter != alignRegionPairs.end(); pairIter++)
	{
		DeNovoFusion deNovoFusion(discordant, anchored, reads, reference, exonRegions, fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength, logFile);

		int pairID = pairIter->first;
		const LocationVec& alignRegionPair = pairIter->second;

		logFile << "id: " << pairID << endl;

		IntegerVec breakPos;
		string sequence;
		int interLen;
		int minCount;
		deNovoFusion.Calculate(alignRegionPair, breakPos, sequence, interLen, minCount);

		seqFile << pairID << "\t" << sequence << "\t" << interLen << "\t" << minCount << endl;
		breakFile << pairID << "\t" << alignRegionPair[0].refName << "\t" << (alignRegionPair[0].strand == PlusStrand ? "+" : "-") << "\t" << breakPos[0] << endl;
		breakFile << pairID << "\t" << alignRegionPair[1].refName << "\t" << (alignRegionPair[1].strand == PlusStrand ? "+" : "-") << "\t" << breakPos[1] << endl;
	}
}
