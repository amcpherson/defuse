/*
 *  clustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "AlignmentIndex.h"
#include "Common.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
#include "FastaIndex.h"
#include "FusionSequence.h"
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
#include "bam.h"
#include "sam.h"

using namespace boost;
using namespace std;

void ReadAlignments(const string& bamFilename, StringVec& referenceNames, CompAlignVec& alignments)
{
	cout << "Reading alignments and creating fragment name index" << endl;
	
	samfile_t* inBamFile = samopen(bamFilename.c_str(), "rb", 0);
	
	referenceNames.clear();
	for (int targetIndex = 0; targetIndex < inBamFile->header->n_targets; targetIndex++)
	{
		referenceNames.push_back(string(inBamFile->header->target_name[targetIndex]));
	}

	bam1_t bamEntry;
	memset(&bamEntry, 0, sizeof(bam1_t));
	while (bam_read1(inBamFile->x.bam, &bamEntry))
	{
		// Split qname into id and end
		// Fragment index encoded in fragment name
		string qname = string((char*)bamEntry.data);
		vector<string> qnameFields;
		split(qnameFields, qname, is_any_of("/"));
		
		CompactAlignment compactAlignment;
		compactAlignment.readID.fragmentIndex = lexical_cast<int>(qnameFields[0]);
		compactAlignment.readID.readEnd = (qnameFields[1] == "1") ? 0 : 1;
		compactAlignment.refStrand.referenceIndex = bamEntry.core.tid;
		compactAlignment.refStrand.strand = (bamEntry.core.flag & 0x10) ? MinusStrand : PlusStrand;
		compactAlignment.region.start = bamEntry.core.pos + 1;
		compactAlignment.region.end = bamEntry.core.pos + bamEntry.core.l_qseq;
		
		alignments.push_back(compactAlignment);
		
		memset(&bamEntry, 0, sizeof(bam1_t));
	}
	
	samclose(inBamFile);
}

void IntepretAlignString(const string& alignString, FusionSequence::Location& alignLocation)
{
	string::size_type colonDividerPos = alignString.find_first_of(":");
	if (colonDividerPos == string::npos || colonDividerPos == 0)
	{
		cerr << "Error: Unable to interpret strand for " << alignString << endl;
		exit(1);
	}
	
	char strand = alignString[colonDividerPos - 1];
	
	if (strand == '+')
	{
		alignLocation.strand = PlusStrand;
	}
	else if (strand == '-')
	{
		alignLocation.strand = MinusStrand;
	}
	else
	{
		cerr << "Error: Unable to interpret strand for " << alignString << endl;
		exit(1);
	}
	
	vector<string> alignFields;
	split(alignFields, alignString, is_any_of("+-:"));
	
	if (alignFields.size() != 4)
	{
		cerr << "Error: Unable to interpret alignment string " << alignString << endl;
		exit(1);
	}
	
	alignLocation.refName = alignFields[0];
	alignLocation.start = lexical_cast<int>(alignFields[2]);
	alignLocation.end = lexical_cast<int>(alignFields[3]);
}

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
	
	try
	{
		TCLAP::CmdLine cmd("Fusion sequence prediction");
		TCLAP::ValueArg<string> discordantBamFilenameArg("d","discordant","Discordant Alignments Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> anchoredBamFilenameArg("a","anchored","Anchored Alignments Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("f","fasta","Reference Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readsPrefixArg("r","reads","Reads Filename Prefix",true,"","string",cmd);
		TCLAP::ValueArg<string> exonRegionsFilenameArg("e","exons","Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","ufrag","Fragment Length Mean",true,0.0,"float",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","sfrag","Fragment Length Standard Deviation",true,0.0,"float",cmd);
		TCLAP::ValueArg<int> minReadLengthArg("m","minread","Minimum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> maxReadLengthArg("x","maxread","Maximum Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> align1Arg("1","align1","Alignment Region 1 (format: chromosomestrand:start-end)",true,"","string",cmd);
		TCLAP::ValueArg<string> align2Arg("2","align2","Alignment Region 2 (format: chromosomestrand:start-end)",true,"","string",cmd);
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
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
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
	
	FusionSequence fusionSequence(discordant, anchored, reads, reference, exonRegions, fragmentLengthMean, fragmentLengthStdDev, minReadLength, maxReadLength);

	FusionSequence::LocationVec alignRegions(2);
	IntepretAlignString(align1String, alignRegions[0]);
	IntepretAlignString(align2String, alignRegions[1]);
	
	string breakSequence;
	fusionSequence.Calculate(alignRegions, breakSequence);

	cout << breakSequence << endl;
}

