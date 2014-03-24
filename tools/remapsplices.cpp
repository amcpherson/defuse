/*
 *  remapsplices.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "Alignments.h"
#include "CompactBreakRegion.h"
#include "DebugCheck.h"
#include "ExonRegions.h"

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
	string exonRegionsFilename;
	string geneTranscriptsFilename;
	string breakRegionsFilename;
	string outputBreakRegionsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Remapping of splice variant break regions tool");
		TCLAP::ValueArg<string> exonRegionsFilenameArg("e","exons","Input Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> geneTranscriptsFilenameArg("g","genetran","Gene Transcripts Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakRegionsFilenameArg("b","breaks","Input Breakpoint Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> outputBreakRegionsFilenameArg("o","output","Output Breakpoint Regions Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		exonRegionsFilename = exonRegionsFilenameArg.getValue();
		geneTranscriptsFilename = geneTranscriptsFilenameArg.getValue();
		breakRegionsFilename = breakRegionsFilenameArg.getValue();
		outputBreakRegionsFilename = outputBreakRegionsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	if (sizeof(long) != 8)
	{
		cerr << "Error: sizeof(long) != 8" << endl;		
		exit(1);
	}
	
	cout << "Reading exon regions and gene transcripts" << endl;

	ExonRegions exonRegions;

	ifstream exonRegionsFile(exonRegionsFilename.c_str());
	if (!exonRegionsFile)
	{
		cout << "Error: unable to read exon regions" << endl;		
		exit(1);
	}
	
	exonRegions.ReadExonRegions(exonRegionsFile);
	
	ifstream geneTranscriptsFile(geneTranscriptsFilename.c_str());
	if (!geneTranscriptsFile)
	{
		cout << "Error: unable to gene transcripts" << endl;		
		exit(1);
	}
	
	exonRegions.ReadGeneTranscripts(geneTranscriptsFile);
	
	cout << "Reading break regions" << endl;
	
	BrRegVec breakRegions;
	NameIndex referenceNames;
	
	ReadBreakRegions(breakRegionsFilename, referenceNames, breakRegions);

	IntegerVec isTranscript(referenceNames.Size(), 0);
	for (int referenceIndex = 0; referenceIndex < referenceNames.Size(); referenceIndex++)
	{
		if (exonRegions.IsTranscript(referenceNames.Get(referenceIndex)))
		{
			isTranscript[referenceIndex] = 1;
		}
	}
	
	ofstream outputBreakRegionsFile(outputBreakRegionsFilename.c_str());
	if (!outputBreakRegionsFile)
	{
		cout << "Error: unable to write to break regions" << endl;		
		exit(1);
	}
	
	for (BrRegVecIter brRegIter = breakRegions.begin(); brRegIter != breakRegions.end(); brRegIter++)
	{
		if (isTranscript[brRegIter->refStrand.referenceIndex])
		{
			const string& transcript = referenceNames.Get(brRegIter->refStrand.referenceIndex);
			string gene;
			
			bool remapStartResult = exonRegions.RemapTranscriptToGene(transcript, brRegIter->start, gene, brRegIter->start);
			bool remapEndResult = exonRegions.RemapTranscriptToGene(transcript, brRegIter->end, gene, brRegIter->end);

			DebugCheck(remapStartResult && remapEndResult);			

			brRegIter->refStrand.referenceIndex = referenceNames.Index(gene);			
		}

		outputBreakRegionsFile << brRegIter->clusterID << "\t" << referenceNames.Get(brRegIter->refStrand.referenceIndex) << "\t" << ((brRegIter->refStrand.strand == PlusStrand) ? "+" : "-") << "\t" << brRegIter->start << "\t" << brRegIter->end << endl;
	}
	
	outputBreakRegionsFile.close();
}

