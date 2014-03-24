/*
 *  estislands.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */


#include "EstCatalog.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <tclap/CmdLine.h>

using namespace std;

int main(int argc, char* argv[])
{
	string estAlignmentsFilename;
	string breakAlignmentsFilename;
	string outputAlignmentsFilename;

	try
	{
		TCLAP::CmdLine cmd("Identify products of cotranscribed regions using est islands");
		TCLAP::ValueArg<string> estAlignmentsFilenameArg("e","ests","EST Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakAlignmentsFilenameArg("b","breaks","Breakpoint Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> outputAlignmentsFilenameArg("o","output","Output Filtered Alignments Filename",true,"","string",cmd);
		cmd.parse(argc,argv);

		estAlignmentsFilename = estAlignmentsFilenameArg.getValue();
		breakAlignmentsFilename = breakAlignmentsFilenameArg.getValue();
		outputAlignmentsFilename = outputAlignmentsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	ifstream estAlignmentsFile(estAlignmentsFilename.c_str());
	if (!estAlignmentsFile.good())
	{
		cerr << "Error: Unable to open est file " << estAlignmentsFilename << endl;
		exit(1);
	}
	
	EstCatalog estCatalog;
	estCatalog.ReadEsts(estAlignmentsFile);
	
	estAlignmentsFile.close();

	ifstream breakAlignmentsFile(breakAlignmentsFilename.c_str());
	if (!breakAlignmentsFile.good())
	{
		cerr << "Error: Unable to open break alignments file " << breakAlignmentsFilename << endl;
		exit(1);
	}

	ofstream outputAlignmentsFile(outputAlignmentsFilename.c_str());
	if (!outputAlignmentsFile.good())
	{
		cerr << "Error: Unable to open output file " << outputAlignmentsFilename << endl;
		exit(1);
	}
	
	estCatalog.FilterContainedInEstIslands(breakAlignmentsFile, outputAlignmentsFile);
	
	breakAlignmentsFile.close();
	outputAlignmentsFile.close();
}

