/*
 *  connectivity.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-04.
 *
 */

#include "Common.h"
#include "QuaterpartiteEdges.h"
#include "DebugCheck.h"
#include "Parsers.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <tclap/CmdLine.h>

using namespace std;

int main(int argc, char* argv[])
{
	string dnaClustersFilename;
	string rnaClustersFilename;
	string overlapsFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Find connected components of the overlapping dna/rna clusters graph");
		TCLAP::ValueArg<string> dnaClustersFilenameArg("d","dna","DNA clusters filename",true,"","string",cmd);
		TCLAP::ValueArg<string> rnaClustersFilenameArg("r","rna","RNA clusters filename",true,"","string",cmd);
		TCLAP::ValueArg<string> overlapsFilenameArg("o","overlaps","DNA-RNA overlaps filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		dnaClustersFilename = dnaClustersFilenameArg.getValue();
		rnaClustersFilename = rnaClustersFilenameArg.getValue();
		overlapsFilename = overlapsFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cerr << "Reading DNA clusters" << endl;
	
	IntegerTable dnaClusters;	
	ReadClusters(dnaClustersFilename, dnaClusters);
	
	cerr << "Reading RNA clusters" << endl;
	
	IntegerTable rnaClusters;	
	ReadClusters(rnaClustersFilename, rnaClusters);
	
	cerr << "Reading DNA/RNA overlaps" << endl;
	
	IntegerPairVec overlaps;
	ReadOverlaps(overlapsFilename, overlaps);
	
	cerr << "Calculating connected components" << endl;

	QuaterpartiteEdges edges(dnaClusters, overlaps, rnaClusters);
	
	IntegerTable dnaClusterComponents;
	IntegerTable rnaClusterComponents;
	edges.CalculateComponents(dnaClusterComponents, rnaClusterComponents);
	
	cerr << "Writing components" << endl;

	for (int componentIndex = 0; componentIndex < dnaClusterComponents.size(); componentIndex++)
	{
		cout << dnaClusterComponents[componentIndex].size() << "\t" << rnaClusterComponents[componentIndex].size() << endl;
	}
}

