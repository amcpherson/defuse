/*
 *  overlapclusters.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "Alignments.h"
#include "CompactBreakRegion.h"
#include "DebugCheck.h"
#include "IMatePairClusterer.h"
#include "MatePairCombinatorial.h"

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

void Print(const CompactBreakRegion& cbr)
{
	cout << cbr.clusterID << "\t" << cbr.refStrand.id << "\t" << cbr.start << "\t" << cbr.end << "\n";
}

typedef unordered_map<IntegerPair,BrRegVec> RefPairBreakRegionsMap;
typedef unordered_map<IntegerPair,BrRegVec>::iterator RefPairBreakRegionsMapIter;

void CreateRefPairBreakRegionsMap(BrRegVec& breakRegions, RefPairBreakRegionsMap& refPairMap)
{
	while (breakRegions.size() != 0)
	{
		CompactBreakRegion breakRegionEnd1 = breakRegions.back(); breakRegions.pop_back();
		CompactBreakRegion breakRegionEnd2 = breakRegions.back(); breakRegions.pop_back();
		
		if (breakRegionEnd1.clusterID != breakRegionEnd2.clusterID)
		{
			cerr << "Error: break regions should be ordered by cluster ID" << endl;		
			exit(1);		
		}
		
		if (breakRegionEnd1.refStrand.id > breakRegionEnd2.refStrand.id)
		{
			swap(breakRegionEnd1,breakRegionEnd2);
		}
		
		refPairMap[IntegerPair(breakRegionEnd1.refStrand.id,breakRegionEnd2.refStrand.id)].push_back(breakRegionEnd1);
		refPairMap[IntegerPair(breakRegionEnd1.refStrand.id,breakRegionEnd2.refStrand.id)].push_back(breakRegionEnd2);
	}
}

int main(int argc, char* argv[])
{
	string breakRegionsAFilename;
	string breakRegionsBFilename;
	string overlapFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Overlap between break regions tool");
		TCLAP::ValueArg<string> breakRegionsAFilenameArg("a","abreaks","Input Breakpoint Regions A Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> breakRegionsBFilenameArg("b","bbreaks","Input Breakpoint Regions B Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> overlapFilenameArg("o","overlap","Output Overlap Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		breakRegionsAFilename = breakRegionsAFilenameArg.getValue();
		breakRegionsBFilename = breakRegionsBFilenameArg.getValue();
		overlapFilename = overlapFilenameArg.getValue();
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
	
	cout << "Reading break regions" << endl;
	
	BrRegVec breakRegionsA;
	BrRegVec breakRegionsB;
	NameIndex referenceNames;

	ReadBreakRegions(breakRegionsAFilename, referenceNames, breakRegionsA);
	ReadBreakRegions(breakRegionsBFilename, referenceNames, breakRegionsB);
	
	if ((breakRegionsA.size() % 2) == 1)
	{
		cerr << "Error: break regions A has uneven number of break regions" << endl;		
		exit(1);		
	}
	
	if ((breakRegionsB.size() % 2) == 1)
	{
		cerr << "Error: break regions B has uneven number of break regions" << endl;		
		exit(1);		
	}
	
	cout << "Creating reference pair break regions map" << endl;

	RefPairBreakRegionsMap refPairBreakRegionsMapA;
	RefPairBreakRegionsMap refPairBreakRegionsMapB;
	
	CreateRefPairBreakRegionsMap(breakRegionsA, refPairBreakRegionsMapA);
	CreateRefPairBreakRegionsMap(breakRegionsB, refPairBreakRegionsMapB);
	
	// Open cluster pairs file
	ofstream overlapFile(overlapFilename.c_str());
	if (!overlapFile)
	{
		cerr << "Error: unable to write to overlap file" << endl;		
		exit(1);
	}
	
	cout << "Iterating reference pairs" << endl;
	
	// Iterate through break regions for each ref pair
	int numOverlaps = 0;
	for (RefPairBreakRegionsMapIter refPairIter = refPairBreakRegionsMapA.begin(); refPairIter != refPairBreakRegionsMapA.end(); refPairIter++)
	{
		const IntegerPair& refPair = refPairIter->first;
		if (refPairBreakRegionsMapB.find(refPair) == refPairBreakRegionsMapB.end())
		{
			continue;
		}
		
		BrRegVec& refPairBreakRegionsA = refPairBreakRegionsMapA[refPair];
		BrRegVec& refPairBreakRegionsB = refPairBreakRegionsMapB[refPair];

		sort(refPairBreakRegionsA.begin(), refPairBreakRegionsA.end(), BreakRegionStartLessThan);
		sort(refPairBreakRegionsB.begin(), refPairBreakRegionsB.end(), BreakRegionStartLessThan);
		
		IntegerPairVec clusterPairs;

		BrRegVecConstIter brRegAIter = refPairBreakRegionsA.begin();
		BrRegVecConstIter brRegBIter = refPairBreakRegionsB.begin();
		while (brRegAIter != refPairBreakRegionsA.end())
		{
			while (brRegBIter != refPairBreakRegionsB.end() && (brRegBIter->refStrand.id < brRegAIter->refStrand.id || brRegBIter->end < brRegAIter->start))
			{
				brRegBIter++;
			}
			
			BrRegVecConstIter brRegBOverlapIter = brRegBIter;
			while (brRegBOverlapIter != refPairBreakRegionsB.end() && (brRegBOverlapIter->refStrand.id == brRegAIter->refStrand.id && brRegBOverlapIter->start <= brRegAIter->end))
			{
				if (brRegBOverlapIter->end >= brRegAIter->start)
				{
					clusterPairs.push_back(IntegerPair(brRegAIter->clusterID,brRegBOverlapIter->clusterID));
				}
				
				brRegBOverlapIter++;
			}
			
			brRegAIter++;
		}
		
		sort(clusterPairs.begin(), clusterPairs.end());
		
		IntegerPairVecConstIter prevPairIter = clusterPairs.begin();
		for (IntegerPairVecConstIter pairIter = clusterPairs.begin(); pairIter != clusterPairs.end(); pairIter++)
		{
			if (pairIter != clusterPairs.begin() && *pairIter == *prevPairIter)
			{
				overlapFile << pairIter->first << "\t" << pairIter->second << endl;
				numOverlaps++;
			}
			
			prevPairIter = pairIter;
		}
	}
	
	overlapFile.close();

	cout << "Output " << numOverlaps << " overlaps" << endl;
}

