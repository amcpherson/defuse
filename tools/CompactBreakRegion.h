/*
 *  CompactBreakRegion.h
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef COMPACTBREAKREGION_H_
#define COMPACTBREAKREGION_H_

#include "Common.h"
#include "Indexer.h"

#include <vector>
#include <string>

using namespace std;

struct CompactBreakRegion
{
	int clusterID;
	RefStrand refStrand;
	int start;
	int end;
};

typedef vector<CompactBreakRegion> BrRegVec;
typedef vector<CompactBreakRegion>::iterator BrRegVecIter;
typedef vector<CompactBreakRegion>::const_iterator BrRegVecConstIter;

bool BreakRegionStartLessThan(const CompactBreakRegion& br1, const CompactBreakRegion& br2);
void ReadBreakRegions(const string& breakRegionsFilename, NameIndex& referenceNames, BrRegVec& breakRegion);


#endif

