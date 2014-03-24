/*
 *  AlignmentIndex.h
 *
 *  Created by Andrew McPherson on 02/09/10.
 *
 */

#ifndef ALIGNMENTINDEX_H_
#define ALIGNMENTINDEX_H_

#include "Common.h"
#include "Indexer.h"

#include "api/BamReader.h"

#include <map>
#include <string>
#include <iostream>

using namespace std;
using namespace BamTools;

class AlignmentIndex
{
public:
	void Open(const string& bamFilename);
	void Find(const string& reference, int strand, int start, int end, NameIndex& referenceNames, CompAlignVec& alignments) const;
	void Close();
	
private:
	mutable BamReader mReader;
};

#endif
