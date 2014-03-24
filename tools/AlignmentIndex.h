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
#include "bam.h"
#include "sam.h"

#include <map>
#include <string>
#include <iostream>

using namespace std;

class AlignmentIndex
{
public:
	AlignmentIndex();
	
	void Open(const string& bamFilename);
	void Find(const string& reference, int strand, int start, int end, NameIndex& referenceNames, CompAlignVec& alignments) const;
	void Close();
	
private:
	samfile_t* mBamFile;
	bam_index_t* mBamIndexFile;
	unordered_map<string,int> mReferenceIndexMap;
};

#endif
