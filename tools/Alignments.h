/*
 *  Alignments.h
 *  linkexons
 *
 *  Created by Andrew McPherson on 09-09-12.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ALIGNMENTS_H_
#define ALIGNMENTS_H_

#include "Common.h"
#include "Indexer.h"
#include "AlignmentStream.h"

#include <vector>
#include <iostream>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

using namespace std;
using namespace boost;

struct Alignment
{
	ReadID readID;
	RefStrand refStrand;
	Region region;
};

class Alignments
{
public:	
	Alignments();
	
	bool Read(IAlignmentStream* alignmentStream);
	
	void RemoveRedundant();
	void RemoveUnPaired();
	void FilterEnsemblTranscripts();
	void IndexByReadID();
	
	int GetNumAlignments() const;
	const Alignment& GetAlignment(int alignmentIndex) const;
	const IntegerVec& GetAlignments(const ReadID& readID) const;
	
	static NameIndex FragmentNames;
	static NameIndex ReferenceNames;
	
private:
	typedef vector<Alignment> AlignmentVec;
	typedef vector<Alignment>::iterator AlignmentVecIter;
	typedef vector<Alignment>::const_iterator AlignmentVecConstIter;
	
	typedef vector<IntegerVec> ReadAlignmentTable;
	typedef vector<IntegerVec>::iterator ReadAlignmentTableIter;
	typedef vector<IntegerVec>::const_iterator ReadAlignmentTableConstIter;
	
	AlignmentVec mAlignments;
	ReadAlignmentTable mReadAlignments[2];
};

#endif
