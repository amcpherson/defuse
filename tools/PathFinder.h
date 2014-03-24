/*
 *  PathFinder.h
 *  findseqpath
 *
 *  Created by Andrew McPherson on 11/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PATHFINDER_H_
#define PATHFINDER_H_

#include "Common.h"

#include <string>

using namespace std;

class PathFinder
{
public:
	PathFinder(int k, int maxLength, bool findAll, int minCount) : mK(k), mMaxLength(maxLength), mFindAll(findAll), mMinCount(minCount) {}

	bool FindPath(const SequenceVec& sequences, const IntegerVec& leftSeq, const IntegerVec& rightSeq);

	const SequenceVec& GetSequences() { return mSequences; }
	const IntegerVec& GetMinCounts() { return mMinCounts; }

private:
	int mK;
	int mMaxLength;
	bool mFindAll;
	int mMinCount;
	SequenceVec mSequences;
	IntegerVec mMinCounts;
};

#endif

