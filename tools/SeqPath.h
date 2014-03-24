/*
 *  SeqPath.h
 *
 *  Created by Andrew McPherson on 09-10-08.
 *
 */

#ifndef SEQPATH_H_
#define SEQPATH_H_

#include "Common.h"

#include <vector>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

class SeqPath
{
public:
	SeqPath();
	
	const char& Get(int index) const;
	const IntegerVec& LookBack(int index) const;
	
	int MergedLength() const;
	bool Empty() const;
	
	void Reverse();

	void InitializePaths(const RegionTable& paths);
	void SetReverseComplement(bool reverseComplement);
	void InitializeSequence(const string& sequence);

private:
	RegionVec mMergedPaths;
	IntegerTable mLookBack;
	string mSequence;
	bool mReverseComplement;
};

#endif

