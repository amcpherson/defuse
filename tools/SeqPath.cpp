/*
 *  SeqPath.cpp
 *
 *  Created by Andrew McPherson on 09-10-08.
 *
 */

#include "SeqPath.h"
#include "Common.h"
#include "DebugCheck.h"

#include <vector>
#include <algorithm>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;


SeqPath::SeqPath()
{
}

void MergeRegions(const RegionVec& regions, RegionVec& merged, unordered_map<int,int>& mergedPositions)
{
	if (regions.size() == 0)
	{
		return;
	}
	
	Region merge;
	merge.start = regions.front().start;
	merge.end = regions.front().end;
	for (RegionVecConstIter regionIter = regions.begin(); regionIter != regions.end(); regionIter++)
	{
		if (regionIter->start > merge.end + 1)
		{
			merged.push_back(merge);
			
			merge.start = regionIter->start;
			merge.end = regionIter->end;
		}
		else
		{
			merge.end = max(merge.end, regionIter->end);
		}
	}

	merged.push_back(merge);
	
	int localPosition = 1;
	for (RegionVecConstIter mergedIter = merged.begin(); mergedIter != merged.end(); mergedIter++)
	{
		for (int globalPosition = mergedIter->start; globalPosition <= mergedIter->end; globalPosition++)
		{
			mergedPositions[globalPosition] = localPosition;
			localPosition++;
		}
	}
}

void CombineExons(const RegionTable& regions, RegionVec& combined)
{
	for (RegionTableConstIter regionsIter = regions.begin(); regionsIter != regions.end(); regionsIter++)
	{
		combined.insert(combined.end(), regionsIter->begin(), regionsIter->end());
	}
}

void SeqPath::InitializePaths(const RegionTable& paths)
{
	mMergedPaths.clear();
	
	RegionVec combinedPaths;
	CombineExons(paths, combinedPaths);
	
	unordered_map<int,int> mergedPositions;
	MergeRegions(combinedPaths, mMergedPaths, mergedPositions);
	
	mLookBack.clear();
	mLookBack.resize(mergedPositions.size());
	
	for (RegionTableConstIter pathIter = paths.begin(); pathIter != paths.end(); pathIter++)
	{
		int previousMergedIndex = -1;
		for (RegionVecConstIter regionIter = pathIter->begin(); regionIter != pathIter->end(); regionIter++)
		{
			for (int globalPosition = regionIter->start; globalPosition <= regionIter->end; globalPosition++)
			{
				int mergedPosition = mergedPositions[globalPosition];
				int mergedIndex = mergedPosition - 1;
				
				if (previousMergedIndex != -1)
				{
					mLookBack[mergedIndex].push_back(previousMergedIndex);
				}
				
				previousMergedIndex = mergedIndex;
			}
		}
	}
	
	for (int mergedIndex = 0; mergedIndex < mLookBack.size(); mergedIndex++)
	{
		unordered_set<int> uniqueIndices(mLookBack[mergedIndex].begin(), mLookBack[mergedIndex].end());
		mLookBack[mergedIndex] = IntegerVec(uniqueIndices.begin(), uniqueIndices.end());
	}
}

void SeqPath::SetReverseComplement(bool reverseComplement)
{
	mReverseComplement = reverseComplement;
}

void SeqPath::InitializeSequence(const string& sequence)
{
	mSequence = "";
	for (RegionVecConstIter regionIter = mMergedPaths.begin(); regionIter != mMergedPaths.end(); regionIter++)
	{
		DebugCheck(regionIter->start >= 1);
		DebugCheck(regionIter->end >= regionIter->start);
		DebugCheck(regionIter->end <= sequence.length());
		mSequence.append(sequence.substr(regionIter->start - 1, regionIter->end - regionIter->start + 1));
	}
	
	if (mReverseComplement)
	{
		Reverse();		
		ReverseComplement(mSequence);
		reverse(mSequence.begin(), mSequence.end());
	}
}

const char& SeqPath::Get(int index) const
{
	return mSequence[index];
}

const IntegerVec& SeqPath::LookBack(int index) const
{
	return mLookBack[index];
}

int SeqPath::MergedLength() const
{
	return mSequence.length();
}

bool SeqPath::Empty() const
{
	return mSequence.empty();
}

void SeqPath::Reverse()
{
	reverse(mSequence.begin(), mSequence.end());
	
	IntegerTable newLookBack(mSequence.size());
	for (int nucleotideIndex = 0; nucleotideIndex < mSequence.size(); nucleotideIndex++)
	{
		int backTo = mSequence.size() - nucleotideIndex - 1;
		for (IntegerVecConstIter lookBackIter = mLookBack[nucleotideIndex].begin(); lookBackIter != mLookBack[nucleotideIndex].end(); lookBackIter++)
		{
			int from = mSequence.size() - *lookBackIter - 1;
			
			newLookBack[from].push_back(backTo);
		}
	}
	
	swap(mLookBack, newLookBack);
}


