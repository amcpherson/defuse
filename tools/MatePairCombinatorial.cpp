/*
 *  MatePairCombinatorial.cpp
 *
 *  Created by Andrew McPherson.
 *
 */

#include "MatePairCombinatorial.h"
#include "Common.h"
#include "DebugCheck.h"

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <boost/unordered_set.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>

using namespace std;
using namespace boost;
using namespace boost::bimaps;


void MatePairCombinatorial::Initialize(const Parameters& parameters)
{
	double precision = normalpdf(0, 0, parameters.fragmentStdDev);
	int x = 1;
	while (x < 10 * parameters.fragmentStdDev && parameters.precision > precision)
	{
		precision += 2 * normalpdf(x, 0, parameters.fragmentStdDev);
		x++;
	}
	
	mMaxFragmentLength = parameters.fragmentMean + x;
	mMaxFragmentLengthDiff = 2 * x;

	mMinClusterSize = parameters.minClusterSize;
	mResolution = parameters.resolution;
}

void MatePairCombinatorial::StrandRemap(const Region& region, int strand, Region& remapped) const
{
	Region remappedTemp;
	
	remappedTemp.start = (strand == PlusStrand) ? region.start : -region.end;
	remappedTemp.end = (strand == PlusStrand) ? region.end : -region.start;
	
	remapped = remappedTemp;
}

void MatePairCombinatorial::CalcBreakRegion(const Region& alignment, Region& breakRegion) const
{
	breakRegion.start = alignment.end + 1;
	breakRegion.end = 2 * alignment.start - alignment.end - 2 + mMaxFragmentLength;
}

void MatePairCombinatorial::CreateBreakRegion(int matePairIndex, const Region& alignment, FeatureRegion& breakRegion) const
{
	Region breakRegionTemp;
	
	CalcBreakRegion(alignment, breakRegionTemp);
	
	breakRegion.matePairIndex = matePairIndex;
	breakRegion.start = breakRegionTemp.start;
	breakRegion.end = breakRegionTemp.end;
}

void MatePairCombinatorial::CreateLengthRange(int matePairIndex, const Region& alignment1, const Region& alignment2, FeatureRegion& breakRegion) const
{
	breakRegion.matePairIndex = matePairIndex;
	breakRegion.start = alignment1.start + alignment2.start;
	breakRegion.end = breakRegion.start + mMaxFragmentLengthDiff - 1;
}

void CountFilter(IntegerPairVec& elements, int minCount)
{
	if (elements.size() < 1)
	{
		return;
	}
	
	sort(elements.begin(), elements.end());
	
	IntegerPairVec newElements;
	IntegerPairVec::iterator currentElementIter = elements.begin();
	int currentCount = 0;
	for (IntegerPairVec::iterator elementIter = elements.begin() + 1; elementIter != elements.end(); elementIter++)
	{
		if (*currentElementIter < *elementIter)
		{
			if (currentCount >= minCount)
			{
				newElements.push_back(*currentElementIter);
			}
			currentCount = 0;
		}
		
		currentElementIter = elementIter;
		currentCount++;
	}
	
	if (currentCount >= minCount)
	{
		newElements.push_back(*currentElementIter);
	}
	
	swap(elements,newElements);
}

void Intersect(const IntegerVec& vec1, const IntegerVec& vec2, IntegerVec& result)
{
	unordered_set<int> set1(vec1.begin(), vec1.end());
	unordered_set<int> intersect;
	for (IntegerVecConstIter element2Iter = vec2.begin(); element2Iter != vec2.end(); element2Iter++)
	{
		if (set1.find(*element2Iter) != set1.end())
		{
			intersect.insert(*element2Iter);
		}
	}
	
	result = IntegerVec(intersect.begin(), intersect.end());	
}

void MatePairCombinatorial::CreateOverlapClusters(FeatureRegionVec& regions, IntegerTable& overlapClusters) const
{
	sort(regions.begin(), regions.end(), FeatureRegionStartLessThan());
	
	multiset<FeatureRegion,FeatureRegionEndLessThan> currentCluster;
	
	int currentRegionIndex = 0;
	while (currentRegionIndex < (int)regions.size())
	{
		bool inserted = false;
		while (currentRegionIndex < (int)regions.size() && (currentCluster.size() == 0 || regions[currentRegionIndex].start <= currentCluster.begin()->end))
		{
			inserted = true;
			currentCluster.insert(regions[currentRegionIndex]);
			currentRegionIndex++;
		}
		
		if (inserted && currentCluster.size() != 0)
		{
			overlapClusters.push_back(IntegerVec());
			for (multiset<FeatureRegion,FeatureRegionEndLessThan>::const_iterator regionIter = currentCluster.begin(); regionIter != currentCluster.end(); regionIter++)
			{
				overlapClusters.back().push_back(regionIter->matePairIndex);
			}
		}
		
		int currentEnd = currentCluster.begin()->end;
		while (!currentCluster.empty() && currentCluster.begin()->end < currentEnd + mResolution)
		{
			currentCluster.erase(currentCluster.begin());
		}
	}
}

void MatePairCombinatorial::CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, MatePair& matePair) const
{
	matePair.id = id;
	StrandRemap(alignment1, strand1, matePair.alignment1);
	StrandRemap(alignment2, strand2, matePair.alignment2);
}

void MatePairCombinatorial::CalcBreakRegion(const Region& alignment, int strand, Region& breakRegion) const
{
	Region remapAlignment;
	StrandRemap(alignment, strand, remapAlignment);
	CalcBreakRegion(remapAlignment, breakRegion);
	StrandRemap(breakRegion, strand, breakRegion);
}

void MatePairCombinatorial::SelectFeatureRegions(const IntegerVec& indices, const FeatureRegionVec& features, FeatureRegionVec& selected) const
{
	for (IntegerVecConstIter iter = indices.begin(); iter != indices.end(); iter++)
	{
		selected.push_back(features[*iter]);
	}
}

void MatePairCombinatorial::RemoveRedundant(IntegerTable& clusters) const
{
	for (IntegerTableIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		sort(clusterIter->begin(), clusterIter->end());
	}
	
	unordered_set<int> redundant;
	for (int cluster1Index = 0; cluster1Index < clusters.size(); cluster1Index++)
	{
		for (int cluster2Index = 0; cluster2Index < clusters.size(); cluster2Index++)
		{
			const IntegerVec& cluster1 = clusters[cluster1Index];
			const IntegerVec& cluster2 = clusters[cluster2Index];
			
			bool clusterIncludes = includes(cluster1.begin(), cluster1.end(), cluster2.begin(), cluster2.end());
			bool clusterEqual = clusterIncludes && cluster1.size() == cluster2.size();
			
			if (clusterIncludes && !clusterEqual)
			{
				redundant.insert(cluster2Index);
			}
			
			if (clusterEqual && cluster1Index < cluster2Index)
			{
				redundant.insert(cluster2Index);				
			}
		}
	}
	
	IntegerTable newClusters;
	for (int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
	{
		if (redundant.find(clusterIndex) == redundant.end())
		{
			newClusters.push_back(clusters[clusterIndex]);
		}
	}
	
	swap(newClusters, clusters);
}

void MatePairCombinatorial::DoClustering(const MatePairVec& matePairs, IntegerTable& clusters) const
{
	if (matePairs.size() < mMinClusterSize)
	{
		return;
	}

	FeatureRegionVec breakRegions1;
	FeatureRegionVec breakRegions2;
	FeatureRegionVec lengthRanges;
	
	for (int matePairIndex = 0; matePairIndex < (int)matePairs.size(); matePairIndex++)
	{
		FeatureRegion breakRegion1;
		CreateBreakRegion(matePairIndex, matePairs[matePairIndex].alignment1, breakRegion1);

		breakRegions1.push_back(breakRegion1);

		FeatureRegion breakRegion2;
		CreateBreakRegion(matePairIndex, matePairs[matePairIndex].alignment2, breakRegion2);
		
		breakRegions2.push_back(breakRegion2);

		FeatureRegion lengthRange;
		CreateLengthRange(matePairIndex, matePairs[matePairIndex].alignment1, matePairs[matePairIndex].alignment2, lengthRange);

		lengthRanges.push_back(lengthRange);
	}
	
	IntegerTable breakRegionClusters1;
	CreateOverlapClusters(breakRegions1, breakRegionClusters1);
cout << breakRegionClusters1.size() << endl;
		
	IntegerTable breakRegionClusters2;
	for (IntegerTableConstIter clusterIter = breakRegionClusters1.begin(); clusterIter != breakRegionClusters1.end(); clusterIter++)
	{		
		FeatureRegionVec clusterBreakRegions2;
		SelectFeatureRegions(*clusterIter, breakRegions2, clusterBreakRegions2);
		
		CreateOverlapClusters(clusterBreakRegions2, breakRegionClusters2);
	}
cout << breakRegionClusters2.size() << endl;
	
	RemoveRedundant(breakRegionClusters2);
cout << breakRegionClusters2.size() << endl;
	
	IntegerTable lengthRangeClusters;
	for (IntegerTableConstIter clusterIter = breakRegionClusters2.begin(); clusterIter != breakRegionClusters2.end(); clusterIter++)
	{		
		FeatureRegionVec clusterLengthRange;
		SelectFeatureRegions(*clusterIter, lengthRanges, clusterLengthRange);
		
		CreateOverlapClusters(clusterLengthRange, lengthRangeClusters);
	}
	
cout << lengthRangeClusters.size() << endl;
	RemoveRedundant(lengthRangeClusters);
cout << lengthRangeClusters.size() << endl;
	
	clusters.insert(clusters.end(), lengthRangeClusters.begin(), lengthRangeClusters.end());	
}

