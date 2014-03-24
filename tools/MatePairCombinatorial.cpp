/*
 *  MatePairCombinatorial.cpp
 *
 *  Created by Andrew McPherson.
 *
 */

#include "MatePairCombinatorial.h"
#include "Clusters.h"
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

void MatePairCombinatorial::CreateOverlapClusters(FeatureRegionVec& regions, Clusters& overlapClusters) const
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
		
		if (inserted)
		{
			overlapClusters.OpenCluster();
			for (multiset<FeatureRegion,FeatureRegionEndLessThan>::const_iterator regionIter = currentCluster.begin(); regionIter != currentCluster.end(); regionIter++)
			{
				overlapClusters.AddToCluster(regionIter->matePairIndex);
			}
			overlapClusters.CloseCluster();
		}
		
		currentCluster.erase(currentCluster.begin());
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

void MatePairCombinatorial::DoClustering(const MatePairVec& matePairs, Clusters& clusters) const
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
	
	Clusters breakRegionClusters1;
	Clusters breakRegionClusters2;
	Clusters lengthRangeClusters;
	
	CreateOverlapClusters(breakRegions1, breakRegionClusters1);
	CreateOverlapClusters(breakRegions2, breakRegionClusters2);
	CreateOverlapClusters(lengthRanges, lengthRangeClusters);
	
	ElementClusterMap clusterElements1;
	ElementClusterMap clusterElements2;
	ElementClusterMap clusterElements3;
	
	breakRegionClusters1.CreateElementClusterMap(clusterElements1);
	breakRegionClusters2.CreateElementClusterMap(clusterElements2);
	lengthRangeClusters.CreateElementClusterMap(clusterElements3);

	ClusterElementMap clusterElements;
	
	// Divide elements into clusters according to 3 features
	for (int matePairIndex = 0; matePairIndex < (int)matePairs.size(); matePairIndex++)
	{
		int matePairID = matePairs[matePairIndex].id;

		const IntegerVec& matePairClusters1 = clusterElements1[matePairIndex];
		const IntegerVec& matePairClusters2 = clusterElements2[matePairIndex];
		const IntegerVec& matePairClusters3 = clusterElements3[matePairIndex];
		
		for (IntegerVecConstIter matePairCluster1Iter = matePairClusters1.begin(); matePairCluster1Iter != matePairClusters1.end(); matePairCluster1Iter++)
		{
			for (IntegerVecConstIter matePairCluster2Iter = matePairClusters2.begin(); matePairCluster2Iter != matePairClusters2.end(); matePairCluster2Iter++)
			{
				for (IntegerVecConstIter matePairCluster3Iter = matePairClusters3.begin(); matePairCluster3Iter != matePairClusters3.end(); matePairCluster3Iter++)
				{
					ClusterTuple clusterTuple;
					clusterTuple.cluster1 = *matePairCluster1Iter;
					clusterTuple.cluster2 = *matePairCluster2Iter;
					clusterTuple.cluster3 = *matePairCluster3Iter;
					
					clusterElements[clusterTuple].push_back(matePairID);
				}
			}
		}
	}
	
	// Add clusters
	Clusters combClusters;
	for (ClusterElementMapConstIter clusterElementsIter = clusterElements.begin(); clusterElementsIter != clusterElements.end(); clusterElementsIter++)
	{
		combClusters.AddCluster(clusterElementsIter->second);
	}
	
	// Remove redundant clusters resulting from the combinatorial algorithm
	combClusters.RemoveRedundant();
	
	// Merge the new clusters into the main set of clusters
	clusters.Merge(combClusters);
}

