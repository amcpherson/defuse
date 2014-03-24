/*
 *  HierarchicalClusterer.cpp
 *
 *  Created by Andrew McPherson on 1/27/10.
 *
 */

#include "HierarchicalClusterer.h"

#include <iostream>
#include <vector>
#include <list>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>

using namespace std;
using namespace boost;
using namespace boost::bimaps;


struct SortedPair
{
	SortedPair(int a, int b) : first(min(a,b)), second(max(a,b)) {}
	
	int first;
	int second;
};

inline bool operator==(const SortedPair& sp1, const SortedPair& sp2)
{
	return sp1.first == sp2.first && sp1.second == sp2.second;
}

inline size_t hash_value(const SortedPair& sp)
{
	size_t seed = 0;
	hash_combine(seed, sp.first);
	hash_combine(seed, sp.second);
    return seed;
}

typedef bimap<multiset_of<double>,unordered_set_of<SortedPair> > DistancesBimap;


void HierarchicalClusterer::DoClustering(const DoubleTable& distances, double threshold)
{
	mClusters.clear();
	
	if (distances.size() < 1)
	{
		return;
	}
	
	IntegerTable clusters;
	DistancesBimap sortedDistances;
	list<int> clusterIndices;
	
	for (int i = 0; i < (int)distances.size(); i++)
	{
		clusters.push_back(IntegerVec());
		clusters.back().push_back(i);
		
		clusterIndices.push_back(i);

		for (int j = i + 1; j < (int)distances.size(); j++)
		{
			sortedDistances.insert(DistancesBimap::value_type(distances[i][j],SortedPair(i,j)));
		}
	}

	while (!sortedDistances.empty() && sortedDistances.left.begin()->first < threshold)
	{
		SortedPair toMerge = sortedDistances.left.begin()->second;
		
		double sizeFirst = clusters[toMerge.first].size();
		double sizeSecond = clusters[toMerge.second].size();
		double sizeMerged = sizeFirst + sizeSecond;
		
		// Create a merged cluster
		int mergedClusterIndex = clusters.size();
		clusters.push_back(IntegerVec());
		clusters.back().insert(clusters.back().end(), clusters[toMerge.first].begin(), clusters[toMerge.first].end());
		clusters.back().insert(clusters.back().end(), clusters[toMerge.second].begin(), clusters[toMerge.second].end());
		
		// Iterate through cluster indices
		list<int>::iterator clusterIndexIter = clusterIndices.begin();
		while (clusterIndexIter != clusterIndices.end())
		{
			int clusterIndex = *clusterIndexIter;
			
			if (clusterIndex == toMerge.first || clusterIndex == toMerge.second)
			{
				// Remove indices of clusters being merged
				clusterIndexIter = clusterIndices.erase(clusterIndexIter);
			}
			else
			{
				// Find distances from each merge cluster to current cluster
				DistancesBimap::right_iterator distToFirstIter = sortedDistances.right.find(SortedPair(toMerge.first,clusterIndex));
				DistancesBimap::right_iterator distToSecondIter = sortedDistances.right.find(SortedPair(toMerge.second,clusterIndex));
				
				double distToFirst = distToFirstIter->second;
				double distToSecond = distToSecondIter->second;
				
				// Calculate distance between current cluster and merged cluster
				double distToMerged = (sizeFirst * distToFirst + sizeSecond * distToSecond) / sizeMerged;
				
				// Erase distances to each merge cluster
				sortedDistances.right.erase(distToFirstIter);
				sortedDistances.right.erase(distToSecondIter);
				
				// Insert distance between current cluster and merged cluster
				sortedDistances.insert(DistancesBimap::value_type(distToMerged,SortedPair(mergedClusterIndex,clusterIndex)));
			
				clusterIndexIter++;
			}
		}
		
		// Add index of merged cluster
		clusterIndices.push_back(mergedClusterIndex);
		
		// Erase the distance between the merged clusters
		sortedDistances.left.erase(sortedDistances.left.begin());
	}
	
	// Add clusters to main cluster list
	for (list<int>::iterator clusterIndexIter = clusterIndices.begin(); clusterIndexIter != clusterIndices.end(); clusterIndexIter++)
	{
		int clusterIndex = *clusterIndexIter;
		
		if (clusters[clusterIndex].size() < 1)
		{
			continue;
		}
		
		mClusters.push_back(IntegerVec());
		swap(mClusters.back(),clusters[clusterIndex]);
	}
}

const IntegerTable& HierarchicalClusterer::GetClusters()
{
	return mClusters;
}

