/*
 *  Clusters.cpp
 *
 *  Created by Andrew McPherson on 09-10-08.
 *
 */

#include "Clusters.h"
#include "Common.h"
#include "DebugCheck.h"

#include <vector>
#include <algorithm>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;


Clusters::Clusters()
{}

void Clusters::OpenCluster()
{
	mCurrentCluster.clear();
}

void Clusters::AddToCluster(int element)
{
	mCurrentCluster.push_back(element);
}

int Clusters::OpenClusterSize() const
{
	return mCurrentCluster.size();
}

int Clusters::CloseCluster()
{
	return AddCluster(mCurrentCluster);
}

int Clusters::AddCluster(const IntegerVec& cluster)
{
	int clusterIndex = mClusters.size();

	mClusters.push_back(cluster);
	
	return clusterIndex;
}

void Clusters::Merge(Clusters& other)
{
	for (int clusterIndex = 0; clusterIndex != other.mClusters.size(); clusterIndex++)
	{
		mClusters.push_back(IntegerVec());
		swap(mClusters.back(), other.mClusters[clusterIndex]);
	}
	other.mClusters.clear();
}

void Clusters::RemoveRedundant()
{
	for (IntegerTableIter clusterIter = mClusters.begin(); clusterIter != mClusters.end(); clusterIter++)
	{
		sort(clusterIter->begin(), clusterIter->end());
	}
	
	unordered_set<int> redundant;
	for (int cluster1Index = 0; cluster1Index < mClusters.size(); cluster1Index++)
	{
		for (int cluster2Index = 0; cluster2Index < mClusters.size(); cluster2Index++)
		{
			const IntegerVec& cluster1 = mClusters[cluster1Index];
			const IntegerVec& cluster2 = mClusters[cluster2Index];
			
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
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		if (redundant.find(clusterIndex) == redundant.end())
		{
			newClusters.push_back(mClusters[clusterIndex]);
		}
	}
	
	swap(newClusters, mClusters);
}

int Clusters::GetNumClusters() const
{
	return mClusters.size();
}

const IntegerTable& Clusters::GetClusters() const
{
	return mClusters;
}

const IntegerVec& Clusters::GetCluster(int clusterIndex) const
{
	return mClusters[clusterIndex];
}

void Clusters::CreateElementClusterMap(ElementClusterMap& elementClusters) const
{
	// Find the maximum element
	int maxElement = -1;
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		const IntegerVec& cluster = mClusters[clusterIndex];
		
		for (IntegerVecConstIter elementIter = cluster.begin(); elementIter != cluster.end(); elementIter++)
		{
			int element = *elementIter;
			
			DebugCheck(element >= 0);
			
			maxElement = max(maxElement, element);
		}
	}
	
	// Resize to fit all elements
	elementClusters.clear();
	elementClusters.resize(maxElement + 1);
	
	// Add elements to map
	for (int clusterIndex = 0; clusterIndex < mClusters.size(); clusterIndex++)
	{
		const IntegerVec& cluster = mClusters[clusterIndex];
		
		for (IntegerVecConstIter elementIter = cluster.begin(); elementIter != cluster.end(); elementIter++)
		{
			int element = *elementIter;
			
			elementClusters[element].push_back(clusterIndex);
		}
	}
}

void Clusters::Clear()
{
	mClusters.clear();
	mCurrentCluster.clear();	
}

