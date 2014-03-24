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
	unordered_set<IntegerVec> nonredundant;
	
	for (int clusterIndex = 0; clusterIndex != mClusters.size(); clusterIndex++)
	{
		sort(mClusters[clusterIndex].begin(), mClusters[clusterIndex].end());
		nonredundant.insert(mClusters[clusterIndex]);
	}
	
	Clear();
	
	for (unordered_set<IntegerVec>::const_iterator clusterIter = nonredundant.begin(); clusterIter != nonredundant.end(); clusterIter++)
	{
		AddCluster(*clusterIter);
	}
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

