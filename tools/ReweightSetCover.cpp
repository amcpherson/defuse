/*
 *  ReweightSetCover.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-09-02.
 *
 */

#include "ReweightSetCover.h"
#include "DebugCheck.h"

#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_set.hpp>

using namespace boost;
using namespace boost::bimaps;
using namespace std;

void CreateSizeVec(const IntegerTable& clusters, IntegerVec& sizes)
{
	for (IntegerTableConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		sizes.push_back(clusterIter->size());
	}
}

void CreateElementClusterMap(const IntegerTable& clusters, IntegerTable& elementClusterMap)
{
	for (int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
	{
		const IntegerVec& cluster = clusters[clusterIndex];
		for (IntegerVecConstIter elementIter = cluster.begin(); elementIter != cluster.end(); elementIter++)
		{
			elementClusterMap[*elementIter].push_back(clusterIndex);
		}
	}
}

ReweightSetCover::ReweightSetCover(const IntegerTable& dnaClusters, const IntegerTable& rnaClusters, const IntegerPairVec& overlap, float dnaWeight, float rnaWeight, float bothWeight)
: mDNAClusters(dnaClusters), mRNAClusters(rnaClusters), mOverlap(overlap), mDNAWeight(dnaWeight), mRNAWeight(rnaWeight), mBothWeight(bothWeight)
{
	// Find max fragment indices
	int maxDNAFragmentIndex = FindMaxElement(dnaClusters);
	int maxRNAFragmentIndex = FindMaxElement(rnaClusters);
	
	// Create dna/rna cluster size vectors
	CreateSizeVec(dnaClusters, mDNAClusterSize);
	CreateSizeVec(rnaClusters, mRNAClusterSize);
	
	// Create dna/rna fragment to cluster maps
	mDNAFragmentsToClusters = IntegerTable(maxDNAFragmentIndex + 1);
	mRNAFragmentsToClusters = IntegerTable(maxRNAFragmentIndex + 1);
	CreateElementClusterMap(dnaClusters, mDNAFragmentsToClusters);
	CreateElementClusterMap(rnaClusters, mRNAFragmentsToClusters);
	
	// Initialize assigned RNA/DNA fragments
	mAssignedDNAFragments = IntegerVec(maxDNAFragmentIndex + 1, 0);
	mAssignedRNAFragments = IntegerVec(maxRNAFragmentIndex + 1, 0);
	
	// Create cluster pairs
	// Create cluster pair sizes vector
	IntegerVec overlapDNACluster(dnaClusters.size(), 0);
	IntegerVec overlapRNACluster(rnaClusters.size(), 0);	
	mClusterPairs = overlap;
	for (IntegerPairVecConstIter overlapIter = overlap.begin(); overlapIter != overlap.end(); overlapIter++)
	{
		overlapDNACluster[overlapIter->first] = 1;
		overlapRNACluster[overlapIter->second] = 1;
	}
	
	// Add DNA clusters with no overlap
	for (int clusterIndex = 0; clusterIndex < overlapDNACluster.size(); clusterIndex++)
	{
		if (!overlapDNACluster[clusterIndex])
		{
			mClusterPairs.push_back(IntegerPair(clusterIndex,-1));
		}
	}
	
	// Add RNA clusters with no overlap
	for (int clusterIndex = 0; clusterIndex < overlapRNACluster.size(); clusterIndex++)
	{
		if (!overlapRNACluster[clusterIndex])
		{
			mClusterPairs.push_back(IntegerPair(-1,clusterIndex));
		}
	}
	
	// Initialize cluster weight map
	// Create dna/rna clusters to cluster pairs map
	mClusterPairEfficiencyMap.clear();
	mDNAClustersToClusterPairs = IntegerTable(dnaClusters.size());
	mRNAClustersToClusterPairs = IntegerTable(rnaClusters.size());
	for (int clusterPairIndex = 0; clusterPairIndex < mClusterPairs.size(); clusterPairIndex++)
	{
		int dnaClusterIndex = mClusterPairs[clusterPairIndex].first;
		int rnaClusterIndex = mClusterPairs[clusterPairIndex].second;
		
		// Add to dna/rna clusters to cluster pairs map
		if (dnaClusterIndex >= 0)
		{
			mDNAClustersToClusterPairs[dnaClusterIndex].push_back(clusterPairIndex);
		}
		if (rnaClusterIndex >= 0)
		{
			mRNAClustersToClusterPairs[rnaClusterIndex].push_back(clusterPairIndex);
		}
		
		// Initialize cluster efficiency
		mClusterPairEfficiencyMap.push_back(FloatBimap::value_type(clusterPairIndex, CalculateEfficiency(clusterPairIndex)));
	}
}

void ReweightSetCover::Calculate(IntegerTable& dnaSolutionClusters, IntegerTable& rnaSolutionClusters)
{
	// Solution clusters indexed similar to input clusters
	dnaSolutionClusters = IntegerTable(mDNAClusters.size());
	rnaSolutionClusters = IntegerTable(mRNAClusters.size());
	
	// Iterate until no more clusters found
	while (mClusterPairEfficiencyMap.size() > 0)
	{
		// Find the next cluster with the most support
		int nextClusterPairIndex = mClusterPairEfficiencyMap.right.rbegin()->second;
		int dnaClusterIndex = mClusterPairs[nextClusterPairIndex].first;
		int rnaClusterIndex = mClusterPairs[nextClusterPairIndex].second;
		
		// All remaining clusters are empty
		if (mClusterPairEfficiencyMap.right.rbegin()->first == 0.0f)
		{
			break;
		}
		
		// Assign remaining dna fragments in this cluster
		if (dnaClusterIndex >= 0)
		{
			Reweight(dnaClusterIndex, mDNAClusters, mDNAFragmentsToClusters, mDNAClustersToClusterPairs, mAssignedDNAFragments, mDNAClusterSize, dnaSolutionClusters);
		}
		
		if (rnaClusterIndex >= 0)
		{
			Reweight(rnaClusterIndex, mRNAClusters, mRNAFragmentsToClusters, mRNAClustersToClusterPairs, mAssignedRNAFragments, mRNAClusterSize, rnaSolutionClusters);
		}
	}		
}

float ReweightSetCover::CalculateEfficiency(int clusterPairIndex)
{
	int dnaClusterIndex = mClusterPairs[clusterPairIndex].first;
	int rnaClusterIndex = mClusterPairs[clusterPairIndex].second;
	
	// Calculate cluster pair weight
	float clusterPairWeight = 0.0;
	if (dnaClusterIndex >= 0 && mDNAClusterSize[dnaClusterIndex] > 0 && rnaClusterIndex >= 0 && mRNAClusterSize[rnaClusterIndex] > 0)
	{
		clusterPairWeight = mBothWeight;
	}
	else if (dnaClusterIndex >= 0 && mDNAClusterSize[dnaClusterIndex] > 0)
	{
		clusterPairWeight = mDNAWeight;
	}
	else if (rnaClusterIndex >= 0 && mRNAClusterSize[rnaClusterIndex] > 0)
	{
		clusterPairWeight = mRNAWeight;
	}
	else
	{
		return 0.0f;
	}
	
	// Calculate cluster pair size
	float clusterPairSize = mDNAClusterSize[dnaClusterIndex] + mRNAClusterSize[rnaClusterIndex];
	
	// Efficiency as elements per cost
	return clusterPairSize / clusterPairWeight;
}

void ReweightSetCover::Reweight(int clusterIndex, const IntegerTable& clusters, const IntegerTable& fragmentsToClusters, const IntegerTable& clustersToClusterPairs, IntegerVec& assignedFragments, IntegerVec& clusterSize, IntegerTable& solutionClusters)
{
	// Iterate through fragments and unassign from other cluster pairs
	for (IntegerVecConstIter fragmentIter = clusters[clusterIndex].begin(); fragmentIter != clusters[clusterIndex].end(); fragmentIter++)
	{
		int fragmentIndex = *fragmentIter;
		if (assignedFragments[fragmentIndex] == 0)
		{
			// Assign to this cluster
			solutionClusters[clusterIndex].push_back(fragmentIndex);
			
			// Flag as assigned
			assignedFragments[fragmentIndex] = 1;
			
			// Iterate through all other clusters containing this fragment
			for (IntegerVecConstIter otherClusterIter = fragmentsToClusters[fragmentIndex].begin(); otherClusterIter != fragmentsToClusters[fragmentIndex].end(); otherClusterIter++)
			{
				int otherClusterIndex = *otherClusterIter;
				
				// Subtract from the number of mate pairs in the cluster
				clusterSize[otherClusterIndex]--;
				DebugCheck(clusterSize[otherClusterIndex] >= 0);
				
				// Iterate through other cluster pairs connecting this cluster
				for (IntegerVecConstIter otherClusterPairIter = clustersToClusterPairs[otherClusterIndex].begin(); otherClusterPairIter != clustersToClusterPairs[otherClusterIndex].end(); otherClusterPairIter++)
				{
					int otherClusterPairIndex = *otherClusterPairIter;
					
					// Recalculate the total score of the cluster
					FloatBimap::left_iterator iter = mClusterPairEfficiencyMap.left.begin() + otherClusterPairIndex;
					DebugCheck(otherClusterPairIndex == iter->first);
					
					bool successful_modify = mClusterPairEfficiencyMap.left.replace_data(iter, CalculateEfficiency(otherClusterPairIndex));
					DebugCheck(successful_modify);
				}
			}
		}
	}
}

