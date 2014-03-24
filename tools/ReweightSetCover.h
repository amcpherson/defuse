/*
 *  ReweightSetCover.h
 *  tools
 *
 *  Created by Andrew McPherson on 10-09-02.
 *
 */

#ifndef REWEIGHTSETCOVER_H_
#define REWEIGHTSETCOVER_H_

#include "Common.h"

#include <vector>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>

using namespace boost;
using namespace boost::bimaps;
using namespace std;

class ReweightSetCover
{
public:
	ReweightSetCover(const IntegerTable& dnaClusters, const IntegerTable& rnaClusters, const IntegerPairVec& overlap, float dnaWeight, float rnaWeight, float bothWeight);
	void Calculate(IntegerTable& dnaSolutionClusters, IntegerTable& rnaSolutionClusters);
	
private:
	float CalculateEfficiency(int clusterPairIndex);	
	void Reweight(int clusterIndex, const IntegerTable& clusters, const IntegerTable& fragmentsToClusters, const IntegerTable& clustersToClusterPairs, IntegerVec& assignedFragments, IntegerVec& clusterSize, IntegerTable& solutionClusters);
	
	typedef bimap<vector_of<int>,multiset_of<float> > FloatBimap;
	
	const IntegerTable& mDNAClusters;
	const IntegerTable& mRNAClusters;
	const IntegerPairVec& mOverlap;
	float mDNAWeight;
	float mRNAWeight;
	float mBothWeight;
	
	IntegerVec mDNAClusterSize;
	IntegerVec mRNAClusterSize;
	IntegerTable mDNAFragmentsToClusters;
	IntegerTable mRNAFragmentsToClusters;
	IntegerVec mAssignedDNAFragments;
	IntegerVec mAssignedRNAFragments;
	IntegerPairVec mClusterPairs;
	FloatBimap mClusterPairEfficiencyMap;
	IntegerTable mDNAClustersToClusterPairs;
	IntegerTable mRNAClustersToClusterPairs;	
};

#endif
