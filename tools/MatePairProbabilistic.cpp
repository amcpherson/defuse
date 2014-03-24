/*
 *  MatePairProbabilistic.cpp
 *
 *  Created by Andrew McPherson.
 *
 */

#include "MatePairProbabilistic.h"
#include "MatePairProbCalc.h"
#include "HierarchicalClusterer.h"

#include <iostream>
#include <vector>

using namespace std;


void MatePairProbabilistic::Initialize(const Parameters& parameters)
{
	mProbCalc.Initialize(parameters.fragmentMean, parameters.fragmentStdDev, parameters.readLength, parameters.precision);
	mMinClusterSize = parameters.minClusterSize;
}

void MatePairProbabilistic::DoClustering(const MatePairVec& matePairs, IntegerTable& clusters) const
{
	if (matePairs.size() < mMinClusterSize)
	{
		return;
	}

	cout << "Clustering " << matePairs.size() << " mate pairs" << endl;

	RegionVec alignments1;
	RegionVec alignments2;
	
	for (int matePairIndex = 0; matePairIndex < (int)matePairs.size(); matePairIndex++)
	{
		alignments1.push_back(matePairs[matePairIndex].alignment1);
		alignments2.push_back(matePairs[matePairIndex].alignment2);		
	}

	DoubleTable distanceMatrix;
	mProbCalc.DistanceMatrix(alignments1, alignments2, distanceMatrix);
	
	double threshold = mProbCalc.GetDistanceThreshold();

	HierarchicalClusterer clusterer;
	clusterer.DoClustering(distanceMatrix,threshold);
	
	const IntegerTable& hierClusters = clusterer.GetClusters();
	
	for (IntegerTableConstIter clusterIter = hierClusters.begin(); clusterIter != hierClusters.end(); clusterIter++)
	{
		if (clusterIter->size() < mMinClusterSize)
		{
			continue;
		}
		
		clusters.OpenCluster();

		for (IntegerVecConstIter elementIter = clusterIter->begin(); elementIter != clusterIter->end(); elementIter++)
		{
			int matePairID = matePairs[*elementIter].id;
			clusters.AddToCluster(matePairID);
		}

		clusters.CloseCluster();
	}
}

