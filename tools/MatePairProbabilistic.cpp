/*
 *  MatePairProbabilistic.cpp
 *
 *  Created by Andrew McPherson.
 *
 */

#include "MatePairProbabilistic.h"
#include "Clusters.h"
#include "MatePairProbCalc.h"
#include "HierarchicalClusterer.h"

#include <iostream>
#include <vector>

using namespace std;


void MatePairProbabilistic::Initialize(const Parameters& parameters)
{
	//KLUDGE
	mProbCalc.Initialize(parameters.fragmentMean, parameters.fragmentStdDev, /*parameters.readLength*/ 50, parameters.precision);
	mMinClusterSize = parameters.minClusterSize;

	double precision = normalpdf(0, 0, parameters.fragmentStdDev);
	int x = 1;
	while (x < 10 * parameters.fragmentStdDev && parameters.precision > precision)
	{
		precision += 2 * normalpdf(x, 0, parameters.fragmentStdDev);
		x++;
	}
	
	mMaxFragmentLength = parameters.fragmentMean + x;
}

void MatePairProbabilistic::StrandRemap(const Region& region, int strand, Region& remapped) const
{
	Region remappedTemp;
	
	remappedTemp.start = (strand == PlusStrand) ? region.start : -region.end;
	remappedTemp.end = (strand == PlusStrand) ? region.end : -region.start;
	
	remapped = remappedTemp;
}

void MatePairProbabilistic::CalcBreakRegion(const Region& alignment, Region& breakRegion) const
{
	breakRegion.start = alignment.end + 1;
	breakRegion.end = 2 * alignment.start - alignment.end - 2 + mMaxFragmentLength;
}

void MatePairProbabilistic::CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, MatePair& matePair) const
{
	matePair.id = id;
	StrandRemap(alignment1, strand1, matePair.alignment1);
	StrandRemap(alignment2, strand2, matePair.alignment2);
}

void MatePairProbabilistic::CalcBreakRegion(const Region& alignment, int strand, Region& breakRegion) const
{
	Region remapAlignment;
	StrandRemap(alignment, strand, remapAlignment);
	CalcBreakRegion(remapAlignment, breakRegion);
	StrandRemap(breakRegion, strand, breakRegion);
}

void MatePairProbabilistic::DoClustering(const MatePairVec& matePairs, Clusters& clusters) const
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

