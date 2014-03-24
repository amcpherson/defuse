/*
 *  IMatePairClusters.h
 *
 *  Created by Andrew McPherson.
 *
 */

#ifndef IMATEPAIRCLUSTERS_H_
#define IMATEPAIRCLUSTERS_H_

#include "Clusters.h"

#include <vector>

using namespace std;

class IMatePairClusterer
{
public:
	struct Parameters
	{
		int fragmentMean;
		int fragmentStdDev;
		double precision;
		int minClusterSize;
	};
	
	struct MatePair
	{
		int id;
		Region alignment1;
		Region alignment2;
	};
	
	typedef vector<MatePair> MatePairVec;
	
	virtual ~IMatePairClusterer() {}
	virtual void Initialize(const Parameters& parameters) = 0;
	virtual void CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, MatePair& matePair) const = 0;
	virtual void CalcBreakRegion(const Region& alignment, int strand, Region& breakRegion) const = 0;
	virtual void DoClustering(const MatePairVec& matePairs, Clusters& clusters) const = 0;
};

#endif
