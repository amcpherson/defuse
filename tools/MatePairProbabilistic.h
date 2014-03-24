/*
 *  MatePairProbabilistic.h
 *
 *  Created by Andrew McPherson.
 *
 */

#ifndef MATEPAIRPROBABILISTIC_H_
#define MATEPAIRPROBABILISTIC_H_

#include "IMatePairClusterer.h"
#include "Clusters.h"
#include "MatePairProbCalc.h"

#include <vector>

using namespace std;

class MatePairProbabilistic : public IMatePairClusterer
{
public:
	void Initialize(const Parameters& parameters);
	void DoClustering(const MatePairVec& matePairs, Clusters& clusters) const;
	
private:
	Clusters mClusters;
	MatePairProbCalc mProbCalc;
	int mMinClusterSize;
};

#endif
