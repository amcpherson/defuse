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
	void CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, MatePair& matePair) const;
	void CalcBreakRegion(const Region& alignment, int strand, Region& breakRegion) const;
	void DoClustering(const MatePairVec& matePairs, Clusters& clusters) const;
	
private:
	void StrandRemap(const Region& region, int strand, Region& remapped) const;
	void CalcBreakRegion(const Region& alignment, Region& breakRegion) const;

	Clusters mClusters;
	MatePairProbCalc mProbCalc;
	int mMinClusterSize;
	int mMaxFragmentLength;
};

#endif
