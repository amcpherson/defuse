/*
 *  MatePairEM.h
 *
 *  Created by Andrew McPherson.
 *
 */

#ifndef MATEPAIREM_H_
#define MATEPAIREM_H_

#include "IMatePairClusterer.h"

#include <vector>

using namespace std;

class MatePairEM : public IMatePairClusterer
{
public:
	void Initialize(const Parameters& parameters);
	void CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, MatePair& matePair) const;
	void CalcBreakRegion(const Region& alignment, int strand, Region& breakRegion) const;
	void DoClustering(const MatePairVec& matePairs, IntegerTable& clusters);
	
private:
	void StrandRemap(const Region& region, int strand, Region& remapped) const;
	void CalcBreakRegion(const Region& alignment, Region& breakRegion) const;
	
	double PairProbability(double x, double y, double u, double a, double b) const;
	double LogLikelihood() const;
	void UpdateResponsibilities();
	bool MaxLikelihood(const vector<double>& R, const vector<double>& RXO, const vector<double>& RYO, double& a, double& b) const;
	void UpdateMixWeights();
	bool SelectKKZ(int k, vector<double>& A, vector<double>& B);
	bool ExpectationMaximization(double& ll);

	double mFragmentMean;
	double mFragmentStdDev;
	double mMaxFragmentLength;
	double mMinProbability;
	double mMinClusterSize;
	double mLambda;
	int mKMeansIter;
	double mTolerance;
	int mKMax;
	
	int mN;
	int mK;
	
	vector<double> mW;
	vector<double> mX;
	vector<double> mY;
	vector<double> mU;
	vector<double> mA;
	vector<double> mB;
	vector<vector<double> > mR;
	
	vector<double> mXO;
	vector<double> mYO;
	vector<vector<double> > mRXO;
	vector<vector<double> > mRYO;
	vector<int> mToXO;
	vector<int> mToYO;
};

#endif
