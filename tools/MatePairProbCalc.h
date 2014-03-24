/*
 *  MatePairProbCalc.h
 *  findtransbreaks
 *
 *  Created by Andrew McPherson on 30/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATEPAIRPROBCALC_H_
#define MATEPAIRPROBCALC_H_

#include "Common.h"

#include <vector>

using namespace std;

class MatePairProbCalc
{
public:
	MatePairProbCalc();
	
	void Initialize(int fragmentMean, int fragmentStdDev, int readLength, double precision);
	double Calculate(const Region& a1, const Region& a2, const Region& b1, const Region& b2) const;
	void DistanceMatrix(const RegionVec& a1, const RegionVec& a2, DoubleTable& matrix) const;
	double GetThreshold() const;
	double GetDistanceThreshold() const;
	
private:
	double CalculateUnNormalized(int d1, int d2) const;

	int mFragmentMean;
	int mFragmentStdDev;
	int mReadLength;
	int mMu;
	double mZ;
	double mThreshold;
};


#endif

