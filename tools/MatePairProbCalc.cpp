/*
 *  MatePairProbCalc.cpp
 *  findtransbreaks
 *
 *  Created by Andrew McPherson on 30/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "MatePairProbCalc.h"
#include "DebugCheck.h"

#include <math.h>
#include <iostream>

using namespace std;


const int numStdDev = 5;

MatePairProbCalc::MatePairProbCalc()
 : mFragmentMean(0), mFragmentStdDev(0), mReadLength(0)
{}

void MatePairProbCalc::Initialize(int fragmentMean, int fragmentStdDev, int readLength, double precision)
{
	mFragmentMean = fragmentMean;
	mFragmentStdDev = fragmentStdDev;
	mReadLength = readLength;
	
	mMu = mFragmentMean - 2 * mReadLength;
	
	double z = 0.0;

	int fragmentMin = mMu - numStdDev * mFragmentStdDev;
	int fragmentMax = mMu + numStdDev * mFragmentStdDev;

	for (int d1 = -fragmentMax; d1 <= fragmentMax; d1++)
	{
		for (int d2 = -fragmentMax; d2 <= fragmentMax; d2++)
		{
			z += CalculateUnNormalized(d1, d2);
		}
	}

	mZ = z;
	
	DoubleVec probs;
	for (int d1 = -fragmentMax; d1 <= fragmentMax; d1++)
	{
		for (int d2 = -fragmentMax - d1; d2 <= fragmentMax - d1; d2++)
		{
			double p = 0.0;
			
			p += CalculateUnNormalized(d1, d2) / mZ;

			probs.push_back(p);
		}
	}
	
	sort(probs.begin(), probs.end());
	reverse(probs.begin(), probs.end());
			
	double probthreshold = 0.0;
	double probsum = 0.0;
	DoubleVecConstIter probIter = probs.begin();
	while (probIter != probs.end() && probsum < precision)
	{
		probthreshold = *probIter;
		probsum += *probIter;
				
		probIter++;
	}
	
	mThreshold = probthreshold;
}

double MatePairProbCalc::CalculateUnNormalized(int d1, int d2) const
{
	int vMin = max(0, max(d1, max(d2, d1 + d2)));
	int vMax = ((int)(mMu + numStdDev * mFragmentStdDev));

	double p = 0.0;
	for (int v = vMin; v <= vMax; v++)
	{
		double pA = normalpdf(v, mMu, mFragmentStdDev);
		double pB = normalpdf(v - d1 - d2, mMu, mFragmentStdDev);

		p += pA * pB;		
	}
	
	return p;
}

double MatePairProbCalc::GetThreshold() const
{
	return mThreshold;
}

double MatePairProbCalc::GetDistanceThreshold() const
{
	return -log(mThreshold);
}

double MatePairProbCalc::Calculate(const Region& a1, const Region& a2, const Region& b1, const Region& b2) const
{
	int d1 = b1.end - a1.end;
	int d2 = b2.end - a2.end;

	return CalculateUnNormalized(d1,d2) / mZ;
}

void MatePairProbCalc::DistanceMatrix(const RegionVec& a1, const RegionVec& a2, DoubleTable& matrix) const
{
	DebugCheck(a1.size() == a2.size());
	
	matrix.resize(a1.size());

	for (int i = 0; i < (int)a1.size(); i++)
	{
		matrix[i].resize(a1.size());

		for (int j = 0; j < (int)a1.size(); j++)
		{
			matrix[i][j] = -log(Calculate(a1[i], a2[i], a1[j], a2[j]) + 0.000000001);
		}
	}
}




