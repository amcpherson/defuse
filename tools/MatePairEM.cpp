/*
 *  MatePairEM.cpp
 *
 *  Created by Andrew McPherson.
 *
 */

#include "MatePairEM.h"
#include "DebugCheck.h"

#include "asa136.H"

#include <iostream>
#include <vector>
#include <numeric>
#include <limits>
#include "asa241.H"

using namespace std;


template <typename t>
void printvec(t begin, t end)
{
	cout.precision(30);
	cout << "c(";
	
	t iter = begin;
	while (iter != end)
	{
		cout << *iter;
		iter++;
		
		if (iter != end)
		{
			cout << ", ";
		}
	}
	
	cout << ")" << endl;
}

void MatePairEM::Initialize(const Parameters& parameters)
{
	mFragmentMean = parameters.fragmentMean;
	mFragmentStdDev = parameters.fragmentStdDev;
	mMinClusterSize = parameters.minClusterSize;

	double x = -parameters.fragmentStdDev * r8_normal_01_cdf_inverse((1 - parameters.precision) / 2);
	mMinProbability = normalpdf(x, 0, parameters.fragmentStdDev);

	mMaxFragmentLength = parameters.fragmentMean + 3 * parameters.fragmentStdDev;

	mKMeansIter = 1000;
	mLambda = 0.1;
	mTolerance = 0.001;
	mKMax = 10;
}

void MatePairEM::CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, MatePair& matePair) const
{
	matePair.id = id;
	StrandRemap(alignment1, strand1, matePair.alignment1);
	StrandRemap(alignment2, strand2, matePair.alignment2);
}

void MatePairEM::CalcBreakRegion(const Region& alignment, int strand, Region& breakRegion) const
{
	Region remapAlignment;
	StrandRemap(alignment, strand, remapAlignment);
	CalcBreakRegion(remapAlignment, breakRegion);
	StrandRemap(breakRegion, strand, breakRegion);
}

void MatePairEM::StrandRemap(const Region& region, int strand, Region& remapped) const
{
	Region remappedTemp;
	
	remappedTemp.start = (strand == PlusStrand) ? region.start : -region.end;
	remappedTemp.end = (strand == PlusStrand) ? region.end : -region.start;
	
	remapped = remappedTemp;
}

void MatePairEM::CalcBreakRegion(const Region& alignment, Region& breakRegion) const
{
	breakRegion.start = alignment.end + 1;
	breakRegion.end = 2 * alignment.start - alignment.end - 2 + (int)mMaxFragmentLength;
}

double MatePairEM::PairProbability(double x, double y, double u, double a, double b) const
{
	return normalpdf(a + b - x - y, u, mFragmentStdDev) * exp(-mLambda * max(0.0, x - a) - mLambda * max(0.0, y - b));
}

double MatePairEM::LogLikelihood() const
{
	vector<vector<double> > exponents(mK);
	for (int j = 0; j < mK; j++)
	{
		exponents[j].resize(mN);
	}
	
	for (int i = 0; i < mN; i++)
	{
		for (int j = 0; j < mK; j++)
		{
			exponents[j][i] = -0.5 * pow((mA[j] + mB[j] - mX[i] - mY[i] - mU[i]) / mFragmentStdDev, 2.0) - mLambda * max(0.0, mX[i] - mA[j]) - mLambda * max(0.0, mY[i] - mB[j]);
		}
	}
	
	double LL = 0.0;
	for (int i = 0; i < mN; i++)
	{
		double maxexp = exponents[0][i];
		for (int j = 1; j < mK; j++)
		{
			maxexp = max(maxexp, exponents[j][i]);
		}
		
		double sum = 0.0;
		for (int j = 0; j < mK; j++)
		{
			sum += mW[j] * exp(exponents[j][i] - maxexp);
		}
				
		if (sum == 0.0)
		{
			LL = -numeric_limits<double>::max();
			break;
		}
		
		LL = LL + log(sum) + maxexp;
	}
	
	return LL;
}

void MatePairEM::UpdateResponsibilities()
{
	vector<vector<double> > exponents(mK);
	for (int j = 0; j < mK; j++)
	{
		exponents[j].resize(mN);
	}
	
	for (int i = 0; i < mN; i++)
	{
		for (int j = 0; j < mK; j++)
		{
			exponents[j][i] = -0.5 * pow((mA[j] + mB[j] - mX[i] - mY[i] - mU[i]) / mFragmentStdDev, 2.0) - mLambda * max(0.0, mX[i] - mA[j]) - mLambda * max(0.0, mY[i] - mB[j]);
		}
	}
	
	for (int i = 0; i < mN; i++)
	{
		int iXO = mToXO[i];
		int iYO = mToYO[i];
		
		double maxexp = exponents[0][i];
		for (int j = 1; j < mK; j++)
		{
			maxexp = max(maxexp, exponents[j][i]);
		}
		
		double norm = 0.0;
		for (int j = 0; j < mK; j++)
		{
			norm += mW[j] * exp(exponents[j][i] - maxexp);
		}
		
		DebugCheck(norm != 0.0);
		for (int j = 0; j < mK; j++)
		{
			mR[j][i] = mW[j] * exp(exponents[j][i] - maxexp) / norm;
			
			mRXO[j][iXO] = mR[j][i];
			mRYO[j][iYO] = mR[j][i];
		}
	}
}

void MatePairEM::UpdateMixWeights()
{
	for (int j = 0; j < mK; j++)
	{
		double NK = accumulate(mR[j].begin(), mR[j].end(), 0.0);
		mW[j] = NK / mN;
	}
}

bool MatePairEM::MaxLikelihood(const vector<double>& R, const vector<double>& RXO, const vector<double>& RYO, double& a, double& b) const
{
	vector<double> SX(mN);
	vector<double> SY(mN);
	
	partial_sum(RXO.begin(), RXO.end(), SX.begin());
	partial_sum(RYO.begin(), RYO.end(), SY.begin());

	int i = 0;
	int j = 0;
	vector<double> CX;
	vector<double> CY;
	vector<double> CS;
	
	CX.push_back(mXO[0]);
	CY.push_back(mYO[0]);
	CS.push_back(0.0);
	
	while (i < SX.size() && j < SY.size())
	{
		if (i + 1 < SX.size() && mXO[i] == mXO[i+1])
		{
			i++;
			continue;
		}
		
		if (j + 1 < SY.size() && mYO[j] == mYO[j+1])
		{
			j++;
			continue;
		}

		if (SX[i] == SY[j])
		{
			CX.push_back(mXO[i]);
			CY.push_back(mYO[j]);
			CS.push_back(SX[i]);
			
			if (i + 1 < SX.size() && j + 1 < SY.size())
			{
				CX.push_back(mXO[i+1]);
				CY.push_back(mYO[j+1]);
				CS.push_back(SX[i]);
			}
			
			i++;
			j++;
		}
		else if (SX[i] < SY[j])
		{
			CX.push_back(mXO[i]);
			CY.push_back(mYO[j]);
			CS.push_back(SX[i]);
			
			if (i + 1 < SX.size())
			{
				CX.push_back(mXO[i+1]);
				CY.push_back(mYO[j]);
				CS.push_back(SX[i]);
			}
			
			i++;
		}
		else
		{
			CX.push_back(mXO[i]);
			CY.push_back(mYO[j]);
			CS.push_back(SY[j]);
			
			if (j + 1 < SY.size())
			{
				CX.push_back(mXO[i]);
				CY.push_back(mYO[j+1]);
				CS.push_back(SY[j]);
			}
			
			j++;
		}
	}
	
	double NK = accumulate(R.begin(), R.end(), 0.0);
	
	if (NK == 0.0)
	{
		return false;
	}

	double RXYU = 0.0;
	for (int i = 0; i < mN; i++)
	{
		RXYU += R[i] * (mX[i] + mY[i] + mU[i]);
	}
	
	vector<double> partial(CS.size());
	for (int i = 0; i < CS.size(); i++)
	{
		partial[i] = (RXYU - NK * (CX[i] + CY[i])) / pow(mFragmentStdDev, 2.0) + mLambda * CS[i];
	}
	
	int minindex = 0;
	while (minindex < partial.size())
	{
		if (partial[minindex] > 0)
		{
			break;
		}

		minindex++;
	}
	
	double aplusb = (RXYU + pow(mFragmentStdDev, 2.0) * mLambda * CS[minindex]) / NK;
	
	if (minindex == 0)
	{
		double min_a = CX[minindex];
		double max_a = aplusb - CY[minindex];
		a = 0.5 * (min_a + max_a);
		b = aplusb - a;
	}
	else if (CS[minindex] != CS[minindex-1])
	{
		a = CX[minindex];
		b = CY[minindex];
	}
	else
	{
		double min_a = max(CX[minindex], aplusb - CY[minindex-1]);
		double max_a = min(CX[minindex-1], aplusb - CY[minindex]);
		a = 0.5 * (min_a + max_a);
		b = aplusb - a;
	}
	
	return true;
}

bool MatePairEM::SelectKKZ(int k, vector<double>& A, vector<double>& B)
{
	DebugCheck(k <= mN);
	
	A.clear();
	B.clear();
	
	double L2max = mX[0] * mY[0];
	int iL2max = 0;
	for (int i = 1; i < mN; i++)
	{
		double L2 = mX[i] * mY[i];
		
		if (L2 > L2max)
		{
			iL2max = i;
			L2max = L2;
		}
	}
	
	A.push_back(mX[iL2max]);
	B.push_back(mY[iL2max]);
		
	while (A.size() < k)
	{
		vector<double> DistMin(mN);
		for (int i = 0; i < mN; i++)
		{
			double minDist = pow(mX[i] - A[0], 2.0) + pow(mY[i] - B[0], 2.0);
			for (int j = 1; j < A.size(); j++)
			{
				double dist = pow(mX[i] - A[j], 2.0) + pow(mY[i] - B[j], 2.0);
				minDist = min(minDist, dist);
			}
			
			DistMin[i] = minDist;
		}
		
		double DistsMax = DistMin[0];
		int iDistsMax = 0;
		for (int i = 0; i < mN; i++)
		{
			if (DistMin[i] > DistsMax)
			{
				DistsMax = DistMin[i];
				iDistsMax = i;
			}
		}
		
		if (DistsMax == 0.0)
		{
			return false;
		}
		
		A.push_back(mX[iDistsMax]);
		B.push_back(mY[iDistsMax]);
	}
	
	return true;
}

bool MatePairEM::ExpectationMaximization(double& ll)
{
	DebugCheck(mK <= mN);

	mR.resize(mK);
	mRXO.resize(mK);
	mRYO.resize(mK);
	for (int j = 0; j < mK; j++)
	{
		mR[j].resize(mN);
		mRXO[j].resize(mN);
		mRYO[j].resize(mN);
	}
	mW.resize(mK);
	mA.resize(mK);
	mB.resize(mK);
	
	if (mK == 1 || mK == mN)
	{
		for (int j = 0; j < mK; j++)
		{
			mR[j] = vector<double>(mN, 1.0 / mK);
			mRXO[j] = vector<double>(mN, 1.0 / mK);
			mRYO[j] = vector<double>(mN, 1.0 / mK);
		}
	}
	else
	{
		vector<double> px;
		vector<double> py;
		if (!SelectKKZ(mK, px, py))
		{
			return false;
		}
		
		vector<double> a;
		a.insert(a.begin(), mX.begin(), mX.end());
		a.insert(a.begin(), mY.begin(), mY.end());
		
		vector<double> c;
		c.insert(c.begin(), px.begin(), px.end());
		c.insert(c.begin(), py.begin(), py.end());
		
		vector<int> ic1(mN);
		vector<int> nc(mK);
		vector<double> wss(mK);
		int ifault;
		
		kmns(&a.front(), mN, 2, &c.front(), mK, &ic1.front(), &nc.front(), mKMeansIter, &wss.front(), &ifault);
		DebugCheck(ifault != 1);
		DebugCheck(ifault != 3);
		
		for (int i = 0; i < mN; i++)
		{
			int iXO = mToXO[i];
			int iYO = mToYO[i];
			
			for (int j = 0; j < mK; j++)
			{
				mR[j][i] = (j == ic1[i] - 1) ? 1 : 0;
				
				mRXO[j][iXO] = mR[j][i];
				mRYO[j][iYO] = mR[j][i];
			}
		}
	}
	
	double lastLikelihood;
	bool lastLikelihoodValid = false;
	while (true)
	{
		for (int j = 0; j < mK; j++)
		{
			double a;
			double b;
			if (MaxLikelihood(mR[j], mRXO[j], mRYO[j], a, b))
			{
				mA[j] = a;
				mB[j] = b;
			}
		}
		
		UpdateMixWeights();
		
		double likelihood = LogLikelihood();
		if (lastLikelihoodValid && abs(likelihood - lastLikelihood) < mTolerance)
		{
			break;
		}
		
		if (lastLikelihoodValid && likelihood == -numeric_limits<double>::max())
		{
			return false;
		}
		
		DebugCheck(!lastLikelihoodValid || (likelihood / lastLikelihood < 1.0000001));
		
		lastLikelihood = likelihood;
		lastLikelihoodValid = true;

		UpdateResponsibilities();
	}
	
	ll = lastLikelihood;
	
	return true;
}

struct MatePairInfo
{
	int i;
	double x;
	double y;
	double u;
};

typedef vector<MatePairInfo> MatePairInfoVec;

bool XGreaterThan(const MatePairInfo& m1, const MatePairInfo& m2)
{
	return (m1.x > m2.x);
}

bool YGreaterThan(const MatePairInfo& m1, const MatePairInfo& m2)
{
	return (m1.y > m2.y);
}





//int xlist[]={-993917, -972421, -960630, -993910, -993918, -994491, -993907, -993915, -993907, -960618, -993906, -993901, -988339, -960618, -987217, -993906, -993918, -993899, -993906, -988317, -993932, -993906, -993910, -993915, -993910, -960618, -993918, -994133, -993925, -991084, -988321, -972421, -991084, -988321, -993920, -993917, -993903, -993906};
//int ylist[]={991037, 975193, 947665, 991035, 991049, 992004, 991035, 991036, 991035, 947653, 991036, 991043, 985018, 947653, 989891, 991036, 991034, 991026, 991036, 985599, 991054, 991036, 991026, 991035, 991035, 947653, 991035, 992071, 991035, 987200, 985012, 975193, 987200, 985012, 991039, 991035, 991039, 991036};

//void MatePairEM::DoClustering(const MatePairVec& matePairs, IntegerTable& clusters)
//{
//	if (matePairs.size() < mMinClusterSize)
//	{
//		return;
//	}
//	
//	mN = sizeof(xlist) / sizeof(int);
//	
//	MatePairInfoVec info(mN);
//	for (int matePairIndex = 0; matePairIndex < mN; matePairIndex++)
//	{
//		info[matePairIndex].i = matePairIndex;
//		info[matePairIndex].x = xlist[matePairIndex];
//		info[matePairIndex].y = ylist[matePairIndex];
//		info[matePairIndex].u = mFragmentMean - 100;
//	}
void MatePairEM::DoClustering(const MatePairVec& matePairs, IntegerTable& clusters)
{
	if (matePairs.size() < mMinClusterSize)
	{
		return;
	}
	
	mN = matePairs.size();
	
	MatePairInfoVec info(mN);
	for (int matePairIndex = 0; matePairIndex < (int)matePairs.size(); matePairIndex++)
	{
		info[matePairIndex].i = matePairIndex;
		info[matePairIndex].x = matePairs[matePairIndex].alignment1.end;
		info[matePairIndex].y = matePairs[matePairIndex].alignment2.end;
		info[matePairIndex].u = mFragmentMean - Length(matePairs[matePairIndex].alignment1) - Length(matePairs[matePairIndex].alignment2);
	}
	
	mX.resize(mN);
	mY.resize(mN);
	mU.resize(mN);
	for (int matePairIndex = 0; matePairIndex < (int)info.size(); matePairIndex++)
	{
		mX[matePairIndex] = info[matePairIndex].x;
		mY[matePairIndex] = info[matePairIndex].y;
		mU[matePairIndex] = info[matePairIndex].u;
	}
	
	sort(info.begin(), info.end(), XGreaterThan);
	
	mXO.resize(mN);
	mToXO.resize(mN);
	for (int sortIndex = 0; sortIndex < (int)info.size(); sortIndex++)
	{
		mXO[sortIndex] = info[sortIndex].x;
		mToXO[info[sortIndex].i] = sortIndex;
	}
	
	sort(info.begin(), info.end(), YGreaterThan);
	
	mYO.resize(mN);
	mToYO.resize(mN);
	for (int sortIndex = 0; sortIndex < (int)info.size(); sortIndex++)
	{
		mYO[sortIndex] = info[sortIndex].y;
		mToYO[info[sortIndex].i] = sortIndex;
	}

	double minBIC;
	bool minBICValid = false;
	int kMinBIC = 1;
	for (mK = 1; mK <= min(mKMax, mN); mK++)
	{
		double logLikelihood;
		if (!ExpectationMaximization(logLikelihood))
		{
			continue;
		}

		double BIC = -2.0 * logLikelihood + mK * 2.0 * log(mN);
		
		if (!minBICValid || BIC < minBIC)
		{
			minBIC = BIC;
			kMinBIC = mK;
			minBICValid = true;
		}
	}
	
	double logLikelihood;
	mK = kMinBIC;
	if (!ExpectationMaximization(logLikelihood))
	{
		cerr << "Error: No clusters" << endl;
		return;
	}
	
	for (int j = 0; j < mK; j++)
	{
		vector<int> cluster;
		for (int i = 0; i < mN; i++)
		{
			double prob = PairProbability(mX[i], mY[i], mU[i], mA[j], mB[j]);
			
			if (prob > mMinProbability)
			{
				cluster.push_back(i);
			}
		}
		
		if (cluster.size() >= mMinClusterSize)
		{
			clusters.push_back(vector<int>());
			swap(clusters.back(), cluster);
		}
	}	
}



