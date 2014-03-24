/*
 *  SimpleAligner.h
 *
 *  Created by Andrew McPherson on 05/09/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SIMPLEALIGNER_H_
#define SIMPLEALIGNER_H_

#include "Common.h"
#include "Matrix.h"

#include <string>
#include <vector>

using namespace std;

class SimpleAligner
{
public:
	SimpleAligner(int matchScore, int misMatchScore, int gapScore);
	
	int Align(const string& reference, const string& sequence);
	
private:
	int mMatchScore;
	int mMisMatchScore;
	int mGapScore;
	
	Matrix<int> mMatrix;
};

#endif

