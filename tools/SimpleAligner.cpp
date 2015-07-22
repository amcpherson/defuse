/*
 *  SimpleAligner.cpp
 *
 *  Created by Andrew McPherson on 05/09/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "SimpleAligner.h"
#include "DebugCheck.h"

#include <iostream>
#include <list>

using namespace std;


SimpleAligner::SimpleAligner(int matchScore, int misMatchScore, int gapScore)
: mMatchScore(matchScore), mMisMatchScore(misMatchScore), mGapScore(gapScore)
{
}

int SimpleAligner::Align(const string& reference, const string& sequence)
{
	int matrixLength = reference.size() + 1;
	int matrixHeight = sequence.size() + 1;

	mMatrix.SetSize(matrixLength, matrixHeight);
	
	int overallMaxScore = 0;
	
	for (int i = 0; i < matrixLength; i++) 
	{
		int seq1Pos = i - 1;
		
		for (int j = 0; j < matrixHeight; j++)
		{
			int seq2Pos = j - 1;
			
			if (j == 0)
			{
				mMatrix(i,j) = 0;
			}
			else if (i == 0)
			{
				mMatrix(i,j) = mMatrix(i,j-1) + mGapScore;
			}
			else
			{
				int matchScore = mMatrix(i-1,j-1) + ((reference[seq1Pos] == sequence[seq2Pos]) ? mMatchScore : mMisMatchScore);
				int gapRefScore = mMatrix(i-1,j) + mGapScore;
				int gapReadScore = mMatrix(i,j-1) + mGapScore;
				int maxScore = max(matchScore,max(gapRefScore,gapReadScore));
				
				overallMaxScore = max(overallMaxScore, maxScore);
				
				mMatrix(i,j) = maxScore;
			}
		}
	}
	
	return overallMaxScore;
}

