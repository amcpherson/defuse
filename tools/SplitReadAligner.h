/*
 *  SplitReadAligner.h
 *  findbreaks
 *
 *  Created by Andrew McPherson on 30/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SPLITREADALIGNER_H_
#define SPLITREADALIGNER_H_

#include "Common.h"
#include "Matrix.h"

#include <string>
#include <vector>

using namespace std;

struct SplitReadAlignment
{
	IntegerPair refSplit;
	IntegerPair readSplit;
	int score;
	int score1;
	int score2;
	IntegerPairVec matches1;
	IntegerPairVec matches2;
};

typedef vector<SplitReadAlignment> SplitReadAlignVec;

class SplitReadAligner
{
public:
	SplitReadAligner(int matchScore, int misMatchScore, int gapScore, bool endGaps, int minSplitScore, const string& reference1, const string& reference2);
	
	void Align(const string& read);
	void GetAlignments(SplitReadAlignVec& alignments, int minScore, bool forceSplits, bool firstOnly, bool backtrace);
	
private:
	void FillMatrix(const string& seq1, const string& seq2, Matrix<int>& matrix, Matrix<Cell>& backTrace);
	void FindMaxRowEntry(const Matrix<int>& matrix, int j, int minAccepted, int& max);
	void FindMaxRowEntry(const Matrix<int>& matrix, int j, int minAccepted, int& max, IntegerVec& columns);
	void BackTracePath(const Matrix<Cell>& backTrace, const Cell& startCell, IntegerPairVec& matches);
	void ReverseMatches(int refLength, int readLength, IntegerPairVec& matches);

	typedef vector<Cell> CellVec;
	typedef vector<Cell>::iterator CellVecIter;
	typedef vector<Cell>::const_iterator CellVecConstIter;
	
	int mMatchScore;
	int mMisMatchScore;
	int mGapScore;
	bool mEndGaps;
	int mMinSplitScore;
	
	string mReference1;
	string mReference2;
	string mRead1;
	string mRead2;
	Matrix<int> mMatrix1;
	Matrix<int> mMatrix2;
	Matrix<Cell> mBackTrace1;
	Matrix<Cell> mBackTrace2;
};

#endif

