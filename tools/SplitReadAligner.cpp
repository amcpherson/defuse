/*
 *  SplitReadAligner.cpp
 *  findbreaks
 *
 *  Created by Andrew McPherson on 30/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "SplitReadAligner.h"
#include "DebugCheck.h"

#include <iostream>
#include <list>

using namespace std;


SplitReadAligner::SplitReadAligner(int matchScore, int misMatchScore, int gapScore, bool endGaps, int minSplitScore, const string& reference1, const string& reference2)
: mMatchScore(matchScore), mMisMatchScore(misMatchScore), mGapScore(gapScore), mEndGaps(endGaps), mMinSplitScore(minSplitScore), mReference1(reference1), mReference2(reference2)
{
	reverse(mReference2.begin(), mReference2.end());	
}

void SplitReadAligner::FillMatrix(const string& seq1, const string& seq2, Matrix<int>& matrix, Matrix<Cell>& backTrace)
{
	int matrixLength = seq1.size() + 1;
	int matrixHeight = seq2.size() + 1;

	matrix.Resize(matrixLength, matrixHeight);
	backTrace.Resize(matrixLength, matrixHeight);

	for (int i = 0; i < matrixLength; i++) 
	{
		int seq1Pos = i - 1;
		
		for (int j = 0; j < matrixHeight; j++)
		{
			int seq2Pos = j - 1;
			
			if (j == 0)
			{
				matrix(i,j) = 0;
			}
			else if (i == 0)
			{
				matrix(i,j) = matrix(i,j-1) + (mEndGaps ? 0 : mGapScore);
				backTrace(i,j) = Cell(i,j-1);
			}
			else
			{
				int matchScore = matrix(i-1,j-1) + ((seq1[seq1Pos] == seq2[seq2Pos]) ? mMatchScore : mMisMatchScore);
				int gapRefScore = matrix(i-1,j) + mGapScore;
				int gapReadScore = matrix(i,j-1) + mGapScore;
				int maxScore = max(matchScore,max(gapRefScore,gapReadScore));
								
				if (matchScore == maxScore)
				{
					backTrace(i,j) = Cell(i-1,j-1);
				}
				
				if (gapRefScore == maxScore)
				{
					backTrace(i,j) = Cell(i-1,j);
				}
				
				if (gapReadScore == maxScore)
				{
					backTrace(i,j) = Cell(i,j-1);
				}
				
				matrix(i,j) = maxScore;
			}
		}
	}				
}

void SplitReadAligner::Align(const string& read)
{
	mRead1 = read;
	
	mRead2 = read;
	reverse(mRead2.begin(), mRead2.end());
	
	FillMatrix(mReference1, mRead1, mMatrix1, mBackTrace1);
	FillMatrix(mReference2, mRead2, mMatrix2, mBackTrace2);
}

void SplitReadAligner::FindMaxRowEntry(const Matrix<int>& matrix, int j, int minAccepted, int& max)
{
	max = 0;
	
	for (int i = 0; i < matrix.GetLength(); i++) 
	{
		if (matrix(i,j) >= minAccepted && matrix(i,j) > max)
		{
			max = matrix(i,j);
		}
	}
}

void SplitReadAligner::FindMaxRowEntry(const Matrix<int>& matrix, int j, int minAccepted, int& max, IntegerVec& columns)
{
	max = 0;
	columns.clear();
	
	for (int i = 0; i < matrix.GetLength(); i++) 
	{
		if (matrix(i,j) >= minAccepted && matrix(i,j) > max)
		{
			max = matrix(i,j);
			columns.clear();
			columns.push_back(i);
		}
		else if (matrix(i,j) >= minAccepted && matrix(i,j) == max)
		{
			columns.push_back(i);
		}
	}
}

void SplitReadAligner::BackTracePath(const Matrix<Cell>& backTrace, const Cell& startCell, IntegerPairVec& matches)
{
	Cell cell = startCell;
	while (cell.j > 0)
	{
		Cell nextCell = backTrace(cell);
		
		int refPos = cell.i - 1;
		int readPos = cell.j - 1;

		if (cell.i - 1 == nextCell.i && cell.j - 1 == nextCell.j)
		{
			matches.push_back(IntegerPair(refPos,readPos));
		}
		
		cell = nextCell;
	}
	
	reverse(matches.begin(), matches.end());
}

void SplitReadAligner::ReverseMatches(int refLength, int readLength, IntegerPairVec& matches)
{
	for (IntegerPairVecIter pairIter = matches.begin(); pairIter != matches.end(); pairIter++)
	{
		pairIter->first = refLength - pairIter->first - 1;
		pairIter->second = readLength - pairIter->second - 1;
	}
	
	reverse(matches.begin(), matches.end());
}

void SplitReadAligner::GetAlignments(SplitReadAlignVec& alignments, int minScore, bool forceSplit, bool firstOnly, bool backTrace)
{
	int maxScore = 0;

	// Look for a non split alignment first unless we are forcing a split alignment
	if (!forceSplit)
	{
		// Look in matrix 1 for alignments that are not split
		for (int i = 0; i < mMatrix1.GetLength(); i++) 
		{
			for (int j = 0; j < mMatrix1.GetHeight(); j++)
			{
				if (mMatrix1(i,j) >= minScore && mMatrix1(i,j) > maxScore)
				{
					maxScore = mMatrix1(i,j);
				}
			}
		}
		
		// Look in matrix 2 for alignments that are not split
		for (int i = 0; i < mMatrix2.GetLength(); i++) 
		{
			for (int j = 0; j < mMatrix2.GetHeight(); j++)
			{
				if (mMatrix2(i,j) >= minScore && mMatrix2(i,j) > maxScore)
				{
					maxScore = mMatrix2(i,j);
				}
			}
		}
		
		// Return nothing if we found an unsplit alignment
		if (maxScore != 0)
		{
			return;
		}
	}
	
	// Find max score for split alignments
	IntegerVec alignedToRef1sMax;
	for (int alignedToRef1 = 0; alignedToRef1 <= mRead1.length(); alignedToRef1++)
	{
		int alignedToRef2 = mRead2.length() - alignedToRef1;
		
		DebugCheck(alignedToRef2 >= 0);
		
		int maxAcrossRow1 = 0;
		FindMaxRowEntry(mMatrix1, alignedToRef1, mMinSplitScore, maxAcrossRow1);
		
		int maxAcrossRow2 = 0;
		FindMaxRowEntry(mMatrix2, alignedToRef2, mMinSplitScore, maxAcrossRow2);
		
		int maxAcrossRows = maxAcrossRow1 + maxAcrossRow2;

		if (maxAcrossRows >= minScore && maxAcrossRows > maxScore)
		{
			maxScore = maxAcrossRows;
			
			alignedToRef1sMax.clear();
			alignedToRef1sMax.push_back(alignedToRef1);
		}
		else if (maxAcrossRows >= minScore && maxAcrossRows == maxScore)
		{
			alignedToRef1sMax.push_back(alignedToRef1);
		}
	}
	
	// Return nothing if we did not find a split alignment
	if (maxScore == 0)
	{
		return;
	}
	
	// Create list of start cells with same maximum score
	CellVec startCells1;
	CellVec startCells2;
	bool done = false;
	for (IntegerVecConstIter alignedToRef1Iter = alignedToRef1sMax.begin(); alignedToRef1Iter != alignedToRef1sMax.end(); alignedToRef1Iter++)
	{
		int alignedToRef1 = *alignedToRef1Iter;
		int alignedToRef2 = mRead2.length() - alignedToRef1;
		
		int j1 = alignedToRef1;
		int j2 = alignedToRef2;
		
		int maxAcrossRow1 = 0;
		IntegerVec columns1;
		FindMaxRowEntry(mMatrix1, alignedToRef1, mMinSplitScore, maxAcrossRow1, columns1);
		
		int maxAcrossRow2 = 0;
		IntegerVec columns2;
		FindMaxRowEntry(mMatrix2, alignedToRef2, mMinSplitScore, maxAcrossRow2, columns2);
		
		for (IntegerVecConstIter column1Iter = columns1.begin(); column1Iter != columns1.end(); column1Iter++)
		{
			int i1 = *column1Iter;
			
			for (IntegerVecConstIter column2Iter = columns2.begin(); column2Iter != columns2.end(); column2Iter++)
			{
				int i2 = *column2Iter;

				startCells1.push_back(Cell(i1,j1));
				startCells2.push_back(Cell(i2,j2));
				
				done = firstOnly;

				if (done) { break; }
			}
			
			if (done) { break; }
		}
		
		if (done) { break; }
	}
	
	// Backtrace each split alignment
	for (int startIndex = 0; startIndex < startCells1.size(); startIndex++)
	{
		const Cell& startCell1 = startCells1[startIndex];
		const Cell& startCell2 = startCells2[startIndex];
		
		SplitReadAlignment splitReadAlignment;
		
		splitReadAlignment.refSplit.first = startCell1.i;
		splitReadAlignment.refSplit.second = mReference2.length() - startCell2.i - 1;
		splitReadAlignment.readSplit.first = startCell1.j;
		splitReadAlignment.readSplit.second = startCell2.j;
		splitReadAlignment.score = maxScore;
		splitReadAlignment.score1 = mMatrix1(startCell1);
		splitReadAlignment.score2 = mMatrix2(startCell2);
		
		DebugCheck(splitReadAlignment.score == splitReadAlignment.score1 + splitReadAlignment.score2);
		
		if (backTrace)
		{
			BackTracePath(mBackTrace1, startCell1, splitReadAlignment.matches1);
			BackTracePath(mBackTrace2, startCell2, splitReadAlignment.matches2);
			ReverseMatches(mReference2.length(), mRead2.length(), splitReadAlignment.matches2);
		}
		
		alignments.push_back(splitReadAlignment);
	}
}

