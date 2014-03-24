/*
 *  Matrix.h
 *  findbreaks
 *
 *  Created by Andrew McPherson on 30/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "Common.h"

#include <vector>

using namespace std;

struct Cell
{
	Cell() : i(0), j(0) {}
	Cell(int i, int j) : i(i), j(j) {}
	int i;
	int j;

	inline bool operator==(const Cell& other)
	{
		return i == other.i && j == other.j;
	}
	
	inline bool operator<(const Cell& other)
	{
		if (i < other.i)
		{
			return true;
		}
		else if (i > other.i)
		{
			return false;
		}
		
		if (j < other.j)
		{
			return true;
		}
		
		return false;
	}
};

template <typename TType>
class Matrix
{
public:
	Matrix() : mLength(0), mHeight(0)
	{}
	
	Matrix(int length, int height) : mLength(length), mHeight(height)
	{
		mMatrix.resize(length * height);
	}
	
	inline TType& operator()(int i, int j)
	{
		return mMatrix[j * mLength + i];		
	}
	
	inline const TType& operator()(int i, int j) const
	{
		return mMatrix[j * mLength + i];		
	}
	
	inline TType& operator()(const Cell& cell)
	{
		return mMatrix[cell.j * mLength + cell.i];
	}
	
	inline const TType& operator()(const Cell& cell) const
	{
		return mMatrix[cell.j * mLength + cell.i];
	}
	
	inline int GetLength() const
	{
		return mLength;
	}

	inline int GetHeight() const
	{
		return mHeight;
	}

	void Clear(const TType& clearValue)
	{
		fill(mMatrix.begin(), mMatrix.end(), clearValue);
	}
	
	void Resize(int length, int height)
	{
		if (length != mLength || height != mHeight)
		{
			mLength = length;
			mHeight = height;
			
			mMatrix.resize(length * height);
		}
	}
	
private:
	int mLength;
	int mHeight;
	vector<TType> mMatrix;
};

#endif

