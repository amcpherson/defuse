/*
 *  Parameters.cpp
 *  findbreaks
 *
 *  Created by Andrew McPherson on 08/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Parameters.h"

Parameters::Parameters(int minFragmentSize, 
					   int maxFragmentSize, 
					   int maxReadLength, 
					   int minReadAlignment,
					   int minClusterSize,
					   const string& outputDirectory)
: 
mMinFragmentSize(minFragmentSize), 
mMaxFragmentSize(maxFragmentSize), 
mMaxReadLength(maxReadLength), 
mMinReadAlignment(minReadAlignment),
mMinClusterSize(minClusterSize), 
mOutputDirectory(outputDirectory)
{}

int Parameters::GetFragmentSizeMean() const
{
	return (mMaxFragmentSize + mMinFragmentSize) / 2;
}

int Parameters::GetMinFragmentSize() const
{
	return mMinFragmentSize;
}

int Parameters::GetMaxFragmentSize() const
{
	return mMaxFragmentSize;
}

int Parameters::GetMaxFragmentSizeDiff() const
{
	return mMaxFragmentSize - mMinFragmentSize + 1;
}

int Parameters::GetMaxReadLength() const
{
	return mMaxReadLength;
}

int Parameters::GetMinReadAlignment() const
{
	return mMinReadAlignment;
}

int Parameters::GetMaxBreakRegionSize() const
{
	return mMaxFragmentSize - 2 * mMinReadAlignment;
}

int Parameters::GetMinClusterSize() const
{
	return mMinClusterSize;
}

const string& Parameters::GetOutputDirectory() const
{
	return mOutputDirectory;
}

