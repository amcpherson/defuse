/*
 *  Parameters.h
 *  findbreaks
 *
 *  Created by Andrew McPherson on 08/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>

using namespace std;

class Parameters
{
public:
	Parameters(int minFragmentSize,
			   int maxFragmentSize,
			   int maxReadLength,
			   int minReadAlignment,
			   int minClusterSize,
			   const string& outputDirectory);
		
	int GetFragmentSizeMean() const;
	int GetMinFragmentSize() const;
	int GetMaxFragmentSize() const;
	int GetMaxFragmentSizeDiff() const;
	int GetMaxReadLength() const;
	int GetMinReadAlignment() const;
	int GetMaxBreakRegionSize() const;
	int GetMinClusterSize() const;
	const string& GetOutputDirectory() const;

private:	
	int mMinFragmentSize;
	int mMaxFragmentSize;
	int mMaxReadLength;
	int mMinReadAlignment;
	int mMinClusterSize;
	string mOutputDirectory;
};

#endif

