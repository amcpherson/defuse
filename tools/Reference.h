/*
 *  Reference.h
 *  transcript_to_sequence
 *
 *  Created by Andrew McPherson on 02/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include "Common.h"

#include <map>
#include <string>
#include <iostream>

using namespace std;

class Reference
{
public:
	Reference(const string& fastaFilename);
	
	void GetSequences(const LocationVec& locations, vector<string>& sequences);
	
private:
	const string mFastaFilename;
};

#endif
