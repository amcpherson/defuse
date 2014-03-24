/*
 *  FastaIndex.h
 *
 *  Created by Andrew McPherson on 02/09/10.
 *
 */

#ifndef FASTAINDEX_H_
#define FASTAINDEX_H_

#include "Common.h"
#include "faidx.h"

#include <string>

using namespace std;

class FastaIndex
{
public:
	FastaIndex();
	
	void Open(const string& fastaFilename);
	void Get(const string& reference, int strand, int& start, int& length, string& sequence) const;
	void Close();
	
private:
	faidx_t* mFaidx;
};

#endif
