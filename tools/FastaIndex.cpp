/*
 *  FastaIndex.cpp
 *
 *  Created by Andrew McPherson on 02/09/10.
 *
 */

#include "FastaIndex.h"

#include <sstream>
#include <iostream>

using namespace std;

FastaIndex::FastaIndex() : mFaidx(0) {}

void FastaIndex::Open(const string& fastaFilename)
{
	mFaidx = fai_load(fastaFilename.c_str());	
}

void FastaIndex::Get(const string& reference, int strand, int& start, int& length, string& sequence) const
{
	stringstream regionStr;
	
	if (length < 0)
	{
		sequence = "";
		return;
	}
	
	if (start < 1)
	{
		length -= 1 - start;
		start = 1;
	}
	
	int end = start + length - 1;
	
	regionStr << reference << ":" << start << "-" << end;
	
	char* seq = fai_fetch(mFaidx, regionStr.str().c_str(), &length);
	
	if (seq == 0)
	{
		cerr << "Error: Unable to find sequence for " << reference << endl;
		exit(1);
	}
	
	sequence = string(seq);
	
	free(seq);
	
	if (strand == MinusStrand)
	{
		ReverseComplement(sequence);
	}
}

void FastaIndex::Close()
{
	fai_destroy(mFaidx);
}

