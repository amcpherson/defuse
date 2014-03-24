/*
 *  FragmentAlignmentsIterator.h
 *
 *  Created by Andrew McPherson on 03/18/11.
 *
 */

#ifndef FRAGMENTALIGNMENSITERATOR_H_
#define FRAGMENTALIGNMENSITERATOR_H_

#include "Common.h"
#include "bam.h"
#include "sam.h"

#include <map>
#include <string>
#include <iostream>

using namespace std;


class FragmentAlignmentsIterator
{
public:
	FragmentAlignmentsIterator();
	~FragmentAlignmentsIterator();
	
	void Open(const string& bamFilename, StringVec& referenceNames);
	bool GetNext(CompAlignVec& alignments);
	void Close();
	
private:
	void ConvertBamAlignment(const bam1_t& bamEntry, CompactAlignment& alignment);
	
	samfile_t* mBamFile;
	bool mGood;
	bam1_t mCurrentEntry;
};

#endif
