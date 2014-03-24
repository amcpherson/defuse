/*
 *  ReadIndex.h
 *
 *  Created by Andrew McPherson on 02/09/10.
 *
 */

#ifndef READINDEX_H_
#define READINDEX_H_

#include "Common.h"

#include <string>
#include <fstream>

using namespace std;

class ReadIndex
{
public:
	void Open(const string& readsPrefix);
	void Find(int fragmentIndex, string& readSeq1, string& readSeq2) const;
	void Find(int fragmentIndex, int readEnd, string& readSeq) const;
	void Close();
	
private:
	void CheckIndexFail(int fragmentIndex) const;
	void CheckFastqFail(int fragmentIndex, int readEnd) const;

	string mReadsIndexFilename;
	string mReadsFastqFilename[2];

	mutable ifstream mReadsIndexFile;
	mutable ifstream mReadsFastqFile[2];
};

#endif
