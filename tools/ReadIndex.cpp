/*
 *  ReadIndex.cpp
 *
 *  Created by Andrew McPherson on 02/09/10.
 *
 */

#include "ReadIndex.h"
#include "Common.h"

#include <string>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

void ReadIndex::Open(const string& readsPrefix)
{
	mReadsIndexFilename = readsPrefix + ".fqi";
	mReadsIndexFile.open(mReadsIndexFilename.c_str(), fstream::in | fstream::binary);
	if (!mReadsIndexFile.good())
	{
		cerr << "Error: Unable to open file " << mReadsIndexFilename << endl;
		exit(1);
	}
	
	mReadsFastqFilename[0] = readsPrefix + ".1.fastq";
	mReadsFastqFilename[1] = readsPrefix + ".2.fastq";
	
	for (int readEnd = 0; readEnd <= 1; readEnd++)
	{
		mReadsFastqFile[readEnd].open(mReadsFastqFilename[readEnd].c_str(), fstream::in | fstream::binary);
		if (!mReadsFastqFile[readEnd].good())
		{
			cerr << "Error: Unable to open file " << mReadsFastqFilename[readEnd] << endl;
			exit(1);
		}
	}
}

void ReadIndex::CheckIndexFail(int fragmentIndex) const
{
	if (mReadsIndexFile.fail())
	{
		cerr << "Error: Failure reading index file " << mReadsIndexFilename << " when searching for fragment " << fragmentIndex << endl;
		exit(1);
	}
}

void ReadIndex::CheckFastqFail(int fragmentIndex, int readEnd) const
{
	if (mReadsFastqFile[readEnd].fail())
	{
		cerr << "Error: Failure reading fastq file " << mReadsFastqFilename[readEnd] << " when searching for fragment " << fragmentIndex << " end " << readEnd << endl;
		exit(1);
	}
}

void ReadIndex::Find(int fragmentIndex, string& readSeq1, string& readSeq2) const
{
	Find(fragmentIndex, 0, readSeq1);
	Find(fragmentIndex, 1, readSeq2);
}

void ReadIndex::Find(int fragmentIndex, int readEnd, string& readSeq) const
{
	long int fragmentFilePos = fragmentIndex * 2 * sizeof(long int) + readEnd * sizeof(long int);
	mReadsIndexFile.seekg(fragmentFilePos);
	CheckIndexFail(fragmentIndex);
	
	long int fastqFilePos;
	mReadsIndexFile.read((char*)&fastqFilePos, sizeof(long int));
	CheckIndexFail(fragmentIndex);
	
	mReadsFastqFile[readEnd].seekg(fastqFilePos);
	CheckFastqFail(fragmentIndex, readEnd);

	string line[4];
	for (int lineIndex = 0; lineIndex < 4; lineIndex++)
	{
		getline(mReadsFastqFile[readEnd], line[lineIndex]);
		CheckFastqFail(fragmentIndex, readEnd);
	}
	
	// Check the read name begins with @
	if (line[0][0] != '@')
	{
		cerr << "Error: Unable to interpret read name " << line[0] << " when searching for fragment " << fragmentIndex << " end " << readEnd << endl;
		exit(1);
	}

	// Split read name into id and end
	string::size_type readEndStart = line[0].find_first_of('/');
	if(readEndStart == std::string::npos && readEndStart + 1 >= line[0].length())
	{
		cerr << "Error: Unable to interpret read name " << line[0] << " when searching for fragment " << fragmentIndex << " end " << readEnd << endl;
		exit(1);
	}
	
	// Read and check read end name
	char readEndName = line[0].substr(readEndStart + 1, 1)[0];
	if(readEndName != '1' && readEndName != '2')
	{
		cerr << "Error: Unable to interpret read end " << line[0] << " when searching for fragment " << fragmentIndex << " end " << readEnd << endl;
		exit(1);
	}
	
	int fastqReadEnd = (readEndName == '1') ? 0 : 1;
	string fastqFragmentName = line[0].substr(1, readEndStart - 1);
	int fastqFragmentIndex = lexical_cast<int>(fastqFragmentName);
	
	// Check fragment index read from the fastq is the same
	if (fastqFragmentIndex != fragmentIndex)
	{
		cerr << "Error: Fragment index mismatch when interpreting " << line[0] << " and searching for fragment " << fragmentIndex << " end " << readEnd << endl;
		exit(1);
	}
	
	// Check read end read from the fastq is the same
	if (fastqReadEnd != readEnd)
	{
		cerr << "Error: Read end mismatch when interpreting " << line[0] << " and searching for fragment " << fragmentIndex << " end " << readEnd << endl;
		exit(1);
	}

	readSeq = line[1];
}

void ReadIndex::Close()
{
	mReadsIndexFile.close();
	mReadsFastqFile[0].close();
	mReadsFastqFile[1].close();
}
