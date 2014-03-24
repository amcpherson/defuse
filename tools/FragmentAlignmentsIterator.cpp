/*
 *  FragmentAlignmentsIterator.cpp
 *
 *  Created by Andrew McPherson on 03/18/11.
 *
 */

#include "FragmentAlignmentsIterator.h"
#include "Common.h"
#include "Indexer.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

FragmentAlignmentsIterator::FragmentAlignmentsIterator() : mBamFile(0), mGood(false)
{
	mCurrentEntry.data = 0;
}

FragmentAlignmentsIterator::~FragmentAlignmentsIterator()
{
	Close();
	free(mCurrentEntry.data);
}

void FragmentAlignmentsIterator::Open(const string& bamFilename, StringVec& referenceNames)
{
	mBamFile = samopen(bamFilename.c_str(), "rb", 0);
	
	if (mBamFile == 0)
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}
	
	referenceNames.clear();
	for (int targetIndex = 0; targetIndex < mBamFile->header->n_targets; targetIndex++)
	{
		referenceNames.push_back(string(mBamFile->header->target_name[targetIndex]));
	}
	
	memset(&mCurrentEntry, 0, sizeof(bam1_t));
	mGood = bam_read1(mBamFile->x.bam, &mCurrentEntry) > 0;
}

bool FragmentAlignmentsIterator::GetNext(CompAlignVec& alignments)
{
	alignments.clear();

	if (!mGood)
	{
		return false;
	}
	
	CompactAlignment firstAlignment;
	ConvertBamAlignment(mCurrentEntry, firstAlignment);
	
	alignments.push_back(firstAlignment);
	
	while (true)
	{
		mGood = bam_read1(mBamFile->x.bam, &mCurrentEntry) > 0;
		if (!mGood)
		{
			break;
		}
		
		if (mCurrentEntry.core.tid < 0)
		{
			continue;
		}
		
		CompactAlignment nextAlignment;
		ConvertBamAlignment(mCurrentEntry, nextAlignment);
		
		if (nextAlignment.readID.fragmentIndex == firstAlignment.readID.fragmentIndex)
		{
			alignments.push_back(nextAlignment);
		}
		else
		{
			break;
		}
	}
	
	return true;
}

void FragmentAlignmentsIterator::Close()
{
	samclose(mBamFile);
	mBamFile = 0;
}

void FragmentAlignmentsIterator::ConvertBamAlignment(const bam1_t& bamEntry, CompactAlignment& alignment)
{
	// Split qname into id and end
	// Fragment index encoded in fragment name
	string qname = string((char*)bamEntry.data);
	vector<string> qnameFields;
	split(qnameFields, qname, is_any_of("/"));
	
	alignment.readID.fragmentIndex = lexical_cast<int>(qnameFields[0]);
	alignment.readID.readEnd = (qnameFields[1] == "1") ? 0 : 1;
	alignment.refStrand.referenceIndex = bamEntry.core.tid;
	alignment.refStrand.strand = (bamEntry.core.flag & 0x10) ? MinusStrand : PlusStrand;
	alignment.region.start = bamEntry.core.pos + 1;
	alignment.region.end = bamEntry.core.pos + bamEntry.core.l_qseq;
}


