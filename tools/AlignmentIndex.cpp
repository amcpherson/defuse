/*
 *  AlignmentIndex.cpp
 *
 *  Created by Andrew McPherson on 02/09/10.
 *
 */

#include "AlignmentIndex.h"
#include "Common.h"
#include "Indexer.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;

AlignmentIndex::AlignmentIndex() : mBamFile(0), mBamIndexFile(0) {}

void AlignmentIndex::Open(const string& bamFilename)
{
	mBamFile = samopen(bamFilename.c_str(), "rb", 0);
	mBamIndexFile = bam_index_load(bamFilename.c_str());

	if (mBamFile == 0)
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}
	
	if (mBamIndexFile == 0)
	{
		cerr << "Error: Unable to open bam index for file " << bamFilename << endl;
		exit(1);
	}
	
	mReferenceIndexMap.clear();
	for (int targetIndex = 0; targetIndex < mBamFile->header->n_targets; targetIndex++)
	{
		mReferenceIndexMap[mBamFile->header->target_name[targetIndex]] = targetIndex;
	}
}

struct FetchData
{
	FetchData(int strand, int referenceIndex, CompAlignVec& alignments) : strand(strand), referenceIndex(referenceIndex), alignments(alignments)
	{}
	
	int strand;
	int referenceIndex;
	CompAlignVec& alignments;
};

int FetchCallback(const bam1_t* b, void* data)
{
	FetchData* fetchData = (FetchData*)data;

	int alignmentStrand = (b->core.flag & 0x10) ? MinusStrand : PlusStrand;
	
	if (alignmentStrand == fetchData->strand)
	{
		// Try to split qname into id and end
		string qname = string((char*)b->data);
		vector<string> qnameFields;
		split(qnameFields, qname, is_any_of("/"));

		// Fragment index encoded in query name
		int fragmentIndex = lexical_cast<int>(qnameFields[0]);
		
		// Read end encoded in query name or flag
		int readEnd;
		if (b->core.flag & 0x40)
		{
			// Read is end 1
			readEnd = 0;
		}
		else if (b->core.flag & 0x80)
		{
			// Read is end 2
			readEnd = 1;
		}
		else if (qnameFields.size() == 2)
		{
			// Read is in query name
			readEnd = (qnameFields[1] == "1") ? 0 : 1;
		}
		else
		{
			cerr << "Error: Unable to determine read end for qname " << qname << endl;
			exit(1);
		}
		
		CompactAlignment compactAlignment;
		compactAlignment.readID.fragmentIndex = fragmentIndex;
		compactAlignment.readID.readEnd = readEnd;
		compactAlignment.refStrand.referenceIndex = fetchData->referenceIndex;
		compactAlignment.refStrand.strand = alignmentStrand;
		compactAlignment.region.start = b->core.pos + 1;
		compactAlignment.region.end = b->core.pos + b->core.l_qseq;
		
		fetchData->alignments.push_back(compactAlignment);
	}
	
	return 0;
}

void AlignmentIndex::Find(const string& reference, int strand, int start, int end, NameIndex& referenceNames, CompAlignVec& alignments) const
{
	if (mReferenceIndexMap.find(reference) == mReferenceIndexMap.end())
	{
		cerr << "Error: Unable to find alignments to reference " << reference << endl;
		exit(1);
	}
	
	if (start > end)
	{
		cerr << "Error: start " << start << " end " << end << " for alignments query" << endl;
		exit(1);
	}
	
	// Bam indexing is 0 based with end being 1 after the last base
	start = max(0, start - 1);
	
	int referenceIndex = referenceNames.Index(reference);
	
	FetchData fetchData(strand,referenceIndex,alignments);
	
	bam_fetch(mBamFile->x.bam, mBamIndexFile, mReferenceIndexMap.find(reference)->second, start, end, &fetchData, FetchCallback);
}

void AlignmentIndex::Close()
{
	bam_index_destroy(mBamIndexFile);
	samclose(mBamFile);
}
