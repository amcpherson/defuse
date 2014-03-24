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

void AlignmentIndex::Open(const string& bamFilename)
{
	if (!mReader.Open(bamFilename))
	{
		cerr << "Error: Unable to open bam file " << bamFilename << endl;
		exit(1);
	}
	
	if (!mReader.LocateIndex())
	{
		cerr << "Error: Unable to open bam index for file " << bamFilename << endl;
		exit(1);
	}
}

void AlignmentIndex::Find(const string& reference, int strand, int start, int end, NameIndex& referenceNames, CompAlignVec& alignments) const
{
	if (start > end)
	{
		cerr << "Error: start " << start << " end " << end << " for alignments query" << endl;
		exit(1);
	}
	
	// Bam indexing is 0 based with end being 1 after the last base
	start = max(0, start - 1);
	
	int referenceIndex = mReader.GetReferenceID(reference);
	
	mReader.SetRegion(referenceIndex, start, referenceIndex, end);
	
	BamAlignment alignment;
	while (mReader.GetNextAlignment(alignment))
	{
		if ((strand == MinusStrand) != alignment.IsReverseStrand())
		{
			continue;
		}
		
		// Try to split qname into id and end
		vector<string> qnameFields;
		split(qnameFields, alignment.Name, is_any_of("/"));
		
		// Fragment index encoded in query name
		int fragmentIndex;
		try
		{
			fragmentIndex = lexical_cast<int>(qnameFields[0]);
		}
		catch (boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::bad_lexical_cast> > e)
		{
			cout << "Failed to interpret fragment:" << endl << qnameFields[0] << endl;
			exit(1);
		}
		
		// Read end encoded in query name or flag
		int readEnd;
		if (alignment.IsFirstMate())
		{
			// Read is end 1
			readEnd = 0;
		}
		else if (alignment.IsSecondMate())
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
			cerr << "Error: Unable to determine read end for qname " << alignment.Name << endl;
			exit(1);
		}
		
		CompactAlignment compactAlignment;
		compactAlignment.readID.fragmentIndex = fragmentIndex;
		compactAlignment.readID.readEnd = readEnd;
		compactAlignment.refStrand.referenceIndex = referenceIndex;
		compactAlignment.refStrand.strand = strand;
		compactAlignment.region.start = alignment.Position + 1;
		compactAlignment.region.end = alignment.GetEndPosition();
		
		alignments.push_back(compactAlignment);
	}
}

void AlignmentIndex::Close()
{
	mReader.Close();
}
