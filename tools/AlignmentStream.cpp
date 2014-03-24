/*
 *  AlignmentStream.cpp
 *
 */

#include "AlignmentStream.h"

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


SamAlignmentStream::SamAlignmentStream(const string& samFilename) : mStream(0), mLineNumber(0)
{
	if (samFilename == "-")
	{
		mStream = &cin;
	}
	else
	{
		mStream = new ifstream(samFilename.c_str());
		
		if (!mStream->good())
		{
			cerr << "Error: Unable to open sam file " << samFilename << endl;
			exit(1);
		}
	}
}

SamAlignmentStream::~SamAlignmentStream()
{
	delete mStream;
}

bool SamAlignmentStream::GetNextAlignment(RawAlignment& alignment)
{
	string line;
	while (getline(*mStream, line))
	{
		mLineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty alignment line " << mLineNumber << endl;
			exit(1);
		}
		
		if (line[0] == '@')
		{
			continue;
		}
		
		vector<string> alignmentFields;
		split(alignmentFields, line, is_any_of("\t"));
		
		if (alignmentFields.size() < 10)
		{
			cerr << "Error: Format error for alignment line " << mLineNumber << endl;
			exit(1);
		}
		
		const string& qname = alignmentFields[0];
		int flag = lexical_cast<int>(alignmentFields[1]);
		const string& rname = alignmentFields[2];
		int pos = lexical_cast<int>(alignmentFields[3]);
		const string& seq = alignmentFields[9];
		
		if (rname == "*")
		{
			continue;
		}
		
		// Interpret strand
		int strand;
		if ((flag & 0x0010) == 0)
		{
			strand = PlusStrand;
		}
		else
		{
			strand = MinusStrand;
		}
		
		string fragment;
		int readEnd;
		
		// Split qname into id and end
		vector<string> qnameFields;
		split(qnameFields, qname, is_any_of("/"));
		if (qnameFields.size() == 2)
		{
			if (qnameFields[1] != "1" && qnameFields[1] != "2")
			{
				cerr << "Error: Unable to interpret qname for alignment line " << mLineNumber << endl;
				exit(1);
			}
			
			fragment = qnameFields[0];
			readEnd = (qnameFields[1] == "1") ? 0 : 1;
		}
		else
		{
			fragment = qname;
			if (flag & 0x0040)
			{
				readEnd = 0;
			}
			else if (flag & 0x0080)
			{
				readEnd = 1;
			}
		}
		
		alignment.fragment = fragment;
		alignment.readEnd = readEnd;
		alignment.reference = rname;
		alignment.strand = strand;
		alignment.region.start = lexical_cast<int>(pos);
		alignment.region.end = alignment.region.start + seq.length() - 1;
		alignment.sequence = seq;
		
		return true;
	}
	
	return false;
}


CompactAlignmentStream::CompactAlignmentStream(const string& alignFilename) : mStream(NULL), mLineNumber(0)
{
	if (alignFilename == "-")
	{
		mStream = &cin;
	}
	else
	{
		mStream = new ifstream(alignFilename.c_str());
		
		if (!mStream->good())
		{
			cerr << "Error: Unable to open alignment file " << alignFilename << endl;
			exit(1);
		}
	}
}

CompactAlignmentStream::~CompactAlignmentStream()
{
	delete mStream;
}

bool CompactAlignmentStream::GetNextAlignment(RawAlignment& alignment)
{
	string line;
	while (getline(*mStream, line))
	{
		mLineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty alignment line " << mLineNumber << endl;
			exit(1);
		}
		
		vector<string> alignmentFields;
		split(alignmentFields, line, is_any_of("\t"));
		
		if (alignmentFields.size() < 6)
		{
			cerr << "Error: Format error for alignment line " << mLineNumber << endl;
			exit(1);
		}
		
		alignment.fragment = alignmentFields[0];
		alignment.readEnd = (alignmentFields[1] == "1") ? 0 : 1;
		alignment.reference = alignmentFields[2];
		alignment.strand = (alignmentFields[3] == "-") ? MinusStrand : PlusStrand;
		alignment.region.start = lexical_cast<int>(alignmentFields[4]);
		alignment.region.end = lexical_cast<int>(alignmentFields[5]);
		
		return true;
	}
	
	return false;	
}



FragmentAlignmentStream::FragmentAlignmentStream(AlignmentStream* alignmentStream) : mAlignmentStream(alignmentStream)
{
	mGood = mAlignmentStream->GetNextAlignment(mNextAlignment);
}

bool FragmentAlignmentStream::GetNextAlignments(RawAlignmentVec& alignments)
{
	if (!mGood)
	{
		return false;
	}
	
	alignments.clear();
	alignments.push_back(mNextAlignment);
	
	while (mGood = mAlignmentStream->GetNextAlignment(mNextAlignment))
	{
		if (alignments.front().fragment != mNextAlignment.fragment)
		{
			break;
		}
		else
		{
			alignments.push_back(mNextAlignment);
		}
	}
	
	return true;
}


