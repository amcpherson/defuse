/*
 *  AlignmentStream.cpp
 *  findbreaks
 *
 *  Created by Andrew McPherson on 12/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "AlignmentStream.h"

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;

IAlignmentStream* IAlignmentStream::Create(const string& alignmentFilename)
{
	string::size_type extStart = alignmentFilename.find_last_of('.');
	string alignmentExt = alignmentFilename.substr(extStart + 1);

	IAlignmentStream* alignmentStream = 0;
	if (alignmentExt == "bwtout")
	{
		alignmentStream = new BowtieAlignmentStream(alignmentFilename);
	}
	else if (alignmentExt == "novo")
	{
		alignmentStream = new NovoAlignmentStream(alignmentFilename);
	}
	else
	{
		cerr << "Error: unrecognized extension " << alignmentExt << endl;
		return 0;
	}
	
	if (alignmentStream->Good())
	{
		return alignmentStream;
	}
	else
	{
		cerr << "Error: unable to open file " << alignmentFilename << endl;
		delete alignmentStream;
		return 0;
	}
}

IAlignmentStream::IAlignmentStream(const string& alignmentFilename) : mStream(alignmentFilename.c_str()), mLineNumber(0)
{
}

bool IAlignmentStream::Good()
{
	return mStream.good();
}

void IAlignmentStream::Reset()
{
	mStream.seekg(0);
}

BowtieAlignmentStream::BowtieAlignmentStream(const string& filename) : IAlignmentStream(filename)
{
}

bool BowtieAlignmentStream::GetNextAlignment(RawAlignment& alignment)
{
	string line;
	
	while (getline(mStream, line))
	{
		mLineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> alignmentFields;
		split(alignmentFields, line, is_any_of("\t"));
		
		if (alignmentFields.size() < 7)
		{
			continue;
		}
		
		alignmentFields.push_back("");
		
		string readName = alignmentFields[0];
		string strandName = alignmentFields[1];
		string referenceName = alignmentFields[2];
		string startPosition = alignmentFields[3];
		string sequence = alignmentFields[4];
		string mismatches = alignmentFields[7];
		
		// Count mismatches
		int numMismatches = count(mismatches.begin(),mismatches.end(),':');
		
		// Interpret strand
		int strand;
		if (strandName == "+")
		{
			strand = PlusStrand;
		}
		else if (strandName == "-")
		{
			strand = MinusStrand;
		}
		else
		{
			cerr << "Error: Unable to interpret strand" << endl;
			return false;
		}
		
		// Reverse complement sequence if aligned to minus strand
		if (strand == MinusStrand)
		{
			ReverseComplement(sequence);
		}
		
		// Split read name into id and end
		string::size_type readEndStart = readName.find_first_of('/');
		if(readEndStart == std::string::npos && readEndStart + 1 >= readName.length())
		{
			cerr << "Error: Unable to interpret read name" << endl;
			return false;
		}
		
		string fragmentName = readName.substr(0, readEndStart);
		char readEndName = readName.substr(readEndStart + 1, 1)[0];
		
		// Check read end
		if (readEndName != '1' && readEndName != '2')
		{
			cerr << "Error: Unable to interpret read end" << endl;
			return false;
		}
	
		int readEnd = (readEndName == '1') ? 0 : 1;
		
		alignment.fragment = fragmentName;
		alignment.readEnd = readEnd;
		alignment.reference = referenceName;
		alignment.strand = strand;
		alignment.region.start = lexical_cast<int>(startPosition) + 1;
		alignment.region.end = alignment.region.start + sequence.length() - 1;
		alignment.sequence = sequence;
		alignment.numMis = numMismatches;

		return true;
	}
	
	return false;	
}

NovoAlignmentStream::NovoAlignmentStream(const string& filename) : IAlignmentStream(filename)
{
}

bool NovoAlignmentStream::GetNextAlignment(RawAlignment& alignment)
{
	string line;
	
	while (getline(mStream, line))
	{
		mLineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		if (line[0] == '#')
		{
			continue;
		}
		
		vector<string> alignmentFields;
		split(alignmentFields, line, is_any_of("\t"));
		
		if (alignmentFields.size() < 13)
		{
			continue;
		}
		
		alignmentFields.push_back("");
		
		string readName = alignmentFields[0];
		string sequence = alignmentFields[2];
		string status = alignmentFields[4];
		string referenceName = alignmentFields[7];
		string startPosition = alignmentFields[8];
		string strandName = alignmentFields[9];
		string mismatches = alignmentFields[13];
		
		// Split read name into id and end
		string::size_type fragmentNameStart = readName.find_first_of('@');
		if(fragmentNameStart == std::string::npos)
		{
			fragmentNameStart = 0;
		}
		else
		{
			fragmentNameStart++;
		}
		
		string::size_type readEndStart = readName.find_first_of('/');
		if(readEndStart == std::string::npos && readEndStart + 1 >= readName.length())
		{
			cerr << "Error: Unable to interpret read name" << endl;
			return false;
		}
		
		string fragmentName = readName.substr(fragmentNameStart, readEndStart - fragmentNameStart);
		char readEndName = readName.substr(readEndStart + 1, 1)[0];
		
		// Check read end
		if (readEndName != '1' && readEndName != '2')
		{
			cerr << "Error: Unable to interpret read end" << endl;
			return false;
		}
		
		// Filter poor quality and non mapping alignments
		if (status == "NM" or status == "QC")
		{
			continue;
		}
		
		// Interpret strand
		int strand;
		if (strandName == "F")
		{
			strand = PlusStrand;
		}
		else if (strandName == "R")
		{
			strand = MinusStrand;
		}
		else
		{
			cerr << "Error: Unable to interpret strand" << endl;
			return false;
		}
		
		// Calculate position
		int numInsertions = count(mismatches.begin(), mismatches.end(), '+');
		int numDeletions = count(mismatches.begin(), mismatches.end(), '-');
		int matchLength = sequence.length() - numInsertions + numDeletions;

		int readEnd = (readEndName == '1') ? 0 : 1;

		alignment.fragment = fragmentName;
		alignment.readEnd = readEnd;
		alignment.reference = referenceName;
		alignment.strand = strand;
		alignment.region.start = lexical_cast<int>(startPosition);
		alignment.region.end = alignment.region.start + matchLength - 1;
		alignment.sequence = sequence;
		alignment.numMis = mismatches.size();
		
		return true;
	}
	
	return false;
}

