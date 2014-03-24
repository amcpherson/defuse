/*
 *  AlignmentStream.h
 *
 */

#ifndef ALIGNMENTSTREAM_H_
#define ALIGNMENTSTREAM_H_

#include "Common.h"

#include <iostream>
#include <fstream>

using namespace std;


struct RawAlignment
{
	string fragment;
	int readEnd;
	string reference;
	int strand;
	Region region;
	string sequence;
};

typedef vector<RawAlignment> RawAlignmentVec;


class AlignmentStream
{
public:
	virtual ~AlignmentStream() {}
	
	virtual bool GetNextAlignment(RawAlignment& alignment) = 0;
};


class SamAlignmentStream : public AlignmentStream
{
public:
	SamAlignmentStream(const string& samFilename);
	~SamAlignmentStream();
	
	bool GetNextAlignment(RawAlignment& alignment);
	
protected:
	istream* mStream;
	int mLineNumber;
};


class CompactAlignmentStream : public AlignmentStream
{
public:
	CompactAlignmentStream(const string& alignFilename);
	~CompactAlignmentStream();
	
	bool GetNextAlignment(RawAlignment& alignment);
	
protected:
	istream* mStream;
	int mLineNumber;
};


class FragmentAlignmentStream
{
public:
	FragmentAlignmentStream(AlignmentStream* alignmentStream);
	
	bool GetNextAlignments(RawAlignmentVec& alignments);
	
protected:
	AlignmentStream* mAlignmentStream;
	RawAlignment mNextAlignment;
	bool mGood;
};


#endif
