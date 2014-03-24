/*
 *  AlignmentStream.h
 *  findbreaks
 *
 *  Created by Andrew McPherson on 12/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
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
	int numMis;
};

class IAlignmentStream
{
public:	
	virtual ~IAlignmentStream() {};
	
	virtual bool Good();
	virtual bool GetNextAlignment(RawAlignment& alignment) = 0;
	virtual void Reset();
	
	static IAlignmentStream* Create(const string& filename);
	
protected:
	IAlignmentStream(const string& filename);

	ifstream mStream;
	int mLineNumber;
};

class BowtieAlignmentStream : public IAlignmentStream
{
public:
	explicit BowtieAlignmentStream(const string& filename);
	virtual bool GetNextAlignment(RawAlignment& alignment);
};

class NovoAlignmentStream : public IAlignmentStream
{
public:
	explicit NovoAlignmentStream(const string& filename);
	virtual bool GetNextAlignment(RawAlignment& alignment);
};

#endif
