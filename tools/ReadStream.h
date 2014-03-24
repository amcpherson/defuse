/*
 *  ReadStream.h
 *  findbreaks
 *
 *  Created by Andrew McPherson on 12/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef READSTREAM_H_
#define READSTREAM_H_

#include "Common.h"

#include <iostream>
#include <fstream>

using namespace std;

struct RawRead
{
	string fragment;
	int readEnd;
	string sequence;
	string quality;
};

class IReadStream
{
public:	
	virtual ~IReadStream() {};
	
	virtual bool Good();
	virtual bool GetNextRead(RawRead& read) = 0;
	
	static IReadStream* Create(const string& filename);

protected:
	IReadStream(const string& filename);
	
	ifstream mStream;
};

class FastqReadStream : public IReadStream
{
public:
	explicit FastqReadStream(const string& filename);	
	virtual bool GetNextRead(RawRead& read);
};

#endif
