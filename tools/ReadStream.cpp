/*
 *  AlignmentStream.cpp
 *  findbreaks
 *
 *  Created by Andrew McPherson on 12/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "ReadStream.h"

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;

IReadStream* IReadStream::Create(const string& readsFilename)
{
	string::size_type extStart = readsFilename.find_last_of('.');
	string readsExt = readsFilename.substr(extStart + 1);

	IReadStream* readStream = 0;
	if (readsExt == "fastq" || readsExt == "fq")
	{
		readStream = new FastqReadStream(readsFilename);
	}
	else
	{
		cerr << "Error: unrecognized extension " << readsExt << endl;
		return 0;
	}
	
	if (readStream->Good())
	{
		return readStream;
	}
	else
	{
		cerr << "Error: unable to open file " << readsFilename << endl;
		delete readStream;
		return 0;
	}
}

IReadStream::IReadStream(const string& filename) : mStream(filename.c_str())
{
}

bool IReadStream::Good()
{
	return mStream.good();
}

FastqReadStream::FastqReadStream(const string& filename) : IReadStream(filename)
{
}

bool FastqReadStream::GetNextRead(RawRead& read)
{
	string line[4];
	int lineIndex = 0;
	
	while (lineIndex < 4 && getline(mStream, line[lineIndex]))
	{
		lineIndex++;
	}
	
	if (lineIndex < 4)
	{
		return false;
	}
	
	// Split read name into id and end
	if (line[0][0] != '@')
	{
		cerr << "Error: Unable to interpret read name " << line[0] << endl;
		return false;
	}

	string::size_type readEndStart = line[0].find_first_of('/');
	if(readEndStart == std::string::npos && readEndStart + 1 >= line[0].length())
	{
		cerr << "Error: Unable to interpret read name " << line[0] << endl;
		return false;
	}
	
	char readEndName = line[0].substr(readEndStart + 1, 1)[0];
	if(readEndName != '1' && readEndName != '2')
	{
		cerr << "Error: Unable to interpret read end " << line[0] << endl;
		return false;
	}
	
	int readEnd = (readEndName == '1') ? 0 : 1;
	
	read.fragment = line[0].substr(1, readEndStart - 1);
	read.readEnd = readEnd;
	read.sequence = line[1];
	read.quality = line[3];
	
	return true;
}

