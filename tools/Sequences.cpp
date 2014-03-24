/*
 *  Sequences.cpp
 *
 *  Created by Andrew McPherson on 07/16/11.
 *
 */

#include "Sequences.h"
#include "DebugCheck.h"

#include <fstream>
#include <algorithm>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

void Sequences::Read(const string& fastaFilename)
{
	ifstream fastaFile(fastaFilename.c_str());
	if (!fastaFile.good())
	{
		cerr << "Error: unable to open file " << fastaFilename << endl;
		exit(1);
	}
	
	string id;
	string sequence;
	
	string line;
	while (getline(fastaFile, line))
	{
		if (line.size() == 0)
		{
			continue;
		}
		
		if (line[0] == '>')
		{
			if (!id.empty())
			{
				mSequences[id] = sequence;
			}
			
			id = line.substr(1);
			sequence = "";
		}
		else
		{
			sequence.append(line);
		}
	}
	
	if (!id.empty())
	{
		mSequences[id] = sequence;
	}
}

void Sequences::Get(const string& id, int start, int end, string& sequence) const
{
	if (mSequences.find(id) == mSequences.end())
	{
		cerr << "Error: Unable to find sequence " << id << endl;
		exit(1);
	}
	
	const string& fullSequence = mSequences.find(id)->second;
	
	int seqStart = max(1, start);
	int prependN = seqStart - start;
	
	int seqEnd = min((int)fullSequence.length(), end);
	int appendN = end - seqEnd;
	
	int seqLength = seqEnd - seqStart + 1;
	
	sequence = string(prependN, 'N') + fullSequence.substr(seqStart - 1, seqLength) + string(appendN, 'N');
}

