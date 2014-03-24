/*
 *  Reference.cpp
 *  transcript_to_sequence
 *
 *  Created by Andrew McPherson on 02/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Alignments.h"
#include "Reference.h"
#include "DebugCheck.h"

#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

Reference::Reference(const string& fastaFilename) : mFastaFilename(fastaFilename)
{}

void Reference::GetSequences(const LocationVec& locations, vector<string>& sequences)
{
	typedef multimap<string,int> RefIndexMMap;
	typedef multimap<string,int>::const_iterator RefIndexMMapIter;
	
	RefIndexMMap refToIndex;
	
	sequences.resize(locations.size());
	
	for (int locationIndex = 0; locationIndex < locations.size(); locationIndex++)
	{
		const Location& location = locations[locationIndex];
		
		const string& referenceName = Alignments::ReferenceNames.Get(location.refStrand.referenceIndex);
		
		refToIndex.insert(make_pair(referenceName,locationIndex));
	}

	ifstream referenceFile(mFastaFilename.c_str());
	if (!referenceFile.good())
	{
		cerr << "Error: unable to open file " << mFastaFilename << endl;
		exit(1);
	}
	
	string referenceName;
	string sequence;
	
	string line;
	while (getline(referenceFile, line))
	{
		if (line.size() == 0)
		{
			continue;
		}
		
		if (line[0] == '>')
		{
			pair<RefIndexMMapIter,RefIndexMMapIter> indexIterRange = refToIndex.equal_range(referenceName);
			for (RefIndexMMapIter indexIter = indexIterRange.first; indexIter != indexIterRange.second; indexIter++)
			{
				int index = indexIter->second;
				const Location& location = locations[index];
				
				if (location.region.end < 1 || location.region.start > (int)sequence.size())
				{
					continue;
				}
				
				int seqStart = max(location.region.start, 1);
				int numPrefixed = max(1 - location.region.start, 0);
				
				int seqEnd = min(location.region.end, (int)sequence.size());
				int numSuffixed = max(location.region.end - (int)sequence.size(), 0);
				
				sequences[index] = string(numPrefixed, 'N');
				sequences[index] += sequence.substr(seqStart - 1, seqEnd - seqStart + 1);
				sequences[index] += string(numSuffixed, 'N');
				
				transform(sequences[index].begin(), sequences[index].end(), sequences[index].begin(), ::toupper);

				DebugCheck(sequences[index].size() == location.region.end - location.region.start + 1);
			}
			
			referenceName = line.substr(1);
			sequence = "";
		}
		else
		{
			sequence.append(line);
		}
	}	

	pair<RefIndexMMapIter,RefIndexMMapIter> indexIterRange = refToIndex.equal_range(referenceName);
	for (RefIndexMMapIter indexIter = indexIterRange.first; indexIter != indexIterRange.second; indexIter++)
	{
		int index = indexIter->second;
		const Location& location = locations[index];
		
		if (location.region.end < 1 || location.region.start > (int)sequence.size())
		{
			continue;
		}
		
		int seqStart = max(location.region.start, 1);
		int numPrefixed = max(1 - location.region.start, 0);
		
		int seqEnd = min(location.region.end, (int)sequence.size());
		int numSuffixed = max(location.region.end - (int)sequence.size(), 0);
		
		sequences[index] = string(numPrefixed, 'N');
		sequences[index] += sequence.substr(seqStart - 1, seqEnd - seqStart + 1);
		sequences[index] += string(numSuffixed, 'N');
		
		transform(sequences[index].begin(), sequences[index].end(), sequences[index].begin(), ::toupper);
		
		DebugCheck(sequences[index].size() == location.region.end - location.region.start + 1);
	}
	
	referenceFile.close();
}

