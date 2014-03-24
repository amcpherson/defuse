/*
 *  EstCatalog.cpp
 *  estislands
 *
 *  Created by Andrew McPherson on 09/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "EstCatalog.h"

#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

void EstCatalog::ReadEsts(istream& in)
{
	string line;
	int lineNumber = 0;
	
	while (getline(in, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		if (line[0] == '#' || !isdigit(line[0]))
		{
			continue;
		}
		
		vector<string> estFields;
		split(estFields, line, is_any_of("\t"));
		
		if (estFields.size() < 18)
		{
			continue;
		}
		
		string chromosome = estFields[14];
	
		if (chromosome.substr(0,3) == "chr")
		{
			chromosome = chromosome.substr(3);
		}

		if (chromosome == "M")
		{
			chromosome = "MT";
		}
		
		Segment alignment;
		alignment.start = lexical_cast<int>(estFields[16]) + 1;
		alignment.end = lexical_cast<int>(estFields[17]);
		
		mChrSegments[chromosome].push_back(alignment);
	}
	
	for (ChrSegVecMapIter chrSegmentsIter = mChrSegments.begin(); chrSegmentsIter != mChrSegments.end(); chrSegmentsIter++)
	{
		SegmentVec& segments = chrSegmentsIter->second;
		SortAndMergeSegments(segments,segments);
	}
}

void EstCatalog::SortAndMergeSegments(SegmentVec& inSegments, SegmentVec& outSegments)
{
	if (inSegments.size() == 0)
	{
		return;
	}
	
	SegmentVec merged;
	
	sort(inSegments.begin(), inSegments.end());
	
	Segment currentSegment = inSegments.front();
	for (SegmentVecConstIter segmentIter = inSegments.begin(); segmentIter != inSegments.end(); segmentIter++)
	{
		if (segmentIter->start > currentSegment.end)
		{
			merged.push_back(currentSegment);
			
			currentSegment = *segmentIter;
		}
		else
		{
			currentSegment.end = max(currentSegment.end, segmentIter->end);
		}
	}
	
	merged.push_back(currentSegment);
	
	swap(merged,outSegments);
}

void EstCatalog::FilterContainedInEstIslands(istream& in, ostream& out)
{
	string line;
	
	while (getline(in, line))
	{
		if (line.length() == 0)
		{
			continue;
		}
		
		if (line[0] == '#' || !isdigit(line[0]))
		{
			continue;
		}
		
		vector<string> alignmentFields;
		split(alignmentFields, line, is_any_of("\t"));
		
		if (alignmentFields.size() < 18)
		{
			continue;
		}
		
		string chromosome = alignmentFields[13];
		
		if (chromosome.substr(0,3) == "chr")
		{
			chromosome = chromosome.substr(3);
		}

		if (chromosome == "M")
		{
			chromosome = "MT";
		}

		Segment searchSegment;
		searchSegment.start = lexical_cast<int>(alignmentFields[15]) + 1;
		searchSegment.end = lexical_cast<int>(alignmentFields[16]);
		
		const SegmentVec& segments = mChrSegments[chromosome];
		
		SegmentVecConstIter estIslandIter = lower_bound(segments.begin(), segments.end(), searchSegment);

		if (estIslandIter != segments.begin())
		{
			estIslandIter--;
		}
		
		bool contained = false;
		while (estIslandIter != segments.end() && estIslandIter->start <= searchSegment.end)
		{
			Segment estIsland = *estIslandIter;
			estIsland.start -= 300;
			estIsland.end += 300;

			if (estIsland.Contains(searchSegment))
			{
				contained = true;
			}

			estIslandIter++;
		}
				
		if (contained)
		{
			out << line << endl;
		}
	}
}

