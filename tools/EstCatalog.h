/*
 *  EstCatalog.h
 *  estislands
 *
 *  Created by Andrew McPherson on 09/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ESTCATALOG_H_
#define ESTCATALOG_H_

#include <vector>
#include <iostream>
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

using namespace std;
using namespace boost;

class EstCatalog
{
public:	
	void ReadEsts(istream& in);
	void FilterContainedInEstIslands(istream& in, ostream& out);	

private:
	struct Segment
	{
		int start;
		int end;

		bool operator<(const Segment& other) const
		{
			return start < other.start;
		}
		
		bool Contains(const Segment& other) const
		{
			return start <= other.start && end >= other.end;
		}

		bool Overlaps(const Segment& other) const
		{
			return !(end < other.start || start > other.end);
		}
	};
		
	typedef vector<Segment> SegmentVec;
	typedef vector<Segment>::iterator SegmentVecIter;
	typedef vector<Segment>::const_iterator SegmentVecConstIter;
	
	typedef unordered_map<string,SegmentVec> ChrSegVecMap;
	typedef unordered_map<string,SegmentVec>::iterator ChrSegVecMapIter;
	typedef unordered_map<string,SegmentVec>::const_iterator ChrSegVecMapConstIter;
		
	void SortAndMergeSegments(SegmentVec& inSegments, SegmentVec& outSegments);
	bool Overlapping(const Segment& a, const Segment& b);
	
	ChrSegVecMap mChrSegments;
};

#endif
