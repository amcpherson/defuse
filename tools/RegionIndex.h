/*
 *  RegionIndex.h
 *  linkexons
 *
 *  Created by Andrew McPherson on 10/08/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef REGION_INDEX_H
#define REGION_INDEX_H

#include "DebugCheck.h"
#include <map>
#include <set>
#include <vector>
#include <assert.h>
#include <iostream>

using namespace std;

template <class TLocatable>
class RegionIndex
{
public:
	RegionIndex()
	{
		Clear();
	}		

	void Insert(TLocatable* locatable)
	{
		typedef vector<TLocatable*> RegionSet;
		typedef map<int,BreakNode> BreakNodeMap;
		
		const Region& region = locatable->GetRegion();
		
		DebugCheck(region.start <= region.end);
		
		// Create a new leftmost node if necessary
		if (mBreakNodes.size() == 0)
		{
			BreakNode node;
			node.start = region.start;
			node.size = -1;
			mBreakNodes.insert(make_pair(node.start,node));
		}
		else if (region.start < mBreakNodes.begin()->second.start)
		{
			BreakNode node;
			node.start = region.start;
			node.size = mBreakNodes.begin()->second.start - region.start;
			mBreakNodes.insert(make_pair(node.start,node));			
		}
		
		BreakNodeMapIter iterStart = mBreakNodes.upper_bound(region.start);
		iterStart--;
		
		BreakNodeMapIter iterEnd = mBreakNodes.upper_bound(region.end);
		iterEnd--;
		
		// Inserting starts at the start node and ends at the end node
		BreakNodeMapIter iterInsertBegin = iterStart;
		BreakNodeMapIter iterInsertEnd = iterEnd;
		iterInsertEnd++;
		
		// Check if the end node needs to be split
		if (iterEnd->second.start + iterEnd->second.size - 1 != region.end)
		{
			// Create a new node for the part of the end node not covered by this region
			BreakNode brokenEnd = iterEnd->second;
			brokenEnd.start = region.end + 1;
		
			// Adjust the new node's size unless it is the infinite end node
			if (brokenEnd.size >= 0)
			{
				brokenEnd.size = iterEnd->second.size - (region.end - iterEnd->second.start + 1);
			}
		
			// Insert the new node
			pair<BreakNodeMapIter,bool> insResultEnd = mBreakNodes.insert(make_pair(brokenEnd.start,brokenEnd));
			assert(insResultEnd.second);
		
			// Adjust the broken end node's size
			iterEnd->second.size = region.end - iterEnd->second.start + 1;
			
			// Inserting ends one before the inserted node
			iterInsertEnd = insResultEnd.first;
		}

		// Check if the start node needs to be split
		if (iterStart->second.start != region.start)
		{
			// Create a new node for the part of the start node covered by this region
			BreakNode brokenStart = iterStart->second;
			brokenStart.start = region.start;
			brokenStart.size = iterStart->second.size - region.start + iterStart->second.start;
		
			// Adjust the broken start node's size
			iterStart->second.size = region.start - iterStart->second.start;
		
			// Insert the new node
			pair<BreakNodeMapIter,bool> insResultStart = mBreakNodes.insert(make_pair(brokenStart.start,brokenStart));
			assert(insResultStart.second);
			
			// Inserting begins at the new node
			iterInsertBegin = insResultStart.first;
		}
		
		// Add the region to the fully overlapping nodes
		for (BreakNodeMapIter iter = iterInsertBegin; iter != iterInsertEnd; iter++)
		{
			iter->second.regions.insert(locatable);
		}
	}

	bool Delete(TLocatable* locatable)
	{
		const Region& region = locatable->GetRegion();

		BreakNodeMapIter nodeIter = mBreakNodes.upper_bound(region.start);
		nodeIter--;
		
		BreakNodeMapIter nodeIterEnd = mBreakNodes.upper_bound(region.end);
		
		int numErased = 0;
		int numNodes = 0;
		while (nodeIter != nodeIterEnd)
		{
			numErased += nodeIter->second.regions.erase(region);
			numNodes++;
			
			nodeIter++;
		}
		
		return (numErased == numNodes);
	}
	
	void Clear()
	{
		mBreakNodes.clear();
	}
	
	void GetOverlapping(const Region& region, set<TLocatable*>& overlapping) const
	{
		if (mBreakNodes.size() < 1)
		{
			return;
		}
		
		if (region.end < mBreakNodes.begin()->second.start)
		{
			return;
		}
		
		Region searchRegion;
		searchRegion.start = max(region.start, mBreakNodes.begin()->second.start);
		searchRegion.end = region.end;
		
		BreakNodeMapConstIter nodeIter = mBreakNodes.upper_bound(searchRegion.start);
		if (nodeIter == mBreakNodes.begin())
		{
			for (BreakNodeMapConstIter iter1 = mBreakNodes.begin(); iter1 != mBreakNodes.end(); iter1++)
			{
				cerr << iter1->second.start << " " << iter1->second.size << endl;
			}
		}
		nodeIter--;
		
		BreakNodeMapConstIter nodeIterEnd = mBreakNodes.upper_bound(searchRegion.end);
		
		while (nodeIter != nodeIterEnd)
		{
			LocatableSetConstIter locatableIter = nodeIter->second.regions.begin();
			while (locatableIter != nodeIter->second.regions.end())
			{
				overlapping.insert(*locatableIter);
				
				locatableIter++;
			}
			
			nodeIter++;
		}
	}
	
	bool Contains(const Region& region) const
	{
		if (mBreakNodes.size() < 1)
		{
			return;
		}
		
		if (region.end < mBreakNodes.begin()->second.start)
		{
			return;
		}
		
		Region searchRegion;
		searchRegion.start = max(region.start, mBreakNodes.begin()->second.start);
		searchRegion.end = region.end;

		BreakNodeMapConstIter nodeIter = mBreakNodes.upper_bound(searchRegion.start);
		nodeIter--;
		
		if (nodeIter->first != searchRegion.start)
		{
			return false;
		}
		
		LocatableSetConstIter locatableIter = nodeIter->second.regions.begin();
		while (locatableIter != nodeIter->second.regions.end())
		{
			TLocatable* locatable = *locatableIter;
			
			if (locatable->GetRegion().start == searchRegion.start && locatable->GetRegion().end == searchRegion.end)
			{
				return true;
			}
			
			locatableIter++;
		}
		
		return false;
	}
	
	void PrintRegions() const
	{
		BreakNodeMapConstIter nodeIter = mBreakNodes.begin();
		BreakNodeMapConstIter nodeIterEnd = mBreakNodes.end();
		while (nodeIter != nodeIterEnd)
		{
			cerr << nodeIter->second.start << " ";
			
			if (nodeIter->second.size < 0)
			{
				cerr << "inf";
			}
			else
			{
				cerr << nodeIter->second.start + nodeIter->second.size - 1;
			}
			
			cerr << ":";
			
			LocatableSetConstIter locatableIter = nodeIter->second.regions.begin();
			while (locatableIter != nodeIter->second.regions.end())
			{
				const TLocatable* locatable = *locatableIter;

				cerr << " " << locatable->GetRegion().start << "-" << locatable->GetRegion().end;
				
				locatableIter++;
			}
			
			cerr << endl;
			
			nodeIter++;
		}
	}
	
private:
	typedef set<TLocatable*> LocatableSet;
	typedef typename LocatableSet::iterator LocatableSetIter;
	typedef typename LocatableSet::const_iterator LocatableSetConstIter;
	
	struct BreakNode
	{
		int start;
		int size;
		LocatableSet regions;
		friend struct CompareFunctor;
	};
	
	typedef map<int,BreakNode> BreakNodeMap;
	typedef typename BreakNodeMap::iterator BreakNodeMapIter;
	typedef typename BreakNodeMap::const_iterator BreakNodeMapConstIter;
	
	BreakNodeMap mBreakNodes;
};

#endif

