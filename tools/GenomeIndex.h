/*
 *  GenomeIndex.h
 *  linkexons
 *
 *  Created by Andrew McPherson on 11/08/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GENOME_INDEX_H
#define GENOME_INDEX_H

#include "RegionIndex.h"

#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

template <class TLocatable>
class GenomeIndex
{
public:
	typedef typename TLocatable::ChromosomeType ChromosomeType;
	
	GenomeIndex() {}
	
	void Insert(TLocatable* locatable)
	{
		pair<RegionMapIter,bool> insRes = mChromosomeRegionIndices.insert(make_pair(locatable->GetChromosome(),(TRegionIndex*)0));
		
		if (insRes.second)
		{
			insRes.first->second = new TRegionIndex();
		}
		
		insRes.first->second->Insert(locatable);
	}
	
	bool Delete(TLocatable* locatable)
	{
		RegionMapIter chrIter = mChromosomeRegionIndices.find(locatable->GetChromosome());
		
		if (chrIter == mChromosomeRegionIndices.end())
		{
			return false;
		}
		
		return chrIter->second->Delete(locatable);		
	}
	
	void Clear()
	{
		for (RegionMapIter chrIter = mChromosomeRegionIndices.begin(); chrIter != mChromosomeRegionIndices.end(); chrIter++)
		{
			chrIter->second->Clear();
		}
	}
	
	void GetOverlapping(const ChromosomeType& chromosome, const Region& region, set<TLocatable*>& overlapping) const
	{
		RegionMapConstIter regionIter = mChromosomeRegionIndices.find(chromosome);
		
		if (regionIter != mChromosomeRegionIndices.end())
		{
			regionIter->second->GetOverlapping(region,overlapping);
		}
	}		
	
	template <class TOtherLocatable>
	void GetOverlapping(const TOtherLocatable& locatable, set<TLocatable*>& regions) const
	{
		GetOverlapping(locatable.GetChromosome(), locatable.GetRegion(), regions);
	}
	
	bool Contains(const ChromosomeType& chromosome, const Region& region) const
	{
		RegionMapConstIter regionIter = mChromosomeRegionIndices.find(chromosome);
		
		if (regionIter != mChromosomeRegionIndices.end())
		{
			return regionIter->second->Contains(region);
		}

		return false;
	}		
	
	template <class TOtherLocatable>
	bool Contains(const TOtherLocatable& locatable) const
	{
		return Contains(locatable.GetChromosome(), locatable.GetRegion());
	}
	
private:
	typedef RegionIndex<TLocatable> TRegionIndex;
	typedef unordered_map<ChromosomeType,TRegionIndex*> RegionMap;
	typedef typename RegionMap::iterator RegionMapIter;
	typedef typename RegionMap::const_iterator RegionMapConstIter;

	RegionMap mChromosomeRegionIndices;
};

#endif
