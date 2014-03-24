/*
 *  Indexer.h
 *  findtransbreaks
 *
 *  Created by Andrew McPherson on 26/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef INDEXER_H_
#define INDEXER_H_

#include "Common.h"

#include <vector>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

template <typename TType>
class Indexer
{
public:
	typedef vector<TType> DataVec;

	Indexer()
	{}
	
	Indexer(const DataVec& data)
	{
		mData = data;
		for (int index = 0; index < mData.size(); index++)
		{
			mLookup[mData[index]] = index;
		}
	}
	
	int Index(const TType& data, bool& inserted)
	{
		pair<IndexMapConstIter,bool> dataInsertResult = mLookup.insert(make_pair(data,mData.size()));

		inserted = dataInsertResult.second;
		int index = dataInsertResult.first->second;
		
		if (inserted)
		{
			mData.push_back(data);
		}

		return index;
	}
	
	int Index(const TType& data)
	{
		bool inserted;
		return Index(data,inserted);
	}

	int GetIndex(const TType& data) const
	{
		return mLookup.find(data)->second;
	}
	
	bool Contains(const TType& data)
	{
		return mLookup.find(data) != mLookup.end();
	}
	
	int Size() const
	{
		return mData.size();
	}
	
	const TType& Get(int index) const
	{
		return mData[index];
	}
	
	const DataVec& Get() const
	{
		return mData;
	}

	void Clear()
	{
		mData.clear();
		mLookup.clear();
	}
	
private:
	typedef unordered_map<TType,int> IndexMap;
	typedef typename unordered_map<TType,int>::const_iterator IndexMapConstIter;
	
	DataVec mData;
	IndexMap mLookup;
};

typedef Indexer<string> NameIndex;
typedef Indexer<IntegerPair> PairIndex;

#endif
