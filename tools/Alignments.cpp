/*
 *  Alignments.cpp
 *  linkexons
 *
 *  Created by Andrew McPherson on 09-09-12.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Alignments.h"
#include "DebugCheck.h"

#include <iostream>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;

bool operator==(const Alignment& a1, const Alignment& a2)
{
	return a1.readID.id == a2.readID.id && a1.refStrand.id == a2.refStrand.id && a1.region.start == a2.region.start && a1.region.end == a2.region.end;
}

size_t hash_value(const Alignment& a)
{
	size_t seed = 0;
	hash_combine(seed, a.readID.id);
	hash_combine(seed, a.refStrand.id);
	hash_combine(seed, a.region.start);
	hash_combine(seed, a.region.end);
	return seed;
}

NameIndex Alignments::FragmentNames;
NameIndex Alignments::ReferenceNames;

Alignments::Alignments()
{
}

bool Alignments::Read(IAlignmentStream* alignmentStream)
{
	RawAlignment rawAlignment;
	while (alignmentStream->GetNextAlignment(rawAlignment))
	{
		Alignment alignment;
		
		alignment.readID.readEnd = rawAlignment.readEnd;
		alignment.readID.fragmentIndex = FragmentNames.Index(rawAlignment.fragment);

		alignment.refStrand.referenceIndex = ReferenceNames.Index(rawAlignment.reference);
		alignment.refStrand.strand = rawAlignment.strand;
		
		alignment.region = rawAlignment.region;
		
		mAlignments.push_back(alignment);
	}
	
	IndexByReadID();
	
	return true;
}

void Alignments::RemoveRedundant()
{
	unordered_set<Alignment> alignmentSet;
	
	AlignmentVec alignments;
	
	for (int alignmentIndex = 0; alignmentIndex < (int)mAlignments.size(); alignmentIndex++)
	{
		if (alignmentSet.insert(mAlignments[alignmentIndex]).second)
		{
			alignments.push_back(mAlignments[alignmentIndex]);
		}
	}
	
	swap(alignments,mAlignments);
	
	IndexByReadID();
}

void Alignments::RemoveUnPaired()
{
	IndexByReadID();
	
	AlignmentVec alignments;
	
	for (int fragmentIndex = 0; fragmentIndex < FragmentNames.Size(); fragmentIndex++)
	{
		ReadID end1;
		end1.fragmentIndex = fragmentIndex;
		end1.readEnd = 0;
		
		const IntegerVec& end1AlignIndices = GetAlignments(end1);
		
		ReadID end2;
		end2.fragmentIndex = fragmentIndex;
		end2.readEnd = 1;
		
		const IntegerVec& end2AlignIndices = GetAlignments(end2);
		
		if (end1AlignIndices.size() == 0 || end2AlignIndices.size() == 0)
		{
			continue;
		}

		for (IntegerVecConstIter end1AlignIndexIter = end1AlignIndices.begin(); end1AlignIndexIter != end1AlignIndices.end(); end1AlignIndexIter++)
		{
			int end1AlignIndex = *end1AlignIndexIter;
			const Alignment& end1Alignment = mAlignments[end1AlignIndex];
			alignments.push_back(end1Alignment);
		}
		
		for (IntegerVecConstIter end2AlignIndexIter = end2AlignIndices.begin(); end2AlignIndexIter != end2AlignIndices.end(); end2AlignIndexIter++)
		{
			int end2AlignIndex = *end2AlignIndexIter;
			const Alignment& end2Alignment = mAlignments[end2AlignIndex];
			alignments.push_back(end2Alignment);
		}
	}

	swap(alignments,mAlignments);
	
	IndexByReadID();
}

void Alignments::FilterEnsemblTranscripts()
{
	IndexByReadID();
	
	AlignmentVec alignments;

	for (int fragmentIndex = 0; fragmentIndex < FragmentNames.Size(); fragmentIndex++)
	{
		ReadID end1;
		end1.fragmentIndex = fragmentIndex;
		end1.readEnd = 0;
		
		const IntegerVec& end1AlignIndices = GetAlignments(end1);
		
		ReadID end2;
		end2.fragmentIndex = fragmentIndex;
		end2.readEnd = 1;
		
		const IntegerVec& end2AlignIndices = GetAlignments(end2);
		
		bool sameGene = false;
		for (IntegerVecConstIter end1AlignIndexIter = end1AlignIndices.begin(); end1AlignIndexIter != end1AlignIndices.end(); end1AlignIndexIter++)
		{
			int end1AlignIndex = *end1AlignIndexIter;
			const Alignment& end1Alignment = mAlignments[end1AlignIndex];
			const string& end1RefName = ReferenceNames.Get(end1Alignment.refStrand.referenceIndex);
			
			string::size_type end1GeneEnd = end1RefName.find_first_of('|');
			string end1GeneName = end1RefName.substr(0, end1GeneEnd);			
			
			for (IntegerVecConstIter end2AlignIndexIter = end2AlignIndices.begin(); end2AlignIndexIter != end2AlignIndices.end(); end2AlignIndexIter++)
			{
				int end2AlignIndex = *end2AlignIndexIter;
				const Alignment& end2Alignment = mAlignments[end2AlignIndex];
				
				const string& end2RefName = ReferenceNames.Get(end2Alignment.refStrand.referenceIndex);
				
				string::size_type end2GeneEnd = end2RefName.find_first_of('|');
				string end2GeneName = end2RefName.substr(0, end2GeneEnd);
				
				if (end1GeneName == end2GeneName)
				{
					sameGene = true;
				}
			}
		}
		
		if (!sameGene)
		{
			for (IntegerVecConstIter end1AlignIndexIter = end1AlignIndices.begin(); end1AlignIndexIter != end1AlignIndices.end(); end1AlignIndexIter++)
			{
				int end1AlignIndex = *end1AlignIndexIter;
				alignments.push_back(mAlignments[end1AlignIndex]);
			}
				
			for (IntegerVecConstIter end2AlignIndexIter = end2AlignIndices.begin(); end2AlignIndexIter != end2AlignIndices.end(); end2AlignIndexIter++)
			{
				int end2AlignIndex = *end2AlignIndexIter;
				alignments.push_back(mAlignments[end2AlignIndex]);
			}
		}
	}
	
	swap(alignments,mAlignments);
	
	IndexByReadID();
}

void Alignments::IndexByReadID()
{
	mReadAlignments[0].clear();
	mReadAlignments[1].clear();
	
	mReadAlignments[0].resize(FragmentNames.Size());
	mReadAlignments[1].resize(FragmentNames.Size());
	
	for (int alignmentIndex = 0; alignmentIndex < (int)mAlignments.size(); alignmentIndex++)
	{
		const ReadID& readID = mAlignments[alignmentIndex].readID;
		mReadAlignments[readID.readEnd][readID.fragmentIndex].push_back(alignmentIndex);
	}
}

int Alignments::GetNumAlignments() const
{
	return (int)mAlignments.size();
}

const Alignment& Alignments::GetAlignment(int alignmentIndex) const
{
	return mAlignments[alignmentIndex];
}

const IntegerVec& Alignments::GetAlignments(const ReadID& readID) const
{
	DebugCheck(mReadAlignments[0].size() > 0);
	DebugCheck(mReadAlignments[1].size() > 0);
	
	return mReadAlignments[readID.readEnd][readID.fragmentIndex];
}

