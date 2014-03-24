/*
 *  PathFinder.cpp
 *  findseqpath
 *
 *  Created by Andrew McPherson on 11/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "PathFinder.h"
#include "DebugCheck.h"

#include <iostream>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

using namespace std;
using namespace boost;

struct Extension
{
	int count;
	int nucleotide;
	int kmerIndex;
};

struct ExtensionCountOrder
{
	bool operator()(const Extension& ext1, const Extension& ext2)
	{
		return ext1.count < ext2.count;
	}
};

typedef vector<Extension> ExtensionVec;
typedef vector<Extension>::iterator ExtensionVecIter;
typedef vector<Extension>::const_iterator ExtensionVecConstIter;
typedef vector<ExtensionVec> ExtensionTable;

typedef unordered_map<string,int> SequenceMap;
typedef unordered_map<string,int>::iterator SequenceMapIter;
typedef unordered_map<string,int>::const_iterator SequenceMapConstIter;

enum Nucleotide
{
	NucleotideA = 0,
	NucleotideC = 1,
	NucleotideT = 2,
	NucleotideG = 3,
	
	NucleotideFirst = NucleotideA,
	NucleotideLast = NucleotideG,
};

char TranslateNucleotide(Nucleotide nucleotide)
{
	switch (nucleotide)
	{
		case NucleotideA: return 'A';
		case NucleotideC: return 'C';
		case NucleotideT: return 'T';
		case NucleotideG: return 'G';
	}
	
	return 'N';
}

void ExtendNucleotide(const string& kmer, Nucleotide nucleotide, string& extendedKmer)
{
	extendedKmer = kmer.substr(1);
	extendedKmer.push_back(TranslateNucleotide(nucleotide));
}

bool PathFinder::FindPath(const SequenceVec& sequences, const IntegerVec& leftSeq, const IntegerVec& rightSeq)
{
	DebugCheck(sequences.size() == leftSeq.size());
	DebugCheck(sequences.size() == rightSeq.size());
	
	mSequences.clear();

	SequenceMap vertexLookup;
	SequenceVec vertices;
	IntegerVec counts;
	ExtensionTable edges;

	IntegerVec leftCounts;
	IntegerVec rightCounts;
	
	// Create vertices
	for (int seqIndex = 0; seqIndex < sequences.size(); seqIndex++)
	{
		const string& sequence = sequences[seqIndex];
		bool left = (leftSeq[seqIndex] != 0);
		bool right = (rightSeq[seqIndex] != 0);

		DebugCheck(sequence.length() >= mK);
		DebugCheck(!left || !right);

		if (sequence.length() < mK)
		{
			continue;
		}
		
		for (int start = 0; start <= sequence.size() - mK; start++)
		{
			string kmer = sequence.substr(start, mK);			

			pair<SequenceMapIter,bool> lookupResult = vertexLookup.insert(make_pair(kmer,0));
			if (lookupResult.second)
			{
				lookupResult.first->second = vertices.size();
				
				vertices.push_back(kmer);
				counts.push_back(0);
				edges.push_back(ExtensionVec());
				
				leftCounts.push_back(0);
				rightCounts.push_back(0);
			}
			
			counts[lookupResult.first->second]++;
			
			if (left)
			{
				leftCounts[lookupResult.first->second]++;
			}
			
			if (right)
			{
				rightCounts[lookupResult.first->second]++;
			}
		}
	}
	
	// Create edges
	for (int kmerIndex = 0; kmerIndex < edges.size(); kmerIndex++)
	{
		const string& kmer = vertices[kmerIndex];
		
		for (int nucleotide = NucleotideFirst; nucleotide <= NucleotideLast; nucleotide++)
		{
			string extendedKmer;
			ExtendNucleotide(kmer, (Nucleotide)nucleotide, extendedKmer);
			
			SequenceMapConstIter extendedSeqIter = vertexLookup.find(extendedKmer);
			if (extendedSeqIter != vertexLookup.end())
			{
				Extension extension;
				extension.kmerIndex = extendedSeqIter->second;
				extension.count = counts[extendedSeqIter->second];
				extension.nucleotide = nucleotide;
				
				edges[kmerIndex].push_back(extension);
			}
		}
		
		sort(edges[kmerIndex].begin(), edges[kmerIndex].end(), ExtensionCountOrder());	
	}
	
	// Find highest count left and right sequences
	int startKmerIndex = -1;
	int startKmerMaxCount = 0;
	int endKmerIndex = -1;
	int endKmerMaxCount = 0;
	for (int kmerIndex = 0; kmerIndex < vertices.size(); kmerIndex++)
	{
		if (leftCounts[kmerIndex] > startKmerMaxCount)
		{
			startKmerMaxCount = leftCounts[kmerIndex];
			startKmerIndex = kmerIndex;
		}

		if (rightCounts[kmerIndex] > endKmerMaxCount)
		{
			endKmerMaxCount = rightCounts[kmerIndex];
			endKmerIndex = kmerIndex;
		}
	}
	
	if (startKmerIndex == -1)
	{
		cerr << "Error: Unable to find start kmer for finding paths" << endl;
		return false;
	}

	if (endKmerIndex == -1)
	{
		cerr << "Error: Unable to find end kmer for finding paths" << endl;
		return false;
	}
	
	// Find paths
	IntegerVec visited(vertices.size(), 0);
	
	SequenceVec pathSeqStack;
	IntegerVec pathKmerIndexStack;
	IntegerVec pathMinCount;
	
	pathSeqStack.push_back(vertices[startKmerIndex]);
	pathKmerIndexStack.push_back(startKmerIndex);
	pathMinCount.push_back(counts[startKmerIndex]);
	visited[startKmerIndex] = 1;
	
	while (!pathSeqStack.empty())
	{
		string currentSeq = pathSeqStack.back();
		int currentKmerIndex = pathKmerIndexStack.back();
		int currentMinCount = pathMinCount.back();
		
		pathSeqStack.pop_back();
		pathKmerIndexStack.pop_back();
		pathMinCount.pop_back();
		
		for (ExtensionVecConstIter extIter = edges[currentKmerIndex].begin(); extIter != edges[currentKmerIndex].end(); extIter++)
		{
			string nextSeq = currentSeq + TranslateNucleotide((Nucleotide)extIter->nucleotide);
			int nextKmerIndex = extIter->kmerIndex;
			int nextMinCount = min(currentMinCount, counts[nextKmerIndex]);
			
			if (nextMinCount < mMinCount)
			{
				continue;
			}
						
			if (nextKmerIndex == endKmerIndex)
			{
				mSequences.push_back(nextSeq);
				mMinCounts.push_back(nextMinCount);
				
				if (!mFindAll)
				{
					return true;
				}
			}
			else if (nextSeq.size() <= mMaxLength && !visited[nextKmerIndex])
			{
				pathSeqStack.push_back(nextSeq);
				pathKmerIndexStack.push_back(nextKmerIndex);
				pathMinCount.push_back(nextMinCount);
				visited[nextKmerIndex] = 1;
			}
		}
	}

	return !mSequences.empty();
}

