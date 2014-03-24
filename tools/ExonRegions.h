/*
 *  ExonRegions.h
 *  findfusions
 *
 *  Created by Andrew McPherson on 10-07-14.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef EXONREGIONS_H_
#define EXONREGIONS_H_

#include "Common.h"

#include <string>
#include <iostream>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

class ExonRegions
{
public:
	bool ReadExonRegions(istream& in);
	bool ReadGeneTranscripts(istream& in);
	
	const StringVec& GetGeneNames() const;
	int GetTranscriptLength(const string& transcript) const;
	
	const string& GetTranscriptGeneName(const string& transcript) const;
	const StringVec& GetGeneTranscriptNames(const string& gene) const;
	
	bool IsTranscript(const string& transcript) const;
	bool RemapTranscriptToGenome(const string& transcript, int position, int& remapped) const;
	bool RemapTranscriptToGene(const string& transcript, int position, string& gene, int& remapped) const;
	bool RemapGenomeToTranscript(const string& transcript, int position, int& remapped) const;
	bool RemapGenomeToTranscripts(const string& gene, int position, StringVec& transcripts, IntegerVec& remapped) const;
	bool RemapTranscriptToTranscript(const string& transcript, int position, const string& remapTranscript, int& remapped) const;
	bool TrimTranscriptRegion(const string& transcript, Region& region) const;
	bool RemapTranscriptToGene(const string& transcript, Region& region, RegionVec& remapped) const;
	
private:
	bool FindOverlap(const Region& region1, const Region& region2, Region& overlap) const;
	void MergeRegions(const RegionVec& regions, RegionVec& merged) const;
	int RegionsLength(const RegionVec& regions) const;
	
	StringVec mGenes;
	unordered_map<string,string> mChromosome;
	unordered_map<string,int> mStrand;
	unordered_map<string,RegionVec> mExons;
	unordered_map<string,int> mLength;
	unordered_map<string,StringVec> mGeneTranscripts;
	unordered_map<string,string> mTranscriptGene;
};

#endif
