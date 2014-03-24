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
	bool Read(istream& in);
	
	const StringVec& GetGenes() const;
	int GetTranscriptLength(const string& transcript) const;
	
	const string& GetTranscriptGene(const string& transcript) const;
	const StringVec& GetGeneTranscripts(const string& gene) const;
	void GetRegionTranscripts(const string& chromosome, const Region& region, StringVec& genes) const;
	
	bool IsTranscript(const string& transcript) const;
	bool RemapTranscriptToGenome(const string& transcript, int strand, int position, string& remapChromosome, int& remapStrand, int& remapPosition) const;
	bool RemapGenomeToTranscript(const string& transcript, int position, int& remapped) const;
	bool RemapGenomeToTranscripts(const string& gene, int position, StringVec& transcripts, IntegerVec& remapped) const;
	bool RemapThroughTranscript(const string& transcript, int position, int strand, int extendMin, int extendMax, int& remapStrand, int& start, int& end) const;
	
private:
	void TransformExons(RegionVec& exons) const;
	bool FindOverlap(const Region& region1, const Region& region2, Region& overlap) const;
	void MergeRegions(const RegionVec& regions, RegionVec& merged) const;
	int RegionsLength(const RegionVec& regions) const;
	
	static int mBinLength;
	
	StringVec mGenes;
	unordered_map<string,string> mChromosome;
	unordered_map<string,int> mStrand;
	unordered_map<string,RegionVec> mExons;
	unordered_map<string,int> mLength;
	unordered_map<string,StringVec> mGeneTranscripts;
	unordered_map<string,string> mTranscriptGene;
	unordered_map<string,RegionVec> mExonsStr[2];
	unordered_map<string,Region> mTranscriptRegion;
	unordered_map<string,unordered_map<int,StringVec> > mTranscriptLookup;
};

#endif

