/*
 *  ExonRegions.cpp
 *  findfusions
 *
 *  Created by Andrew McPherson on 10-07-14.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "ExonRegions.h"

#include <boost/unordered_set.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

bool ExonRegions::ReadExonRegions(istream& in)
{
	string line;
	int lineNumber = 0;
	
	unordered_set<string> uniqueGeneNames;
	
	while (getline(in, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> exonRegionsFields;
		split(exonRegionsFields, line, is_any_of("\t"));
		
		if (exonRegionsFields.size() < 5)
		{
			continue;
		}
		
		string reference = exonRegionsFields[0];
		string chromosome = exonRegionsFields[1];
		string strand = exonRegionsFields[2];
		
		RegionVec exons;
		
		for (int fieldIndex = 4; fieldIndex < exonRegionsFields.size(); fieldIndex += 2)
		{
			Region exon;
			exon.start = lexical_cast<int>(exonRegionsFields[fieldIndex - 1]);
			exon.end = lexical_cast<int>(exonRegionsFields[fieldIndex]);
			
			exons.push_back(exon);
		}
		
		string::size_type geneEnd = reference.find_first_of('|');
		string geneName = reference.substr(0, geneEnd);
		
		int strandEnum;
		if (strand == "+")
		{
			strandEnum = PlusStrand;
		}
		else if (strand == "-")
		{
			strandEnum = MinusStrand;
		}
		else
		{
			cerr << "Error: Unable to intepret strand " << strand << endl;
			exit(1);			
		}
		
		int exonsLength = RegionsLength(exons);
		
		mChromosome[reference] = chromosome;
		mStrand[reference] = strandEnum;
		mExons[reference] = exons;
		mLength[reference] = exonsLength;
		
		uniqueGeneNames.insert(geneName);
		mGeneTranscripts[geneName].push_back(reference);
		mTranscriptGene[reference] = geneName;
	}
	
	mGenes = StringVec(uniqueGeneNames.begin(), uniqueGeneNames.end());

	return true;
}

bool ExonRegions::ReadGeneTranscripts(istream& in)
{
	mGeneTranscripts.clear();
	mTranscriptGene.clear();

	unordered_set<string> uniqueGeneNames;

	string line;
	int lineNumber = 0;
	
	while (getline(in, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> geneTranscriptFields;
		split(geneTranscriptFields, line, is_any_of("\t"));
		
		if (geneTranscriptFields.size() < 3)
		{
			continue;
		}
		
		string gene = geneTranscriptFields[0];
		string transcript = geneTranscriptFields[1];
		int length = lexical_cast<int>(geneTranscriptFields[2]);
		
		uniqueGeneNames.insert(gene);
		mGeneTranscripts[gene].push_back(transcript);
		mTranscriptGene[transcript] = gene;
		mLength[transcript] = length;
	}
	
	mGenes = StringVec(uniqueGeneNames.begin(), uniqueGeneNames.end());
	
	return true;
}

const StringVec& ExonRegions::GetGeneNames() const
{
	return mGenes;
}

int ExonRegions::GetTranscriptLength(const string& transcript) const
{
	if (mLength.find(transcript) == mLength.end())
	{
		cerr << "Error: Data mismatch, unable to find length for transcript " << transcript << endl;
		exit(1);
	}
	
	return mLength.find(transcript)->second;
}

const string& ExonRegions::GetTranscriptGeneName(const string& transcript) const
{
	if (mTranscriptGene.find(transcript) == mTranscriptGene.end())
	{
		cerr << "Error: Data mismatch, unable to find gene name for transcript " << transcript << endl;
		exit(1);
	}
	
	return mTranscriptGene.find(transcript)->second;
}

const StringVec& ExonRegions::GetGeneTranscriptNames(const string& gene) const
{
	if (mGeneTranscripts.find(gene) == mGeneTranscripts.end())
	{
		cerr << "Error: Data mismatch, unable to find transcript names for gene " << gene << endl;
		exit(1);
	}

	return mGeneTranscripts.find(gene)->second;
}

bool ExonRegions::IsTranscript(const string& transcript) const
{
	if (mTranscriptGene.find(transcript) != mTranscriptGene.end())
	{
		return true;
	}
	
	return false;
}

bool ExonRegions::FindOverlap(const Region& region1, const Region& region2, Region& overlap) const
{
	if (region1.start >= region2.end || region1.end <= region2.start)
	{
		return false;
	}
	
	overlap.start = max(region1.start, region2.start);
	overlap.end = min(region1.end, region2.end);
	
	return true;
}

void ExonRegions::MergeRegions(const RegionVec& regions, RegionVec& merged) const
{
	if (regions.size() == 0)
	{
		return;
	}
	
	Region merge;
	merge.start = regions.front().start;
	merge.end = regions.front().end;
	for (RegionVecConstIter regionIter = regions.begin(); regionIter != regions.end(); regionIter++)
	{
		if (regionIter->start > merge.end + 1)
		{
			merged.push_back(merge);
			
			merge.start = regionIter->start;
			merge.end = regionIter->end;
		}
		else
		{
			merge.end = max(merge.end, regionIter->end);
		}
	}

	merged.push_back(merge);
}

int ExonRegions::RegionsLength(const RegionVec& regions) const
{
	int length = 0;
	for (RegionVecConstIter regionIter = regions.begin(); regionIter != regions.end(); regionIter++)
	{
		length += regionIter->end - regionIter->start + 1;
	}
	
	return length;
}

bool ExonRegions::RemapTranscriptToGenome(const string& transcript, int position, int& remapped) const
{
	if (mExons.find(transcript) == mExons.end() || mExons.find(transcript)->second.size() == 0)
	{
		cerr << "Error: Data mismatch, unable to find transcript " << transcript << endl;
		exit(1);
	}
	
	const RegionVec& transcriptExons = mExons.find(transcript)->second;
	int transcriptLength = mLength.find(transcript)->second;
	int transcriptStrand = mStrand.find(transcript)->second;
	
	// Strand remapping of transcript position
	if (transcriptStrand == MinusStrand)
	{
		position = transcriptLength - position + 1;
	}
	
	remapped = -1;
	int localOffset = 0;
	for (RegionVecConstIter exonIter = transcriptExons.begin(); exonIter != transcriptExons.end(); exonIter++)
	{
		int exonLength = exonIter->end - exonIter->start + 1;
		
		// Local position of current exon
		Region localExon;
		localExon.start = localOffset + 1;
		localExon.end = localOffset + exonLength;
		
		// Test if the position is within or before this exon
		if (position <= localExon.end)
		{
			remapped = position - localExon.start + exonIter->start;
			return true;
		}
		
		localOffset += exonLength;
	}

	remapped = position - transcriptLength + transcriptExons.back().end;
	return true;
}

bool ExonRegions::RemapTranscriptToGene(const string& transcript, int position, string& gene, int& remapped) const
{
	int genomeRemapped;
	if (!RemapTranscriptToGenome(transcript, position, genomeRemapped))
	{
		return false;
	}
	
	gene = GetTranscriptGeneName(transcript);
	if (gene == "" || mExons.find(gene) == mExons.end())
	{
		return false;
	}
	
	int geneStrand = mStrand.find(gene)->second;
	
	if (mExons.find(gene) == mExons.end())
	{
		return false;
	}
	
	const RegionVec& geneRegions = mExons.find(gene)->second;
	
	if (geneRegions.size() != 1)
	{
		return false;
	}
	
	if (geneStrand == PlusStrand)
	{
		remapped = genomeRemapped - geneRegions[0].start + 1;
	}
	else if (geneStrand == MinusStrand)
	{
		remapped = geneRegions[0].end - genomeRemapped + 1;
	}
	else
	{		
		return false;
	}

	
	return true;
}

bool ExonRegions::RemapGenomeToTranscript(const string& transcript, int position, int& remapped) const
{
	if (mExons.find(transcript) == mExons.end() || mExons.find(transcript)->second.size() == 0)
	{
		cerr << "Error: Data mismatch, unable to find transcript " << transcript << endl;
		exit(1);
	}
	
	const RegionVec& transcriptExons = mExons.find(transcript)->second;
	int transcriptLength = mLength.find(transcript)->second;
	int transcriptStrand = mStrand.find(transcript)->second;
	
	int localOffset = 0;
	for (RegionVecConstIter exonIter = transcriptExons.begin(); exonIter != transcriptExons.end(); exonIter++)
	{
		int exonLength = exonIter->end - exonIter->start + 1;
		
		if (position >= exonIter->start && position <= exonIter->end)
		{
			remapped = position - exonIter->start + localOffset + 1;
			
			// Strand remapping of remapped position
			if (transcriptStrand == MinusStrand)
			{
				remapped = transcriptLength - remapped + 1;
			}
			
			return true;
		}
		
		localOffset += exonLength;
	}
	
	return false;
}

bool ExonRegions::RemapGenomeToTranscripts(const string& gene, int position, StringVec& transcripts, IntegerVec& remapped) const
{
	transcripts = GetGeneTranscriptNames(gene);
	
	for (StringVecConstIter transcriptNameIter = transcripts.begin(); transcriptNameIter != transcripts.end(); transcriptNameIter++)
	{
		const string& transcript = *transcriptNameIter;
		
		int transcriptRemapped;
		if (RemapGenomeToTranscript(transcript, position, transcriptRemapped))
		{
			transcripts.push_back(transcript);
			remapped.push_back(transcriptRemapped);
		}
	}
	
	return true;
}

bool ExonRegions::RemapTranscriptToTranscript(const string& transcript, int position, const string& remapTranscript, int& remapped) const
{
	int genomeRemapped;
	if (!RemapTranscriptToGenome(transcript, position, genomeRemapped))
	{
		return false;
	}

	if (!RemapGenomeToTranscript(remapTranscript, genomeRemapped, remapped))
	{
		return false;
	}
	
	return true;
}

bool ExonRegions::TrimTranscriptRegion(const string& transcript, Region& region) const
{
	if (mLength.find(transcript) == mLength.end())
	{
		cerr << "Error: Data mismatch, unable to find transcript " << transcript << endl;
		exit(1);
	}
	
	int transcriptLength = mLength.find(transcript)->second;
	
	region.start = max(1,region.start);
	region.end = min(transcriptLength,region.end);
	
	return true;
}

bool ExonRegions::RemapTranscriptToGene(const string& transcript, Region& region, RegionVec& remapped) const
{
	if (mExons.find(transcript) == mExons.end() || mExons.find(transcript)->second.size() == 0)
	{
		cerr << "Error: Data mismatch, unable to find transcript " << transcript << endl;
		exit(1);
	}
	
	const RegionVec& transcriptExons = mExons.find(transcript)->second;
	int transcriptLength = mLength.find(transcript)->second;
	int transcriptStrand = mStrand.find(transcript)->second;

	const string& gene = GetTranscriptGeneName(transcript);
	if (gene == "" || mExons.find(gene) == mExons.end())
	{
		cerr << "Error: Data mismatch, unable to find gene for transcript " << transcript << endl;
		exit(1);
	}
	
	int geneLength = mLength.find(gene)->second;
	int geneStrand = mStrand.find(gene)->second;
	
	if (mExons.find(gene) == mExons.end() || mExons.find(gene)->second.size() == 0)
	{
		cerr << "Error: Data mismatch, unable to find region for gene " << gene << endl;
		exit(1);
	}
	
	const RegionVec& geneExons = mExons.find(gene)->second;

	// Strand remapping of start and end
	if (transcriptStrand == MinusStrand)
	{
		region.start = transcriptLength - region.start + 1;
		region.end = transcriptLength - region.end + 1;
	}	
	
	// Remap transcript regions to genome
	RegionVec genomeRegions;
	int localOffset = 0;
	for (RegionVecConstIter exonIter = transcriptExons.begin(); exonIter != transcriptExons.end(); exonIter++)
	{
		int exonLength = exonIter->end - exonIter->start + 1;
		
		// Local position of current exon
		Region localExon;
		localExon.start = localOffset + 1;
		localExon.end = localOffset + exonLength;
		
		// Test for overlap between current exon and region
		Region overlap;
		if (FindOverlap(region, localExon, overlap))
		{
			// Calculate projection of region from this exon onto the genome
			Region genomeRegion;
			genomeRegion.start = exonIter->start + overlap.start - localOffset - 1;
			genomeRegion.end = exonIter->start + overlap.end - localOffset - 1;
			
			genomeRegions.push_back(genomeRegion);
		}
		
		localOffset += exonLength;
	}
	
	// Remap genome regions to gene
	remapped.clear();
	for (RegionVecConstIter genomeRegionIter = genomeRegions.begin(); genomeRegionIter != genomeRegions.end(); genomeRegionIter++)
	{
		Region geneRegion;
		geneRegion.start = genomeRegionIter->start - geneExons[0].start + 1;
		geneRegion.end = genomeRegionIter->end - geneExons[0].start + 1;
		
		// Strand remapping of start and end
		if (geneStrand == MinusStrand)
		{
			geneRegion.start = geneLength - geneRegion.start + 1;
			geneRegion.end = geneLength - geneRegion.end + 1;
		}
		
		remapped.push_back(geneRegion);
	}
	
	// Strand remapping, reverse exons
	if (transcriptStrand == MinusStrand)
	{
		reverse(remapped.begin(), remapped.end());
	}
	
	return true;
}


