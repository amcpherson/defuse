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

int ExonRegions::mBinLength = 100000;

bool ExonRegions::Read(istream& in)
{
	string line;
	int lineNumber = 0;
	
	unordered_set<string> uniqueGenes;
	
	while (getline(in, line))
	{
		lineNumber++;
					
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> exonRegionsFields;
		split(exonRegionsFields, line, is_any_of("\t"));
		
		if (exonRegionsFields.size() < 6)
		{
			continue;
		}
		
		string gene = exonRegionsFields[0];
		string transcript = exonRegionsFields[1];
		string chromosome = exonRegionsFields[2];
		string strand = exonRegionsFields[3];
		
		RegionVec exons;
		for (int fieldIndex = 5; fieldIndex < exonRegionsFields.size(); fieldIndex += 2)
		{
			try 
			{
				Region exon;
				exon.start = lexical_cast<int>(exonRegionsFields[fieldIndex - 1]);
				exon.end = lexical_cast<int>(exonRegionsFields[fieldIndex]);
				
				exons.push_back(exon);
			}
			catch (boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::bad_lexical_cast> > e)
			{
				cout << "Failed to interpret exon:" << endl << line << endl;
				exit(1);
			}
		}
		
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
		
		mChromosome[transcript] = chromosome;
		mStrand[transcript] = strandEnum;
		mExons[transcript] = exons;
		mLength[transcript] = exonsLength;
		
		uniqueGenes.insert(gene);
		mGeneTranscripts[gene].push_back(transcript);
		mTranscriptGene[transcript] = gene;
		
		mExonsStr[PlusStrand][transcript] = exons;
		mExonsStr[MinusStrand][transcript] = exons;
		TransformExons(mExonsStr[MinusStrand][transcript]);
		
		mTranscriptRegion[transcript].start = mExons[transcript].front().start;
		mTranscriptRegion[transcript].end = mExons[transcript].back().end;
		
		int startBin = mTranscriptRegion[transcript].start / mBinLength;
		int endBin = mTranscriptRegion[transcript].end / mBinLength;		
		for (int bin = startBin; bin <= endBin; bin++)
		{
			mTranscriptLookup[chromosome][bin].push_back(transcript);
		}
	}
	
	mGenes = StringVec(uniqueGenes.begin(), uniqueGenes.end());

	return true;
}

void ExonRegions::TransformExons(RegionVec& exons) const
{
	for (RegionVecIter exonIter = exons.begin(); exonIter != exons.end(); exonIter++)
	{
		exonIter->start = -exonIter->start;
		exonIter->end = -exonIter->end;
		swap(exonIter->start, exonIter->end);
	}
	
	reverse(exons.begin(), exons.end());
}

const StringVec& ExonRegions::GetGenes() const
{
	return mGenes;
}

void ExonRegions::GetRegionTranscripts(const string& chromosome, const Region& region, StringVec& transcripts) const
{
	if (mTranscriptLookup.find(chromosome) == mTranscriptLookup.end())
	{
		cerr << "Error: Data mismatch, invalid chromosome " << chromosome << endl;
		exit(1);
	}
	
	unordered_set<string> uniqueTranscripts;
	
	int startBin = region.start / mBinLength;
	int endBin = region.end / mBinLength;		
	for (int bin = startBin; bin <= endBin; bin++)
	{
		unordered_map<int,StringVec>::const_iterator lookupIter = mTranscriptLookup.find(chromosome)->second.find(bin);
		if (lookupIter == mTranscriptLookup.find(chromosome)->second.end())
		{
			continue;
		}
		
		for (StringVecConstIter transcriptIter = lookupIter->second.begin(); transcriptIter != lookupIter->second.end(); transcriptIter++)
		{
			if (Overlap(mTranscriptRegion.find(*transcriptIter)->second,region))
			{
				uniqueTranscripts.insert(*transcriptIter);
			}
		}
	}
	
	transcripts.insert(transcripts.end(), uniqueTranscripts.begin(), uniqueTranscripts.end());
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

const string& ExonRegions::GetTranscriptGene(const string& transcript) const
{
	if (mTranscriptGene.find(transcript) == mTranscriptGene.end())
	{
		cerr << "Error: Data mismatch, unable to find gene for transcript " << transcript << endl;
		exit(1);
	}
	
	return mTranscriptGene.find(transcript)->second;
}

const StringVec& ExonRegions::GetGeneTranscripts(const string& gene) const
{
	if (mGeneTranscripts.find(gene) == mGeneTranscripts.end())
	{
		cerr << "Error: Data mismatch, unable to find transcripts for gene " << gene << endl;
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

bool ExonRegions::RemapTranscriptToGenome(const string& transcript, int strand, int position, string& remapChromosome, int& remapStrand, int& remapPosition) const
{
	if (mExons.find(transcript) == mExons.end() || mExons.find(transcript)->second.size() == 0)
	{
		cerr << "Error: Data mismatch, unable to find transcript " << transcript << endl;
		exit(1);
	}
	
	const RegionVec& transcriptExons = mExons.find(transcript)->second;
	int transcriptLength = mLength.find(transcript)->second;
	int transcriptStrand = mStrand.find(transcript)->second;
	
	remapChromosome = mChromosome.find(transcript)->second;
	remapStrand = (transcriptStrand == strand) ? PlusStrand : MinusStrand;
	
	// Strand remapping of transcript position
	if (transcriptStrand == MinusStrand)
	{
		position = transcriptLength - position + 1;
	}
	
	remapPosition = -1;
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
			remapPosition = position - localExon.start + exonIter->start;
			return true;
		}
		
		localOffset += exonLength;
	}

	remapPosition = position - transcriptLength + transcriptExons.back().end;
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
/*
bool ExonRegions::RemapThroughTranscript(const string& transcript, int position, int strand, int min, int max, int& start, int& end) const
{
	if (mExons.find(transcript) == mExons.end() || mExons.find(transcript)->second.size() == 0)
	{
		cerr << "Error: Data mismatch, unable to find transcript " << transcript << endl;
		exit(1);
	}
	
	const RegionVec& transcriptExons = mExons.find(transcript)->second;
	int transcriptLength = mLength.find(transcript)->second;
	int transcriptStrand = mStrand.find(transcript)->second;
	
	if (position < transcriptExons.front().start && strand == PlusStrand)
	{
		start = min - transcriptExons.front().start + position;
		end = max - transcriptExons.front().start + position;
	}
	else if (position > transcriptExons.back().end && strand == MinusStrand)
	{
		start = transcriptLength + position - transcriptExons.back().end - max;
		end = transcriptLength + position - transcriptExons.back().end - min;
	}
	else
	{
		int localOffset = 0;
		RegionVecConstIter lastExonIter = transcriptExons.end();
		for (RegionVecConstIter exonIter = transcriptExons.begin(); exonIter != transcriptExons.end(); exonIter++)
		{
			int exonLength = exonIter->end - exonIter->start + 1;
			
			if (position < exonIter->start && lastExonIter != transcriptExons.end())
			{
				if (strand == PlusStrand)
				{
					start = localOffset + min - exonIter->start + position;
					end = localOffset + max - exonIter->start + position;
				}
				else
				{
					start = localOffset + position - lastExonIter->end - max;
					end = localOffset + position - lastExonIter->end - min;
				}
			}
			else if (position >= exonIter->start && position <= exonIter->end)
			{
				int remapped = position - exonIter->start + localOffset + 1;
				
				if (strand == PlusStrand)
				{
					start = remapped + min;
					end = remapped + max;
				}
				else
				{
					start = remapped - max;
					end = remapped - min;
				}
			}
			
			lastExonIter = exonIter;
			localOffset += exonLength;
		}
	}
	
	if (end < 1 || start > transcriptLength)
	{
		return false;
	}

	// Strand remapping of remapped position
	if (transcriptStrand == MinusStrand)
	{
		start = transcriptLength - start + 1;
		end = transcriptLength - end + 1;
		swap(start,end);
	}
	
	return true;
}
*/

bool ExonRegions::RemapThroughTranscript(const string& transcript, int position, int strand, int extendMin, int extendMax, int& remapStrand, int& start, int& end) const
{
	if (mExons.find(transcript) == mExons.end() || mExons.find(transcript)->second.size() == 0)
	{
		cerr << "Error: Data mismatch, unable to find transcript " << transcript << endl;
		exit(1);
	}
	
	const RegionVec& transcriptExons = mExonsStr[strand].find(transcript)->second;
	int transcriptLength = mLength.find(transcript)->second;
	int transcriptStrand = mStrand.find(transcript)->second;
	
	remapStrand = (strand == transcriptStrand) ? PlusStrand : MinusStrand;
	
	// Transform into strand space
	int strandPosition = (strand == PlusStrand) ? position : -position;
	
	if (strandPosition > transcriptExons.back().end)
	{
		return false;
	}
	
	int localOffset = 0;
	for (RegionVecConstIter exonIter = transcriptExons.begin(); exonIter != transcriptExons.end(); exonIter++)
	{
		int exonLength = exonIter->end - exonIter->start + 1;
		
		if (strandPosition <= exonIter->end)
		{
			int exonRelativeStart = strandPosition - exonIter->start + extendMin + 1;
			int exonRelativeEnd = strandPosition - exonIter->start + extendMax + 1;
			
			if (exonRelativeEnd < 1)
			{
				return false;
			}
			
			start = max(1,exonRelativeStart) + localOffset;
			end = max(1,exonRelativeEnd) + localOffset;
			
			break;
		}
		
		localOffset += exonLength;
	}
	
	// Check the region is in the transcript boundaries
	if (end < 1 || start > transcriptLength)
	{
		return false;
	}
	
	// Strand remapping of remapped position
	if (strand != transcriptStrand)
	{
		start = transcriptLength - start + 1;
		end = transcriptLength - end + 1;
		swap(start,end);
	}
	
	return true;
}

bool ExonRegions::RemapGenomeToTranscripts(const string& gene, int position, StringVec& transcripts, IntegerVec& remapped) const
{
	transcripts = GetGeneTranscripts(gene);
	
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

