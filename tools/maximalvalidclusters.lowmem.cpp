/*
 *  clustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Alignments.h"
#include "AlignmentIndex.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
#include "FragmentAlignmentsIterator.h"
#include "IMatePairClusterer.h"
#include "MatePairEM.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "bam.h"
#include "sam.h"

using namespace boost;
using namespace std;


struct RefBin
{
	RefBin() {}
	RefBin(int ref, int bin) : ref(ref), bin(bin) {}
	
	int ref;
	int bin;
	
	bool operator<(const RefBin& other) const
	{
		if (ref < other.ref)
		{
			return true;
		}
		else if (ref > other.ref)
		{
			return false;
		}
		
		return (bin < other.bin);
	}
};

inline bool operator==(const RefBin& rb1, const RefBin& rb2)
{
	return rb1.ref == rb2.ref && rb1.bin == rb2.bin;
}

inline size_t hash_value(const RefBin& rb)
{
	size_t seed = 0;
	hash_combine(seed, rb.ref);
	hash_combine(seed, rb.bin);
	return seed;
}

typedef vector<RefBin> RefBinVec;
typedef vector<RefBin>::const_iterator RefBinVecConstIter;

typedef unordered_set<RefBin> RefBinSet;
typedef unordered_set<RefBin>::const_iterator RefBinSetConstIter;

typedef unordered_map<RefBin,IntegerVec> RefBinVecMap;
typedef unordered_map<RefBin,IntegerVec>::const_iterator RefBinVecMapConstIter;

typedef pair<RefBin,RefBin> RefBinPair;
typedef unordered_map<RefBinPair,IntegerPairVec> RefBinPairVecMap;
typedef unordered_map<RefBinPair,IntegerPairVec>::const_iterator RefBinPairVecMapConstIter;

typedef unordered_map<RefBinPair,int> RefBinPairMap;
typedef unordered_map<RefBinPair,int>::const_iterator RefBinPairMapConstIter;

typedef unordered_map<int,CompAlignVec> CompAlignVecMap;
typedef unordered_map<int,CompAlignVec>::const_iterator CompAlignVecMapConstIter;

class Binning
{
public:
	Binning(int binLength, int regionExtend) : mBinLength(binLength), mRegionExtend(regionExtend)
	{}
	
	void GetBins(const Region& region, IntegerVec& bins)
	{
		int startBin = (region.start - mRegionExtend) / mBinLength;
		int endBin = (region.end + mRegionExtend) / mBinLength;
		
		for (int bin = startBin; bin <= endBin; bin++)
		{
			bins.push_back(bin);
		}
	}
	
	void GetBinRegion(int bin, Region& region)
	{
		region.start = bin * mBinLength - mRegionExtend;
		region.end = (bin + 1) * mBinLength - 1 + mRegionExtend;
	}
	
private:
	int mBinLength;
	int mRegionExtend;
};

bool CheckConcordant(const CompAlignVec& alignments, int minFusionRange)
{
	CompAlignVec filtered;
	
	// Check for concordance
	RefBinSet concordanceBins[2];
	Binning concordanceBinning(minFusionRange, 0);
	for (CompAlignVecConstIter alignIter = alignments.begin(); alignIter != alignments.end(); alignIter++)
	{
		const CompactAlignment& alignment = *alignIter;
		
		// Bin end alignments according to bins with length bounded by minimum fusion distance
		IntegerVec rangeBins;
		concordanceBinning.GetBins(alignment.region, rangeBins);
		
		// Insert alignments for concordance check
		for (IntegerVecConstIter rangeBinIter = rangeBins.begin(); rangeBinIter != rangeBins.end(); rangeBinIter++)
		{
			RefBin concordanceBin(alignment.refStrand.referenceIndex,*rangeBinIter);				
			concordanceBins[alignment.readID.readEnd].insert(concordanceBin);
		}
	}
	
	// Classify as concordant if both ends aligning to same region		
	for (RefBinSetConstIter concordanceBinIter = concordanceBins[0].begin(); concordanceBinIter != concordanceBins[0].end(); concordanceBinIter++)
	{
		if (concordanceBins[1].find(*concordanceBinIter) != concordanceBins[1].end())
		{
			return true;
		}
	}
	
	return false;
}

void AddBinPairs(const CompAlignVec& alignments, Binning& binning, RefBinPairMap& binPairs)
{
	// Create binned end alignments
	RefBinVec binnedEndAlignments[2];
	for (CompAlignVecConstIter alignIter = alignments.begin(); alignIter != alignments.end(); alignIter++)
	{
		const CompactAlignment& alignment = *alignIter;
		
		IntegerVec clusterBins;
		binning.GetBins(alignment.region, clusterBins);
		for (IntegerVecConstIter clusterBinIter = clusterBins.begin(); clusterBinIter != clusterBins.end(); clusterBinIter++)
		{
			RefBin refBin(alignment.refStrand.id,*clusterBinIter);
			binnedEndAlignments[alignment.readID.readEnd].push_back(refBin);
		}
	}
	
	// Create binned alignments
	for (RefBinVecConstIter refBin1Iter = binnedEndAlignments[0].begin(); refBin1Iter != binnedEndAlignments[0].end(); refBin1Iter++)
	{
		for (RefBinVecConstIter refBin2Iter = binnedEndAlignments[1].begin(); refBin2Iter != binnedEndAlignments[1].end(); refBin2Iter++)
		{
			RefBin refBin1 = *refBin1Iter;
			RefBin refBin2 = *refBin2Iter;
			
			if (refBin2 < refBin1)
			{
				swap(refBin1,refBin2);
			}
			
			RefBinPair refBinPair(refBin1,refBin2);
			
			if (binPairs.find(refBinPair) == binPairs.end())
			{
				binPairs[refBinPair] = 0;
			}
			
			binPairs[refBinPair]++;
		}
	}	
}

void GetFragmentAlignments(CompAlignVec& alignments, IntegerVecMap& fragments)
{
	for (int alignIndex = 0; alignIndex < alignments.size(); alignIndex++)
	{
		fragments[alignments[alignIndex].readID.fragmentIndex].push_back(alignIndex);		
	}
}

void FilterUnmatched(const IntegerVecMap& fragments1, IntegerVecMap& fragments2)
{
	IntegerVecMap matchedFragments2;
	for (IntegerVecMapIter fragment2Iter = fragments2.begin(); fragment2Iter != fragments2.end(); fragment2Iter++)
	{
		if (fragments1.find(fragment2Iter->first) != fragments1.end())
		{
			swap(matchedFragments2[fragment2Iter->first], fragment2Iter->second);
		}
	}
	
	swap(matchedFragments2, fragments2);
}

void FilterOverlapping(IntegerVecMap& fragments, const CompAlignVec& alignments, int minFusionRange)
{
	Binning overlapBinning(minFusionRange, 0);
	
	for (IntegerVecMapIter fragmentIter = fragments.begin(); fragmentIter != fragments.end(); fragmentIter++)
	{
		IntegerVec filtered;
		
		// Select a set of non overlapping end alignments
		RefBinSet alignmentBins[2];
		IntegerVec nonOverlapping;
		for (IntegerVecConstIter alignIter = fragmentIter->second.begin(); alignIter != fragmentIter->second.end(); alignIter++)
		{
			const CompactAlignment& alignment = alignments[*alignIter];
			
			IntegerVec rangeBins;
			overlapBinning.GetBins(alignment.region, rangeBins);
			
			// Check for a potential overlapping alignment of one end
			bool overlappingAlignment = false;
			for (IntegerVecConstIter rangeBinIter = rangeBins.begin(); rangeBinIter != rangeBins.end(); rangeBinIter++)
			{
				RefBin refIdRangeBin(alignment.refStrand.id,*rangeBinIter);				
				if (alignmentBins[alignment.readID.readEnd].find(refIdRangeBin) != alignmentBins[alignment.readID.readEnd].end())
				{
					overlappingAlignment = true;
				}
			}
			
			// Randomly select the first encountered of potential overlapping ends
			if (!overlappingAlignment)
			{
				for (IntegerVecConstIter rangeBinIter = rangeBins.begin(); rangeBinIter != rangeBins.end(); rangeBinIter++)
				{
					RefBin refIdRangeBin(alignment.refStrand.id,*rangeBinIter);				
					alignmentBins[alignment.readID.readEnd].insert(refIdRangeBin);
				}
				
				filtered.push_back(*alignIter);
			}
		}
		
		swap(filtered, fragmentIter->second);
	}
}

void GetAlignPairs(const IntegerVecMap& fragments1, const IntegerVecMap& fragments2, PairIndex& alignPairs)
{
	for (IntegerVecMapConstIter fragment1Iter = fragments1.begin(); fragment1Iter != fragments1.end(); fragment1Iter++)
	{
		const IntegerVec& alignments1 = fragment1Iter->second;
		const IntegerVec& alignments2 = fragments2.find(fragment1Iter->first)->second;
		
		for (IntegerVecConstIter align1Iter = alignments1.begin(); align1Iter != alignments1.end(); align1Iter++)
		{
			for (IntegerVecConstIter align2Iter = alignments2.begin(); align2Iter != alignments2.end(); align2Iter++)
			{
				alignPairs.Index(IntegerPair(*align1Iter, *align2Iter));
			}
		}
	}
}

void OutputClusterMember(ostream& out, int clusterID, int clusterEnd, const CompactAlignment& alignment, const StringVec& referenceNames)
{
	out << clusterID << "\t";
	out << clusterEnd << "\t";
	out << alignment.readID.fragmentIndex << "\t";
	out << alignment.readID.readEnd << "\t";
	out << referenceNames[alignment.refStrand.referenceIndex] << "\t";
	out << ((alignment.refStrand.strand == PlusStrand) ? "+" : "-") << "\t";
	out << alignment.region.start << "\t";
	out << alignment.region.end << endl;
}

int main(int argc, char* argv[])
{
	string readSortedBamFilename;
	string alignSortedBamFilename;
	string clustersFilename;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	double precision;
	int minClusterSize;

	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> readSortedBamFilenameArg("r","readsorted","Read Sorted Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> alignSortedBamFilenameArg("a","alignsorted","Alignment Sorted Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Output Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","fragmentmean","Fragment Length Mean",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","fragmentstddev","Fragment Length Standard Deviation",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> precisionArg("p","precision","Precision",true,-1,"double",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		cmd.parse(argc,argv);

		readSortedBamFilename = readSortedBamFilenameArg.getValue();
		alignSortedBamFilename = alignSortedBamFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		precision = precisionArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int minFusionRange = (int)(fragmentLengthMean + 10 * fragmentLengthStdDev);
	const int binLength = 10000;
		
	cout << "Finding pairs of reference sequences connected by pairs of alignments" << endl;
		
	FragmentAlignmentsIterator fragmentsIter;
	StringVec referenceNames;
	
	fragmentsIter.Open(readSortedBamFilename, referenceNames);

	NameIndex refNameIndex(referenceNames);

	AlignmentIndex alignmentIndex;
	alignmentIndex.Open(alignSortedBamFilename);
	
	RefBinPairMap binPairs;
	Binning binning(binLength, minFusionRange);
	CompAlignVec alignments;
	while (fragmentsIter.GetNext(alignments))
	{
		if (CheckConcordant(alignments, minFusionRange))
		{
			continue;
		}
		
		AddBinPairs(alignments, binning, binPairs);
	}
	fragmentsIter.Close();
	
	// Initialize clusterer
	cout << "Initializing clusterer" << endl;
	MatePairEM clusterer;
	IMatePairClusterer::Parameters clustererParams;
	clustererParams.fragmentMean = fragmentLengthMean;
	clustererParams.fragmentStdDev = fragmentLengthStdDev;
	clustererParams.precision = precision;
	clustererParams.minClusterSize = minClusterSize;
	clustererParams.resolution = 1;
	clusterer.Initialize(clustererParams);
	
	// Open output clusters file
	ofstream clustersFile(clustersFilename.c_str());
	if (!clustersFile)
	{
		cerr << "Error: unable to write to clusters file" << endl;		
		exit(1);
	}
	
	cout << "Creating clusters" << endl;
	
	int clusterID = 0;
	for (RefBinPairMapConstIter refBinIter = binPairs.begin(); refBinIter != binPairs.end(); refBinIter++)
	{
		RefBin refBin1 = refBinIter->first.first;
		RefBin refBin2 = refBinIter->first.second;
		
		if (refBinIter->second < minClusterSize)
		{
			continue;
		}
		
		RefStrand refStrand1;
		RefStrand refStrand2;
		Region region1;
		Region region2;
		
		cout << refBin1.ref <<"\t"<< refBin1.bin << "\t" << refBin2.ref << "\t" <<refBin2.bin << endl;
		
		refStrand1.id = refBin1.ref;
		refStrand2.id = refBin2.ref;
		
		binning.GetBinRegion(refBin1.bin, region1);		
		binning.GetBinRegion(refBin2.bin, region2);
		
		const string& refName1 = referenceNames[refStrand1.referenceIndex];
		const string& refName2 = referenceNames[refStrand2.referenceIndex];
		
		cout << referenceNames.size() << endl;;
		
		CompAlignVec alignments1;
		CompAlignVec alignments2;
		
		alignmentIndex.Find(refName1, refStrand1.strand, region1.start, region1.end, refNameIndex, alignments1);
		alignmentIndex.Find(refName2, refStrand2.strand, region2.start, region2.end, refNameIndex, alignments2);
		
		IntegerVecMap fragments1;
		IntegerVecMap fragments2;
		
		GetFragmentAlignments(alignments1, fragments1);
		GetFragmentAlignments(alignments2, fragments2);
		
		FilterUnmatched(fragments1, fragments2);
		FilterUnmatched(fragments2, fragments1);
		
		FilterOverlapping(fragments1, alignments1, minFusionRange);
		FilterOverlapping(fragments2, alignments2, minFusionRange);
		
		PairIndex alignPairs;
		GetAlignPairs(fragments1, fragments2, alignPairs);
		
		IMatePairClusterer::MatePairVec matePairs;
		RegionVec breakRegions1;
		RegionVec breakRegions2;
		IntegerVec fragmentIndices;

		cout << fragments1.size() << endl;
		cout << fragments2.size() << endl;
		cout << alignments1.size() << endl;
		cout << alignments2.size() << endl;
		cout << alignPairs.Size() << endl;
		
		for (int alignPairIndex = 0; alignPairIndex < alignPairs.Size(); alignPairIndex++)
		{
			int alignIndex1 = alignPairs.Get(alignPairIndex).first;
			int alignIndex2 = alignPairs.Get(alignPairIndex).second;
			
			const CompactAlignment& alignment1 = alignments1[alignIndex1];
			const CompactAlignment& alignment2 = alignments2[alignIndex2];
			
			DebugCheck(alignment1.readID.fragmentIndex == alignment2.readID.fragmentIndex);
			
			IMatePairClusterer::MatePair matePair;
			clusterer.CreateMatePair(alignPairIndex, alignment1.region, alignment1.refStrand.strand,
									 alignment2.region, alignment2.refStrand.strand, matePair);
			
			matePairs.push_back(matePair);
			
			Region breakRegion1;
			Region breakRegion2;
			clusterer.CalcBreakRegion(alignment1.region, alignment1.refStrand.strand, breakRegion1);
			clusterer.CalcBreakRegion(alignment2.region, alignment2.refStrand.strand, breakRegion2);
			
			breakRegions1.push_back(breakRegion1);
			breakRegions2.push_back(breakRegion2);
			
			fragmentIndices.push_back(alignment1.readID.fragmentIndex);
		}
		
		cout << matePairs.size() << endl;
		
		IntegerTable matePairClusters;
		clusterer.DoClustering(matePairs, matePairClusters);
		
		for (int clusterIndex = 0; clusterIndex < matePairClusters.size(); clusterIndex++)
		{
			const IntegerVec& cluster = matePairClusters[clusterIndex];
			
			if (cluster.size() < minClusterSize)
			{
				continue;
			}
			
			unordered_set<int> clusterFragmentIndices;
			for (IntegerVecConstIter elementIter = cluster.begin(); elementIter != cluster.end(); elementIter++)
			{
				IntegerPair alignPair = alignPairs.Get(*elementIter);
				int fragmentIndex = fragmentIndices[*elementIter];
				
				// For fragments with multiple alignments supporting the same cluster, randomly select one alignment
				if (clusterFragmentIndices.find(fragmentIndex) != clusterFragmentIndices.end())
				{
					continue;
				}
				
				clusterFragmentIndices.insert(fragmentIndex);
				
				int alignIndex1 = alignPair.first;
				int alignIndex2 = alignPair.second;
				
				const CompactAlignment& alignment1 = alignments1[alignIndex1];
				const CompactAlignment& alignment2 = alignments2[alignIndex2];
				
				OutputClusterMember(clustersFile, clusterID, 0, alignment1, referenceNames);
				OutputClusterMember(clustersFile, clusterID, 1, alignment2, referenceNames);
			}
			
			clusterID++;
		}
	}
		
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

