/*
 *  clustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "AlignmentStream.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
#include "IMatePairClusterer.h"
#include "Indexer.h"
#include "MatePairEM.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


struct RefBinPacked
{
	RefBinPacked() {}
	RefBinPacked(const RefStrand& refStrand, int bin) : referenceIndex(refStrand.referenceIndex), strand(refStrand.strand), bin(bin)
	{
		if (refStrand.referenceIndex >= (1<<18))
		{
			cout << refStrand.referenceIndex << endl;
			cout << (1<<18) << endl;
			cerr << "Packing failed, too many reference sequences" << endl;
			exit(1);
		}

		if (bin >= (1<<13))
		{
			cout << bin << endl;
			cout << (1<<13) << endl;
			cerr << "Packing failed, chromosome too large" << endl;
			exit(1);
		}
	}
	
	union
	{
		struct
		{
			unsigned int referenceIndex : 18;
			unsigned int strand : 1;
			unsigned int bin : 13;
		};
		
		unsigned int id;
	};
	
	bool operator<(const RefBinPacked& other) const
	{
		return id < other.id;
	}
};

inline bool operator==(const RefBinPacked& rb1, const RefBinPacked& rb2)
{
	return rb1.id == rb2.id;
}

inline size_t hash_value(const RefBinPacked& rb)
{
	size_t seed = 0;
	hash_combine(seed, rb.id);
	return seed;
}

struct AlignmentPacked
{
	ReadID readID;
	unsigned short int relativeStart;
	unsigned short int relativeEnd;
};

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

typedef vector<AlignmentPacked> AlignmentPackedVec;
typedef vector<AlignmentPacked>::const_iterator AlignmentPackedVecConstIter;

typedef unordered_set<RefBin> RefBinSet;
typedef unordered_set<RefBin>::const_iterator RefBinSetConstIter;

typedef unordered_map<RefBin,IntegerVec> RefBinVecMap;
typedef unordered_map<RefBin,IntegerVec>::const_iterator RefBinVecMapConstIter;

typedef unordered_map<RefBinPacked,AlignmentPackedVec> RefBinAlnMap;
typedef unordered_map<RefBinPacked,AlignmentPackedVec>::const_iterator RefBinAlnMapConstIter;

typedef pair<RefBinPacked,RefBinPacked> RefBinPackedPair;
typedef pair<AlignmentPackedVec,AlignmentPackedVec> AlignmentPackedVecPair;
typedef unordered_map<RefBinPackedPair,AlignmentPackedVecPair> RefBinPairMap;
typedef unordered_map<RefBinPackedPair,AlignmentPackedVecPair>::const_iterator RefBinPairMapConstIter;

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
	
	int CalcBinRelativePos(int bin, int pos) const
	{
		return pos - bin * mBinLength + mBinLength / 2;
	}

	int CalcAbsolutePos(int bin, int pos) const
	{
		return pos + bin * mBinLength - mBinLength / 2;
	}
	
private:
	int mBinLength;
	int mRegionExtend;
};

void PackAlignment(const Binning& binning, int bin, const CompactAlignment& alignment, AlignmentPacked& packed)
{
	int relativeStart = binning.CalcBinRelativePos(bin, alignment.region.start);
	int relativeEnd = binning.CalcBinRelativePos(bin, alignment.region.end);
	
	DebugCheck(relativeStart >= 0);
	DebugCheck(relativeEnd >= 0);
	DebugCheck(relativeStart < (1<<16));
	DebugCheck(relativeEnd < (1 << 16));

	packed.readID = alignment.readID;
	packed.relativeStart = relativeStart;
	packed.relativeEnd = relativeEnd;
}

void UnpackAlignment(const Binning& binning, const RefBinPacked& refBin, const AlignmentPacked& packed, CompactAlignment& alignment)
{
	alignment.refStrand.referenceIndex = refBin.referenceIndex;
	alignment.refStrand.strand = refBin.strand;
	alignment.readID = packed.readID;
	alignment.region.start = binning.CalcAbsolutePos(refBin.bin, packed.relativeStart);
	alignment.region.end = binning.CalcAbsolutePos(refBin.bin, packed.relativeEnd);
}

void UnpackAlignments(const Binning& binning, const RefBinPacked& refBin, const AlignmentPackedVec& packed, CompAlignVec& alignments)
{
	alignments.resize(packed.size());
	for (int packedIndex = 0; packedIndex < packed.size(); packedIndex++)
	{
		UnpackAlignment(binning, refBin, packed[packedIndex], alignments[packedIndex]);
	}
}

bool CheckConcordant(const CompAlignVec& alignments, int minFusionRange)
{
	CompAlignVec filtered;
	
	// Check for concordance
	RefBinSet concordanceBins[2];
	Binning concordanceBinning(minFusionRange, minFusionRange);
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
	RefBinAlnMap binnedEndAlignments[2];
	for (CompAlignVecConstIter alignIter = alignments.begin(); alignIter != alignments.end(); alignIter++)
	{
		const CompactAlignment& alignment = *alignIter;
		
		IntegerVec clusterBins;
		binning.GetBins(alignment.region, clusterBins);
		for (IntegerVecConstIter clusterBinIter = clusterBins.begin(); clusterBinIter != clusterBins.end(); clusterBinIter++)
		{
			RefBinPacked refBin(alignment.refStrand,*clusterBinIter);
			AlignmentPacked alignmentPacked;
			PackAlignment(binning, *clusterBinIter, alignment, alignmentPacked);

			binnedEndAlignments[alignment.readID.readEnd][refBin].push_back(alignmentPacked);
		}
	}
	
	// Create binned alignments
	for (RefBinAlnMapConstIter refBin1Iter = binnedEndAlignments[0].begin(); refBin1Iter != binnedEndAlignments[0].end(); refBin1Iter++)
	{
		for (RefBinAlnMapConstIter refBin2Iter = binnedEndAlignments[1].begin(); refBin2Iter != binnedEndAlignments[1].end(); refBin2Iter++)
		{
			RefBinPacked refBin1 = refBin1Iter->first;
			RefBinPacked refBin2 = refBin2Iter->first;
			
			if (refBin1 < refBin2)
			{
				RefBinPackedPair refBinPair(refBin1,refBin2);
				
				binPairs[refBinPair].first.insert(binPairs[refBinPair].first.end(), refBin1Iter->second.begin(), refBin1Iter->second.end());
				binPairs[refBinPair].second.insert(binPairs[refBinPair].second.end(), refBin2Iter->second.begin(), refBin2Iter->second.end());
			}
			else
			{
				RefBinPackedPair refBinPair(refBin2,refBin1);
				
				binPairs[refBinPair].first.insert(binPairs[refBinPair].first.end(), refBin2Iter->second.begin(), refBin2Iter->second.end());
				binPairs[refBinPair].second.insert(binPairs[refBinPair].second.end(), refBin1Iter->second.begin(), refBin1Iter->second.end());
			}
		}
	}	
}

void GetFragmentAlignments(const CompAlignVec& alignments, IntegerVecMap& fragments)
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
	string alignmentsFilename;
	string clustersFilename;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	double precision;
	int minClusterSize;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> alignmentsFilenameArg("a","align","Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Output Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","fragmentmean","Fragment Length Mean",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","fragmentstddev","Fragment Length Standard Deviation",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> precisionArg("p","precision","Precision",true,-1,"double",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		cmd.parse(argc,argv);

		alignmentsFilename = alignmentsFilenameArg.getValue();
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
	const int binLength = 1<<15;
		
	cout << "Finding pairs of reference sequences connected by pairs of alignments" << endl;
	
	AlignmentStream* alignmentStream = new CompactAlignmentStream(alignmentsFilename);
	FragmentAlignmentStream fragmentAlignmentStream(alignmentStream);
	
	NameIndex refNameIndex;
	
	RefBinPairMap binPairs;
	Binning binning(binLength, minFusionRange);
	RawAlignmentVec rawAlignments;
	while (fragmentAlignmentStream.GetNextAlignments(rawAlignments))
	{
		CompAlignVec alignments;
		for (int alignIndex = 0; alignIndex < rawAlignments.size(); alignIndex++)
		{
			alignments.push_back(CompactAlignment());
			alignments.back().readID.fragmentIndex = lexical_cast<int>(rawAlignments[alignIndex].fragment);
			alignments.back().readID.readEnd = rawAlignments[alignIndex].readEnd;
			alignments.back().refStrand.referenceIndex = refNameIndex.Index(rawAlignments[alignIndex].reference);
			alignments.back().refStrand.strand = rawAlignments[alignIndex].strand;
			alignments.back().region = rawAlignments[alignIndex].region;
		}
		
		if (CheckConcordant(alignments, minFusionRange))
		{
			continue;
		}
		
		AddBinPairs(alignments, binning, binPairs);
	}
	
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
		RefBinPacked refBin1 = refBinIter->first.first;
		RefBinPacked refBin2 = refBinIter->first.second;
		
		if (refBinIter->second.first.size() < minClusterSize || refBinIter->second.second.size() < minClusterSize)
		{
			continue;
		}
				
		CompAlignVec alignments1;
		CompAlignVec alignments2;
		
		UnpackAlignments(binning, refBin1, refBinIter->second.first, alignments1);
		UnpackAlignments(binning, refBin2, refBinIter->second.second, alignments2);
		
		IntegerVecMap fragments1;
		IntegerVecMap fragments2;
		
		GetFragmentAlignments(alignments1, fragments1);
		GetFragmentAlignments(alignments2, fragments2);
		
		FilterUnmatched(fragments1, fragments2);
		FilterUnmatched(fragments2, fragments1);
		
		FilterOverlapping(fragments1, alignments1, minFusionRange);
		FilterOverlapping(fragments2, alignments2, minFusionRange);
		
		if (fragments1.size() < minClusterSize || fragments2.size() < minClusterSize)
		{
			continue;
		}
		
		PairIndex alignPairs;
		GetAlignPairs(fragments1, fragments2, alignPairs);
		
		IMatePairClusterer::MatePairVec matePairs;
		RegionVec breakRegions1;
		RegionVec breakRegions2;
		IntegerVec fragmentIndices;
		
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
				
				OutputClusterMember(clustersFile, clusterID, 0, alignment1, refNameIndex.Get());
				OutputClusterMember(clustersFile, clusterID, 1, alignment2, refNameIndex.Get());
			}

			clusterID++;
		}
	}
		
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

