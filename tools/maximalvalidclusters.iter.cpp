/*
 *  clustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Alignments.h"
#include "DebugCheck.h"
#include "ExonRegions.h"
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

class FragmentAlignmentsIterator
{
public:
	FragmentAlignmentsIterator() : mBamFile(0), mGood(false)
	{}
	
	~FragmentAlignmentsIterator()
	{
		Close();
	}
	
	void Open(const string& bamFilename, StringVec& referenceNames)
	{
		mBamFile = samopen(bamFilename.c_str(), "rb", 0);
		
		referenceNames.clear();
		for (int targetIndex = 0; targetIndex < mBamFile->header->n_targets; targetIndex++)
		{
			referenceNames.push_back(string(mBamFile->header->target_name[targetIndex]));
		}
		
		bam1_t bamEntry;
		memset(&bamEntry, 0, sizeof(bam1_t));
		mGood = bam_read1(mBamFile->x.bam, &bamEntry) > 0;
		
		if (mGood)
		{
			ConvertBamAlignment(bamEntry, mCurrentAlignment);		
		}
		
		bam1_t bamEntry;
		memset(&bamEntry, 0, sizeof(bam1_t));
		while (bam_read1(mBamFile->x.bam, &bamEntry) > 0)
		{
			// Split qname into id and end
			// Fragment index encoded in fragment name
			string qname = string((char*)bamEntry.data);
			vector<string> qnameFields;
			split(qnameFields, qname, is_any_of("/"));
			
			int alignmentIndex = alignments.size();
			
			CompactAlignment alignment;
			alignment.readID.fragmentIndex = lexical_cast<int>(qnameFields[0]);
			alignment.readID.readEnd = (qnameFields[1] == "1") ? 0 : 1;
			alignment.refStrand.referenceIndex = bamEntry.core.tid;
			alignment.refStrand.strand = (bamEntry.core.flag & 0x10) ? MinusStrand : PlusStrand;
			alignment.region.start = bamEntry.core.pos + 1;
			alignment.region.end = bamEntry.core.pos + bamEntry.core.l_qseq;
			
			alignments.push_back(alignment);
			
			fragments.resize(max((int)fragments.size(), (int)alignment.readID.fragmentIndex + 1));
			fragments[alignment.readID.fragmentIndex].push_back(alignmentIndex);
			
			memset(&bamEntry, 0, sizeof(bam1_t));
		}
		
		samclose(mBamFile);		
	}

	bool GetNext(CompAlignVec& alignments)
	{
		if (!mGood)
		{
			return false;
		}
		
		alignments.push_back(mCurrentAlignment);
		
		while (true)
		{
			bam1_t bamEntry;
			memset(&bamEntry, 0, sizeof(bam1_t));
			mGood = bam_read1(mBamFile->x.bam, &bamEntry) > 0;
			if (!mGood)
			{
				break;
			}
			
			CompactAlignment nextAlignment;
			ConvertBamAlignment(bamEntry, nextAlignment);
			
			if (nextAlignment.readID.fragmentIndex == mCurrentAlignment.readID.fragmentIndex)
			{
				alignments.push_back(nextAlignment);
			}
			else
			{
				mCurrentAlignment = nextAlignment;
				break;
			}
		}
		
		return true;
	}
	
	void Close(const string& bamFilename, StringVec& referenceNames, CompAlignVec& alignments, IntegerTable& fragments)
	{
		samclose(mBamFile);		
	}
	
private:
	
	void ConvertBamAlignment(const bam1_t& bamEntry, CompactAlignment& alignment)
	{
		// Split qname into id and end
		// Fragment index encoded in fragment name
		string qname = string((char*)bamEntry.data);
		vector<string> qnameFields;
		split(qnameFields, qname, is_any_of("/"));
		
		alignment.readID.fragmentIndex = lexical_cast<int>(qnameFields[0]);
		alignment.readID.readEnd = (qnameFields[1] == "1") ? 0 : 1;
		alignment.refStrand.referenceIndex = bamEntry.core.tid;
		alignment.refStrand.strand = (bamEntry.core.flag & 0x10) ? MinusStrand : PlusStrand;
		alignment.region.start = bamEntry.core.pos + 1;
		alignment.region.end = bamEntry.core.pos + bamEntry.core.l_qseq;
	}
	
	samfile_t* mBamFile;
	bool mGood;
	CompactAlignment mCurrentAlignment;
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

typedef unordered_set<RefBin> RefBinSet;
typedef unordered_set<RefBin>::const_iterator RefBinSetConstIter;

typedef unordered_map<RefBin,IntegerVec> RefBinVecMap;
typedef unordered_map<RefBin,IntegerVec>::const_iterator RefBinVecMapConstIter;

typedef pair<RefBin,RefBin> RefBinPair;
typedef unordered_map<RefBinPair,IntegerPairVec> RefBinPairVecMap;
typedef unordered_map<RefBinPair,IntegerPairVec>::const_iterator RefBinPairVecMapConstIter;

typedef pair<RefBin,int> RefBinAlign;
typedef vector<RefBinAlign> RefBinAlignVec;
typedef vector<RefBinAlign>::const_iterator RefBinAlignVecConstIter;

void GetBins(const Region& region, IntegerVec& bins, int binLength, int regionRange)
{
	int startBin = (region.start - regionRange) / binLength;
	int endBin = (region.end + regionRange) / binLength;
	
	for (int bin = startBin; bin <= endBin; bin++)
	{
		bins.push_back(bin);
	}
}

int main(int argc, char* argv[])
{
	string discordantBamFilename;
	string clustersFilename;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	double precision;
	int minClusterSize;
	int resolution;

	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> discordantBamFilenameArg("d","discordant","Discordant Alignments Bam Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Output Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","fragmentmean","Fragment Length Mean",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","fragmentstddev","Fragment Length Standard Deviation",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> precisionArg("a","precision","Precision",true,-1,"double",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> resolutionArg("r","resolution","Resolution",false,1,"integer",cmd);
		cmd.parse(argc,argv);

		discordantBamFilename = discordantBamFilenameArg.getValue();
		clustersFilename = clustersFilenameArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		precision = precisionArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
		resolution = resolutionArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	const int fusionRangeMin = (int)(fragmentLengthMean + 10 * fragmentLengthStdDev);
	const int binLength = 100000;
	
	cout << "Reading discordant alignments" << endl;
	
	StringVec referenceNames;
	CompAlignVec alignments;
	IntegerTable fragments;
	ReadAlignments(discordantBamFilename, referenceNames, alignments, fragments);

	cout << "Read " << alignments.size() << " alignments" << endl;
	
	cout << "Finding pairs of reference sequences connected by pairs of alignments" << endl;
	
	RefBinPairVecMap binnedAlignments;
	for (int fragmentIndex = 0; fragmentIndex < fragments.size(); fragmentIndex++)
	{
		// Check for concordance
		RefBinSet concordanceBins[2];
		for (IntegerVecConstIter alignIter = fragments[fragmentIndex].begin(); alignIter != fragments[fragmentIndex].end(); alignIter++)
		{
			int alignIndex = *alignIter;
			const CompactAlignment& alignment = alignments[alignIndex];
			
			// Bin end alignments according to bins with length bounded by minimum fusion distance
			IntegerVec rangeBins;
			GetBins(alignment.region, rangeBins, fusionRangeMin, 0);
			
			// Insert alignments for concordance check
			for (IntegerVecConstIter rangeBinIter = rangeBins.begin(); rangeBinIter != rangeBins.end(); rangeBinIter++)
			{
				RefBin concordanceBin(alignment.refStrand.referenceIndex,*rangeBinIter);				
				concordanceBins[alignment.readID.readEnd].insert(concordanceBin);
			}
		}
		
		// Classify as concordant if both ends aligning to same region		
		bool concordant = false;
		for (RefBinSetConstIter concordanceBinIter = concordanceBins[0].begin(); concordanceBinIter != concordanceBins[0].end(); concordanceBinIter++)
		{
			if (concordanceBins[1].find(*concordanceBinIter) != concordanceBins[1].end())
			{
				concordant = true;
			}
		}
		if (concordant)
		{
			continue;
		}
		
		// Select a set of non overlapping end alignments
		RefBinSet alignmentBins[2];
		IntegerVec nonOverlapping;
		for (IntegerVecConstIter alignIter = fragments[fragmentIndex].begin(); alignIter != fragments[fragmentIndex].end(); alignIter++)
		{
			int alignIndex = *alignIter;
			const CompactAlignment& alignment = alignments[alignIndex];
			
			IntegerVec rangeBins;
			GetBins(alignment.region, rangeBins, fusionRangeMin, 0);
			
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
				nonOverlapping.push_back(alignIndex);
			}
		}
			
		// Create binned end alignments
		RefBinAlignVec binnedEndAlignments[2];
		for (IntegerVecConstIter alignIter = nonOverlapping.begin(); alignIter != nonOverlapping.end(); alignIter++)
		{
			int alignIndex = *alignIter;
			const CompactAlignment& alignment = alignments[alignIndex];
			
			IntegerVec clusterBins;
			GetBins(alignment.region, clusterBins, binLength, fusionRangeMin);
			for (IntegerVecConstIter clusterBinIter = clusterBins.begin(); clusterBinIter != clusterBins.end(); clusterBinIter++)
			{
				RefBinAlign refBinAlign(RefBin(alignment.refStrand.id,*clusterBinIter),alignIndex);
				binnedEndAlignments[alignment.readID.readEnd].push_back(refBinAlign);
			}
		}
		
		// Create binned alignments
		for (RefBinAlignVecConstIter refBinAlign1Iter = binnedEndAlignments[0].begin(); refBinAlign1Iter != binnedEndAlignments[0].end(); refBinAlign1Iter++)
		{
			for (RefBinAlignVecConstIter refBinAlign2Iter = binnedEndAlignments[1].begin(); refBinAlign2Iter != binnedEndAlignments[1].end(); refBinAlign2Iter++)
			{
				RefBin refBin1 = refBinAlign1Iter->first;
				RefBin refBin2 = refBinAlign2Iter->first;
				
				int alignIndex1 = refBinAlign1Iter->second;
				int alignIndex2 = refBinAlign2Iter->second;
				
				if (refBin2 < refBin1)
				{
					swap(refBin1,refBin2);
					swap(alignIndex1,alignIndex2);
				}
				
				binnedAlignments[RefBinPair(refBin1,refBin2)].push_back(IntegerPair(alignIndex1,alignIndex2));						
			}
		}
	}
	
	// Initialize clusterer
	cout << "Initializing clusterer" << endl;
	MatePairEM clusterer;
	IMatePairClusterer::Parameters clustererParams;
	clustererParams.fragmentMean = fragmentLengthMean;
	clustererParams.fragmentStdDev = fragmentLengthStdDev;
	clustererParams.precision = precision;
	clustererParams.minClusterSize = minClusterSize;
	clustererParams.resolution = resolution;
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
	for (RefBinPairVecMapConstIter alignPairIter = binnedAlignments.begin(); alignPairIter != binnedAlignments.end(); alignPairIter++)
	{
		const RefBinPair& refBinPair = alignPairIter->first;
		const IntegerPairVec& alignPairs = alignPairIter->second;
		
		if (alignPairs.size() < minClusterSize)
		{
			continue;
		}
		
		RefStrand refStrand1;
		RefStrand refStrand2;

		refStrand1.id = refBinPair.first.ref;
		refStrand2.id = refBinPair.second.ref;
		
		string refName1 = referenceNames[refStrand1.referenceIndex];
		string refName2 = referenceNames[refStrand2.referenceIndex];
		
		string strand1 = (refStrand1.strand == PlusStrand) ? "+" : "-";
		string strand2 = (refStrand2.strand == PlusStrand) ? "+" : "-";
		
		IMatePairClusterer::MatePairVec matePairs;
		RegionVec breakRegions1;
		RegionVec breakRegions2;
		IntegerVec fragmentIndices;
		
		for (int alignPairIndex = 0; alignPairIndex < alignPairs.size(); alignPairIndex++)
		{
			int alignIndex1 = alignPairs[alignPairIndex].first;
			int alignIndex2 = alignPairs[alignPairIndex].second;
			
			const CompactAlignment& alignment1 = alignments[alignIndex1];
			const CompactAlignment& alignment2 = alignments[alignIndex2];
			
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
				IntegerPair alignPair = alignPairs[*elementIter];
				int fragmentIndex = fragmentIndices[*elementIter];
				
				// For fragments with multiple alignments supporting the same cluster, randomly select one alignment
				if (clusterFragmentIndices.find(fragmentIndex) != clusterFragmentIndices.end())
				{
					continue;
				}
				
				clusterFragmentIndices.insert(fragmentIndex);
				
				int alignIndex1 = alignPair.first;
				int alignIndex2 = alignPair.second;
				
				clustersFile << clusterID << "\t" << "0" << "\t" << fragmentIndex << "\t" << alignIndex1 << endl;
				clustersFile << clusterID << "\t" << "1" << "\t" << fragmentIndex << "\t" << alignIndex2 << endl;
			}
			
			clusterID++;
		}
	}
	
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

