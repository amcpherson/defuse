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
#include "MatePairCombinatorial.h"

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

struct IDAlignment
{
	int id;
	CompactAlignment alignment;
};

typedef vector<IDAlignment> IDAlignmentVec;
typedef vector<IDAlignment>::iterator IDAlignmentVecIter;
typedef vector<IDAlignment>::const_iterator IDAlignmentVecConstIter;

bool IDAlignmentFragmentLessThan(const IDAlignment& a1, const IDAlignment& a2)
{
	return (a1.alignment.readID.fragmentIndex < a2.alignment.readID.fragmentIndex);
}

bool IDAlignmentIDLessThan(const IDAlignment& a1, const IDAlignment& a2)
{
	return (a1.id < a2.id);
}

void ReadAlignments(const string& bamFilename, StringVec& referenceNames, IDAlignmentVec& alignments)
{
	samfile_t* inBamFile = samopen(bamFilename.c_str(), "rb", 0);
	
	referenceNames.clear();
	for (int targetIndex = 0; targetIndex < inBamFile->header->n_targets; targetIndex++)
	{
		referenceNames.push_back(string(inBamFile->header->target_name[targetIndex]));
	}

	bam1_t bamEntry;
	memset(&bamEntry, 0, sizeof(bam1_t));
	while (bam_read1(inBamFile->x.bam, &bamEntry) > 0)
	{
		// Split qname into id and end
		// Fragment index encoded in fragment name
		string qname = string((char*)bamEntry.data);
		vector<string> qnameFields;
		split(qnameFields, qname, is_any_of("/"));
		
		IDAlignment idAlignment;
		idAlignment.id = alignments.size();
		idAlignment.alignment.readID.fragmentIndex = lexical_cast<int>(qnameFields[0]);
		idAlignment.alignment.readID.readEnd = (qnameFields[1] == "1") ? 0 : 1;
		idAlignment.alignment.refStrand.referenceIndex = bamEntry.core.tid;
		idAlignment.alignment.refStrand.strand = (bamEntry.core.flag & 0x10) ? MinusStrand : PlusStrand;
		idAlignment.alignment.region.start = bamEntry.core.pos + 1;
		idAlignment.alignment.region.end = bamEntry.core.pos + bamEntry.core.l_qseq;
		
		alignments.push_back(idAlignment);
		
		memset(&bamEntry, 0, sizeof(bam1_t));
	}
	
	samclose(inBamFile);
}

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

typedef unordered_map<RefBin,int> RefBinMap;
typedef unordered_map<RefBin,int>::const_iterator RefBinMapConstIter;

void GetConcordanceBins(const Region& region, int binSpacing, IntegerVec& bins)
{
	int startBin = region.start / binSpacing;
	int endBin = region.end / binSpacing;
	
	if (startBin == endBin)
	{
		int bin1 = startBin;
		int bin2 = startBin + 1;
		
		bins.push_back(bin1);
		bins.push_back(bin2);
	}
	else
	{
		DebugCheck(startBin + 1 == endBin);
		
		int bin = startBin + 1;
		
		bins.push_back(bin);
	}	
}

int main(int argc, char* argv[])
{
	string discordantBamFilename;
	string geneTranscriptsFilename;
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
		TCLAP::ValueArg<string> geneTranscriptsFilenameArg("g","genetran","Gene Transcripts Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Output Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","fragmentmean","Fragment Length Mean",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","fragmentstddev","Fragment Length Standard Deviation",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> precisionArg("a","precision","Precision",true,-1,"double",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> resolutionArg("r","resolution","Resolution",false,1,"integer",cmd);
		cmd.parse(argc,argv);

		discordantBamFilename = discordantBamFilenameArg.getValue();
		geneTranscriptsFilename = geneTranscriptsFilenameArg.getValue();
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
	
	const int concordanceBinSpacing = (int)(fragmentLengthMean + 10 * fragmentLengthStdDev);
	
	cout << "Reading discordant alignments" << endl;
	
	StringVec referenceNames;
	IDAlignmentVec alignments;
	ReadAlignments(discordantBamFilename, referenceNames, alignments);

	cout << "Read " << alignments.size() << " alignments" << endl;
	
	cout << "Reading gene transcripts" << endl;
	
	ExonRegions exonRegions;
	
	ifstream geneTranscriptsFile(geneTranscriptsFilename.c_str());
	if (!geneTranscriptsFile)
	{
		cout << "Error: unable to gene transcripts" << endl;		
		exit(1);
	}
	
	exonRegions.ReadGeneTranscripts(geneTranscriptsFile);
		
	unordered_map<string,int> referenceNameIndices;
	for (int referenceIndex = 0; referenceIndex < referenceNames.size(); referenceIndex++)
	{
		referenceNameIndices[referenceNames[referenceIndex]] = referenceIndex;
	}

	int numReferenceNames = referenceNames.size();
	IntegerVec transcriptGene(numReferenceNames);
	for (int referenceIndex = 0; referenceIndex < numReferenceNames; referenceIndex++)
	{
		const string& referenceName = referenceNames[referenceIndex];
		if (exonRegions.IsTranscript(referenceName))
		{
			const string& geneName = exonRegions.GetTranscriptGeneName(referenceName);
			
			if (referenceNameIndices.find(geneName) == referenceNameIndices.end())
			{
				referenceNameIndices[geneName] = referenceNameIndices.size();
			}
			
			transcriptGene[referenceIndex] = referenceNameIndices[geneName];
		}
		else
		{
			transcriptGene[referenceIndex] = referenceIndex;
		}
	}
	
	cout << "Sorting alignments by fragment index" << endl;

	sort(alignments.begin(), alignments.end(), IDAlignmentFragmentLessThan);
	
	cout << "Finding pairs of reference sequences connected by pairs of alignments" << endl;
	
	unordered_map<IntegerPair,IntegerPairVec> refPairAlignPairMap;
	IntegerVec currentEndRefIndex[2];
	RefBinMap currentEndRefBinSelAln[2];
		
	int currentFragmentIndex = -1;
	for (IDAlignmentVecConstIter alignIter = alignments.begin(); alignIter != alignments.end(); alignIter++)
	{
		int alignIndex = alignIter->id;
		const CompactAlignment& alignment = alignIter->alignment;

		if ((alignIter + 1) == alignments.end() || alignment.readID.fragmentIndex != currentFragmentIndex)
		{
			// Skip concordant alignments
			sort(currentEndRefIndex[0].begin(), currentEndRefIndex[0].end());
			sort(currentEndRefIndex[1].begin(), currentEndRefIndex[1].end());			
			IntegerVecConstIter intersectEndIter = set_intersection(currentEndRefIndex[0].begin(), currentEndRefIndex[0].end(), currentEndRefIndex[1].begin(), currentEndRefIndex[1].end(), currentEndRefIndex[0].begin());
			if (intersectEndIter == currentEndRefIndex[0].begin())
			{
				// Create list of unique align index pairs
				unordered_set<IntegerPair> uniqueAlignIndexPairs;
				for (RefBinMapConstIter align1Iter = currentEndRefBinSelAln[0].begin(); align1Iter != currentEndRefBinSelAln[0].end(); align1Iter++)
				{
					for (RefBinMapConstIter align2Iter = currentEndRefBinSelAln[1].begin(); align2Iter != currentEndRefBinSelAln[1].end(); align2Iter++)
					{
						int refID1 = align1Iter->first.ref;
						int refID2 = align2Iter->first.ref;
						
						int alignIndex1 = align1Iter->second;
						int alignIndex2 = align2Iter->second;
						
						if (!uniqueAlignIndexPairs.insert(IntegerPair(alignIndex1,alignIndex2)).second)
						{
							continue;
						}
						
						if (refID1 > refID2)
						{
							swap(refID1,refID2);
							swap(alignIndex1,alignIndex2);
						}
						
						refPairAlignPairMap[IntegerPair(refID1,refID2)].push_back(IntegerPair(alignIndex1,alignIndex2));						
					}
				}
			}
			
			for (int end = 0; end < 2; end++)
			{
				currentEndRefIndex[end].clear();
				currentEndRefBinSelAln[end].clear();
			}
		}
		
		if (alignIndex < alignments.size())
		{
			// Add gene reference index for concordance check
			int geneReferenceIndex = transcriptGene[alignment.refStrand.referenceIndex];
			currentEndRefIndex[alignment.readID.readEnd].push_back(geneReferenceIndex);
			
			// Each read end can only have one alignment to each binned region
			IntegerVec concordanceBins;
			GetConcordanceBins(alignment.region, concordanceBinSpacing, concordanceBins);
			
			for (IntegerVecConstIter concBinIter = concordanceBins.begin(); concBinIter != concordanceBins.end(); concBinIter++)
			{
				RefBin refStrandBin(alignment.refStrand.id,*concBinIter);				
				currentEndRefBinSelAln[alignment.readID.readEnd][refStrandBin] = alignIndex;
			}			
		}
		
		currentFragmentIndex = alignment.readID.fragmentIndex;
	}
	
	// Initialize clusterer
	cout << "Initializing clusterer" << endl;
	MatePairCombinatorial clusterer;
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
	
	cout << "Sorting alignments by alignment index" << endl;

	sort(alignments.begin(), alignments.end(), IDAlignmentIDLessThan);
	
	cout << "Creating clusters" << endl;
	
	int clusterID = 0;
	for (unordered_map<IntegerPair,IntegerPairVec>::const_iterator refPairAlignPairIter = refPairAlignPairMap.begin(); refPairAlignPairIter != refPairAlignPairMap.end(); refPairAlignPairIter++)
	{
		const IntegerPair& refPair = refPairAlignPairIter->first;
		const IntegerPairVec& alignPairs = refPairAlignPairIter->second;
		
		if (alignPairs.size() < minClusterSize)
		{
			continue;
		}
		
		RefStrand refStrand1;
		RefStrand refStrand2;

		refStrand1.id = refPair.first;
		refStrand2.id = refPair.second;
		
		string refName1 = referenceNames[refStrand1.referenceIndex];
		string refName2 = referenceNames[refStrand2.referenceIndex];
		
		string strand1 = (refStrand1.strand == PlusStrand) ? "+" : "-";
		string strand2 = (refStrand2.strand == PlusStrand) ? "+" : "-";
		
		bool refIsGene1 = ("ENSG" == refName1.substr(0,4));
		bool refIsGene2 = ("ENSG" == refName2.substr(0,4));

		if (!refIsGene1 && !refIsGene2)
		{
			continue;
		}

		IMatePairClusterer::MatePairVec matePairs;
		RegionVec breakRegions1;
		RegionVec breakRegions2;
		IntegerVec fragmentIndices;
		
		for (int alignPairIndex = 0; alignPairIndex < alignPairs.size(); alignPairIndex++)
		{
			int alignIndex1 = alignPairs[alignPairIndex].first;
			int alignIndex2 = alignPairs[alignPairIndex].second;
			
			const CompactAlignment& alignment1 = alignments[alignIndex1].alignment;
			const CompactAlignment& alignment2 = alignments[alignIndex2].alignment;
			
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
				
				clustersFile << clusterID << "\t" << fragmentIndex << "\t" << alignIndex1 << "\t" << alignIndex2 << endl;
			}
			
			clusterID++;
		}
	}
	
	clustersFile.close();
	
	cout << "Created " << clusterID << " clusters" << endl;
}

