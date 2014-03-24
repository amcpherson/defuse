/*
 *  FusionPredictions.h
 *  findbreaks
 *
 *  Created by Andrew McPherson on 28/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FUSIONPREDICTIONS_H_
#define FUSIONPREDICTIONS_H_

#include "Alignments.h"
#include "Parameters.h"
#include "ReadStream.h"
#include "ExonRegions.h"
#include "Reference.h"

#include <vector>

using namespace std;

class IAlignmentStream;
 
class FusionPredictions
{
public:
	FusionPredictions(const Alignments& alignments, const Parameters& parameters, const IntegerPairTable& clusters, const ExonRegions& exonRegions);
	
	void FilterLowCoverage();
	void ClearClusterReads();
	void FindSpanningClusterReads();
	void FindHangingClusterReads(IAlignmentStream* alignments);
	void RetrieveReadSequences(IReadStream* reads);
	void RetrieveClusterSequences(Reference& reference);
	void RetrievePartnerSequences(Reference& reference);
	void FindBreaks();
	
	int GetNumClusters();
	
	void PrintClusters();
	
private:
	void InitializeClusterBreaks();
	bool GetClusterRefStrand(int clusterIndex, int clusterEnd, RefStrand& refStrand);
	bool CalculateClusterAlignRegion(int clusterIndex, int clusterEnd, int& alignRegionStart, int& alignRegionLength);
	bool CalculateClusterBreakRegion(int clusterIndex, int clusterEnd, int& breakRegionStart, int& breakRegionLength);

	void CalculateBreakRegion(int strand, const Region& alignRegion, const Parameters& parameters, Region& breakRegion) const;
	bool CalculateClusterBreakRegions(int clusterIndex, RefStrand& refStrand1, Region& breakRegion1, RefStrand& refStrand2, Region& breakRegion2) const;
	bool CalculateClusterAlignRegions(int clusterIndex, RefStrand& refStrand1, Region& alignRegion1, RefStrand& refStrand2, Region& alignRegion2) const;
	void CalculateAlignBreakSequence(int clusterIndex, int breakPos1, int breakPos2, string& alignBreak1, string& alignBreak2) const;
	void CalculateAlignBreakSequence(int clusterIndex, string& alignBreak1, string& alignBreak2) const;
	bool CalculateClusterCoverageStat(int clusterIndex, float& coverageStat1, float& coverageStat2) const;
	bool CalculateClusterReads(int clusterIndex, IntegerVec& readIDs1, IntegerVec& readIDs2) const;
	
	bool FindBreakDeNovo(int clusterIndex, ostream& breakOut, string& breakSequence, int& minCount, double& fragmentMean) const;
	bool FindBreakSplitReads(int clusterIndex, ostream& breakOut, string& breakSequence, int& splitReadCount, double& splitPosAvg, double& splitMinAvg, double& fragmentMean) const;
	
	const Alignments& mAlignments;
	const Parameters& mParameters;
	IntegerPairTable mClusters;
	const ExonRegions& mExonRegions;

	typedef unordered_map<IntegerPair,int> ClusterPairMap;
	typedef unordered_map<IntegerPair,int>::iterator ClusterPairMapIter;
	
	typedef unordered_map<int,string> SequenceMap;
	typedef unordered_map<int,string>::iterator SequenceMapIter;
	typedef unordered_map<int,string>::const_iterator SequenceMapConstIter;
	
	typedef unordered_map<int,int> IntegerMap;
	typedef unordered_map<int,int>::iterator IntegerMapIter;
	typedef unordered_map<int,int>::const_iterator IntegerMapConstIter;
	
	typedef vector<IntegerPair> ClusterPairVec;
	typedef vector<IntegerPair>::const_iterator ClusterPairVecConstIter;
	
	typedef vector<string> SequenceVec;
	
	typedef vector<double> DoubleVec;
	
	struct ClusterBreak
	{
		int clusterIndex;
		int end;
		RefStrand refStrand;
		Region region;
		
		typedef unsigned int ChromosomeType;
		
		unsigned int GetChromosome()
		{
			return refStrand.id;
		}
		
		const Region& GetRegion()
		{
			return region;
		}
	};
	
	typedef vector<ClusterBreak> ClusterBreakVec;
	
	const static int mNumBreakPadding;
	
	IntegerTable mClusterReads;
	IntegerTable mClusterReadRevComp;
	LocationVec mClusterSequencePos[2];
	LocationVec mPartnerSequencePos[2];
	DoubleVec mClusterFragmentMean;
	SequenceVec mClusterSequences[2];
	SequenceVec mClusterBreakSequences;
	IntegerVec mSplitReadCount;
	SequenceVec mPartnerSequences[2];
	SequenceMap mReadSequences;
	DoubleVec mSplitPosAvg;
	DoubleVec mSplitMinAvg;
	IntegerMap mReadMisMatches;
	ClusterBreakVec mClusterBreaks;
};

#endif
