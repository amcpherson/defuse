/*
 *  FusionMatePairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "FusionMatePairs.h"
#include "Alignments.h"
#include "AlignmentStream.h"
#include "DebugCheck.h"

#include <map>
#include <list>
#include <fstream>
#include <sstream>
#include <boost/unordered_set.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>

using namespace std;
using namespace boost;
using namespace boost::bimaps;


FusionMatePairs::FusionMatePairs(const Alignments& alignments, const IMatePairClusterer* clusterer, int minClusterSize)
: mAlignments(alignments), mClusterer(clusterer), mMinClusterSize(minClusterSize)
{
}

void MaximizeParsimony(const IntegerTable& clusters, IntegerTable& solution)
{
	// Find max element size
	// Check that no elements are negative
	// Initialize solution
	int maxElement = -1;
	solution.clear();
	for (IntegerTableConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		solution.push_back(IntegerVec(clusterIter->size(), 0));
		
		for (IntegerVecConstIter elementIter = clusterIter->begin(); elementIter != clusterIter->end(); elementIter++)
		{
			DebugCheck(*elementIter >= 0);
			if (*elementIter < 0)
			{
				return;
			}
			
			maxElement = max(maxElement, *elementIter);
		}
	}
	
	solution.resize(clusters.size());
	
	typedef bimap<vector_of<int>,multiset_of<int> > SizeBimap;
	
	SizeBimap clusterSizeMap;
	IntegerVec clusterSizes;
	IntegerTable elementsToClusters(maxElement + 1);
	IntegerVec assignedElements(maxElement + 1, 0);
	
	// Create a table lookup from elements to clusters
	// Create a list of cluster sizes
	for (int clusterIndex = 0; clusterIndex < clusters.size(); clusterIndex++)
	{
		const IntegerVec& cluster = clusters[clusterIndex];
		
		clusterSizes.push_back(cluster.size());
		
		for (IntegerVecConstIter elementIter = cluster.begin(); elementIter != cluster.end(); elementIter++)
		{
			int element = *elementIter;
			elementsToClusters[element].push_back(clusterIndex);
		}
		
		clusterSizeMap.push_back(SizeBimap::value_type(clusterIndex,cluster.size()));		
	}
	
	// Iterate until no more clusters found
	while (clusterSizeMap.size() > 0)
	{
		// Find the next cluster with the most support
		int nextClusterIndex = clusterSizeMap.right.rbegin()->second;
		int nextClusterSize = clusterSizes[nextClusterIndex];
		
		DebugCheck(clusterSizeMap.right.rbegin()->first == nextClusterSize);
		
		// Remaining clusters are empty
		if (nextClusterSize == 0)
		{
			break;
		}
				
		// Iterate through alignment pairs of this cluster
		for (int elementIndex = 0; elementIndex < clusters[nextClusterIndex].size(); elementIndex++)
		{
			int element = clusters[nextClusterIndex][elementIndex];
			
			// Add the mate pair to the new cluster if it hasnt been assigned
			if (assignedElements[element] == 0)
			{
				// Assign to solution cluster
				solution[nextClusterIndex][elementIndex] = 1;
				
				// Flag the mate pair as assigned
				assignedElements[element] = 1;
				
				// Subtract alignment score from all clusters that used to contain this fragment
				for (IntegerVecConstIter clusterIter = elementsToClusters[element].begin(); clusterIter != elementsToClusters[element].end(); clusterIter++)
				{
					int clusterIndex = *clusterIter;
					
					// Subtract from the number of mate pairs in the cluster
					clusterSizes[clusterIndex]--;
					
					// Subtract from the total score of the cluster
					SizeBimap::left_iterator iter = clusterSizeMap.left.begin() + clusterIndex;
					
					DebugCheck(clusterIndex == iter->first);
					
					bool successful_modify = clusterSizeMap.left.replace_data(iter, iter->second - 1);
					DebugCheck(successful_modify);
				}
			}
		}
	}
}

void FusionMatePairs::CreateMatePairClusters(IntegerPairTable& clusters)
{
	unordered_map<IntegerPair,IntegerVec> refPairAlignIndices;

	Indexer<IntegerPair> alignPairIndexer;
	
	cout << "Extracting candidate mate pair alignments" << endl;

	for (int fragmentIndex = 0; fragmentIndex < Alignments::FragmentNames.Size(); fragmentIndex++)
	{
		ReadID readEnd1;
		readEnd1.fragmentIndex = fragmentIndex;
		readEnd1.readEnd = 0;
		
		const IntegerVec& alignIndicesEnd1 = mAlignments.GetAlignments(readEnd1);
		
		ReadID readEnd2;
		readEnd2.fragmentIndex = fragmentIndex;
		readEnd2.readEnd = 1;
		
		const IntegerVec& alignIndicesEnd2 = mAlignments.GetAlignments(readEnd2);

		unordered_set<string> geneNames1;
		for (IntegerVecConstIter alignIndIter = alignIndicesEnd1.begin(); alignIndIter != alignIndicesEnd1.end(); alignIndIter++)
		{
			int alignIndex = *alignIndIter;
			const Alignment& alignment = mAlignments.GetAlignment(alignIndex);
			int referenceIndex = alignment.refStrand.referenceIndex;
			const string& referenceName = Alignments::ReferenceNames.Get(referenceIndex);
			string::size_type geneNameEnd = referenceName.find_first_of('|');
			string geneName = referenceName.substr(0, geneNameEnd);

			geneNames1.insert(geneName);
		}
		
		unordered_set<string> geneNames2;
		for (IntegerVecConstIter alignIndIter = alignIndicesEnd2.begin(); alignIndIter != alignIndicesEnd2.end(); alignIndIter++)
		{
			int alignIndex = *alignIndIter;
			const Alignment& alignment = mAlignments.GetAlignment(alignIndex);
			int referenceIndex = alignment.refStrand.referenceIndex;
			const string& referenceName = Alignments::ReferenceNames.Get(referenceIndex);
			string::size_type geneNameEnd = referenceName.find_first_of('|');
			string geneName = referenceName.substr(0, geneNameEnd);

			geneNames2.insert(geneName);
		}

		if (geneNames1.size() * geneNames2.size() > 100)
		{
			continue;
		}
				
		for (IntegerVecConstIter alignIndIter1 = alignIndicesEnd1.begin(); alignIndIter1 != alignIndicesEnd1.end(); alignIndIter1++)
		{
			int alignIndex1 = *alignIndIter1;
			const Alignment& alignment1 = mAlignments.GetAlignment(alignIndex1);
			RefStrand refStrand1;
			refStrand1.id = alignment1.refStrand.id;
			
			for (IntegerVecConstIter alignIndIter2 = alignIndicesEnd2.begin(); alignIndIter2 != alignIndicesEnd2.end(); alignIndIter2++)
			{
				int alignIndex2 = *alignIndIter2;
				const Alignment& alignment2 = mAlignments.GetAlignment(alignIndex2);
				RefStrand refStrand2;
				refStrand2.id = alignment2.refStrand.id;
				
				if (refStrand1.referenceIndex == refStrand2.referenceIndex)
				{
					continue;
				}
				
				IntegerPair refPair;
				IntegerPair alignPair;
				
				if (refStrand1.id < refStrand2.id)
				{
					refPair = IntegerPair(refStrand1.id, refStrand2.id);
					alignPair = IntegerPair(alignIndex1, alignIndex2);
				}
				else
				{
					refPair = IntegerPair(refStrand2.id, refStrand1.id);
					alignPair = IntegerPair(alignIndex2, alignIndex1);
				}
				
				int alignPairIndex = alignPairIndexer.Index(alignPair);
									
				refPairAlignIndices[refPair].push_back(alignPairIndex);
			}
		}
	}

	cout << "Clustering for each pair of connected transcripts" << endl;
	
	Clusters pairIndexClusters;
	for (unordered_map<IntegerPair,IntegerVec>::iterator refPairIter = refPairAlignIndices.begin(); refPairIter != refPairAlignIndices.end(); refPairIter++)
	{
		const IntegerPair& refStrandPair = refPairIter->first;
		const IntegerVec& alignPairIndices = refPairIter->second;
		
		if (alignPairIndices.size() < mMinClusterSize)
		{
			continue;
		}
		
		RefStrand referenceStrand1;
		RefStrand referenceStrand2;
		
		referenceStrand1.id = refStrandPair.first;
		referenceStrand2.id = refStrandPair.second;
		
		const string& partner1 = Alignments::ReferenceNames.Get(referenceStrand1.referenceIndex);
		const string& partner2 = Alignments::ReferenceNames.Get(referenceStrand2.referenceIndex);
		
		IMatePairClusterer::MatePairVec matePairs;
		
		for (IntegerVecConstIter alignPairIndexIter = alignPairIndices.begin(); alignPairIndexIter != alignPairIndices.end(); alignPairIndexIter++)
		{
			int alignPairIndex = *alignPairIndexIter;
			
			IMatePairClusterer::MatePair matePair;
			
			matePair.id = alignPairIndex;

			int alignIndex1 = alignPairIndexer.Get(alignPairIndex).first;
			int alignIndex2 = alignPairIndexer.Get(alignPairIndex).second;
			
			const Alignment& alignment1 = mAlignments.GetAlignment(alignIndex1);
			const Alignment& alignment2 = mAlignments.GetAlignment(alignIndex2);
						
			matePair.alignment1.start = (alignment1.refStrand.strand == PlusStrand) ? alignment1.region.start : -alignment1.region.end;
			matePair.alignment1.end = (alignment1.refStrand.strand == PlusStrand) ? alignment1.region.end : -alignment1.region.start;
									
			matePair.alignment2.start = (alignment2.refStrand.strand == PlusStrand) ? alignment2.region.start : -alignment2.region.end;
			matePair.alignment2.end = (alignment2.refStrand.strand == PlusStrand) ? alignment2.region.end : -alignment2.region.start;
			
			matePairs.push_back(matePair);
			
			DebugCheck(alignment1.readID.fragmentIndex == alignment2.readID.fragmentIndex);
		}
		
		mClusterer->DoClustering(matePairs, pairIndexClusters);
	}

	// Prepare fragment index clusters for maximum parsimony solver
	IntegerTable fragmentIndexClusters;
	for (int clusterIndex = 0; clusterIndex < pairIndexClusters.GetNumClusters(); clusterIndex++)
	{
		const IntegerVec& pairIndexCluster = pairIndexClusters.GetCluster(clusterIndex);
		
		fragmentIndexClusters.push_back(IntegerVec());
		for (IntegerVecConstIter alignPairIndexIter = pairIndexCluster.begin(); alignPairIndexIter != pairIndexCluster.end(); alignPairIndexIter++)
		{
			int alignPairIndex = *alignPairIndexIter;
			
			int alignIndex1 = alignPairIndexer.Get(alignPairIndex).first;
			int alignIndex2 = alignPairIndexer.Get(alignPairIndex).second;
			
			const Alignment& alignment1 = mAlignments.GetAlignment(alignIndex1);
			const Alignment& alignment2 = mAlignments.GetAlignment(alignIndex2);
			
			DebugCheck(alignment1.readID.fragmentIndex == alignment2.readID.fragmentIndex);
			
			fragmentIndexClusters.back().push_back(alignment1.readID.fragmentIndex);
		}
	}
	
	// Find maximum parsimony solution
	IntegerTable maxParsSolution;
	MaximizeParsimony(fragmentIndexClusters, maxParsSolution);
	
	// Create maximum parsimony clusters
	IntegerTable maxParsClusters;
	for (int clusterIndex = 0; clusterIndex < pairIndexClusters.GetNumClusters(); clusterIndex++)
	{
		const IntegerVec& pairIndexCluster = pairIndexClusters.GetCluster(clusterIndex);

		maxParsClusters.push_back(IntegerVec());
		for (int elementIndex = 0; elementIndex < pairIndexCluster.size(); elementIndex++)
		{
			if (maxParsSolution[clusterIndex][elementIndex] == 0)
			{
				continue;
			}
			
			maxParsClusters.back().push_back(pairIndexCluster[elementIndex]);
		}
	}
	
	// Output clusters
	for (int clusterIndex = 0; clusterIndex < maxParsClusters.size(); clusterIndex++)
	{
		const IntegerVec& pairIndexCluster = maxParsClusters[clusterIndex];
		
		if (pairIndexCluster.size() < mMinClusterSize)
		{
			continue;
		}
		
		clusters.push_back(IntegerPairVec());
		for (IntegerVecConstIter alignPairIndexIter = pairIndexCluster.begin(); alignPairIndexIter != pairIndexCluster.end(); alignPairIndexIter++)
		{
			int alignPairIndex = *alignPairIndexIter;

			int alignIndex1 = alignPairIndexer.Get(alignPairIndex).first;
			int alignIndex2 = alignPairIndexer.Get(alignPairIndex).second;

			clusters.back().push_back(IntegerPair(alignIndex1,alignIndex2));
		}
	}
}

