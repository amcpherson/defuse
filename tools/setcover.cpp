/*
 *  setcover.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-24.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Parsers.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <tclap/CmdLine.h>
#include <boost/algorithm/string.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_set.hpp>

using namespace boost;
using namespace boost::bimaps;
using namespace std;


void SetCover(const IntegerTable& clusters, IntegerTable& solution)
{
	// Solution contains all clusters, indexed similarly, initially empty
	solution = IntegerTable(clusters.size());

	// Find max element size
	int maxElement = FindMaxElement(clusters);
	
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
				solution[nextClusterIndex].push_back(element);
				
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

int main(int argc, char* argv[])
{
	string clustersFilename;
	int minClusterSize;
	string outClustersFilename;
	
	try
	{
		TCLAP::CmdLine cmd("Set cover for maximum parsimony");
		TCLAP::ValueArg<string> clustersFilenameArg("c","clusters","Clusters Filename",true,"","string",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> outClustersFilenameArg("o","outclust","Output Clusters Filename",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		clustersFilename = clustersFilenameArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
		outClustersFilename = outClustersFilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	cout << "Reading clusters" << endl;
	
	IntegerTable clusters;	
	ReadClusters(clustersFilename, clusters);
	
	cout << "Calculating set cover solution" << endl;
	IntegerTable solutionClusters;
	SetCover(clusters, solutionClusters);
	
	cout << "Writing out clusters" << endl;
	WriteClusters(clustersFilename, outClustersFilename, solutionClusters, minClusterSize);
}

