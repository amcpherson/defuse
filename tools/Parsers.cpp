/*
 *  Parsers.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-09-02.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "Parsers.h"
#include "Common.h"
#include "DebugCheck.h"

#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_set.hpp>
#include "bam.h"
#include "sam.h"

using namespace boost;
using namespace std;

void ReadClusters(const string& clustersFilename, IntegerTable& clusters)
{
	// Open clusters file
	ifstream clustersFile(clustersFilename.c_str());
	if (!clustersFile)
	{
		cerr << "Error: unable to read from clusters file " << clustersFilename << endl;		
		exit(1);
	}
	
	// Parse file contents
	string line;
	int lineNumber = 0;
	while (getline(clustersFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty clusters line " << lineNumber << " of " << clustersFilename << endl;
			exit(1);
		}
		
		vector<string> clusterFields;
		split(clusterFields, line, is_any_of("\t"));
		
		if (clusterFields.size() < 3)
		{
			cerr << "Error: Format error for clusters line " << lineNumber << " of " << clustersFilename << endl;
			exit(1);
		}
		
		int clusterID = lexical_cast<int>(clusterFields[0]);
		int clusterEnd = lexical_cast<int>(clusterFields[1]);
		int fragmentIndex = lexical_cast<int>(clusterFields[2]);

		// Only read in cluster end 0
		if (clusterEnd != 0)
		{
			continue;
		}		

		if (clusterID < 0)
		{
			cerr << "Error: Invalid cluster ID for line " << lineNumber << " of " << clustersFilename << endl;
			exit(1);
		}
		
		clusters.resize(max(clusterID + 1, (int)clusters.size()));
		clusters[clusterID].push_back(fragmentIndex);
	}
	
	clustersFile.close();
}

void WriteClusters(const string& inClustersFilename, const string& outClustersFilename, const IntegerTable& clusters, int minClusterSize)
{
	typedef unordered_set<int> IntegerSet;
	typedef vector<IntegerSet> IntegerSetVec;
	
	// Create a hash of fragment indices for each cluster
	IntegerSetVec clusterFragments(clusters.size());
	for (int clusterID = 0; clusterID < clusters.size(); clusterID++)
	{
		if (clusters[clusterID].size() < minClusterSize)
		{
			continue;
		}
		
		for (IntegerVecConstIter fragmentIndexIter = clusters[clusterID].begin(); fragmentIndexIter != clusters[clusterID].end(); fragmentIndexIter++)
		{
			clusterFragments[clusterID].insert(*fragmentIndexIter);
		}
	}
	
	// Open output clusters file
	ofstream outClustersFile(outClustersFilename.c_str());
	if (!outClustersFile)
	{
		cerr << "Error: unable to write to clusters file " << outClustersFilename << endl;		
		exit(1);
	}
	
	// Open input clusters file
	ifstream inClustersFile(inClustersFilename.c_str());
	if (!inClustersFile)
	{
		cerr << "Error: unable to read from clusters file " << inClustersFilename << endl;		
		exit(1);
	}
	
	// Parse file contents
	string line;
	int lineNumber = 0;
	while (getline(inClustersFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty clusters line " << lineNumber << " of " << outClustersFilename << endl;
			exit(1);
		}
		
		vector<string> clusterFields;
		split(clusterFields, line, is_any_of("\t"));
		
		if (clusterFields.size() < 3)
		{
			cerr << "Error: Format error for clusters line " << lineNumber << " of " << outClustersFilename << endl;
			exit(1);
		}
		
		int clusterID = lexical_cast<int>(clusterFields[0]);
		int fragmentIndex = lexical_cast<int>(clusterFields[2]);
		
		if (clusterID < 0)
		{
			cerr << "Error: Invalid cluster ID for line " << lineNumber << " of " << outClustersFilename << endl;
			exit(1);
		}
		
		// Output cluster line if the fragment index is in this output cluster
		if (clusterFragments[clusterID].find(fragmentIndex) != clusterFragments[clusterID].end())
		{
			outClustersFile << line << endl;
		}
	}
	
	outClustersFile.close();
	inClustersFile.close();
}

void ReadOverlaps(const string& overlapsFilename, IntegerPairVec& overlaps)
{
	// Open overlaps file
	ifstream overlapsFile(overlapsFilename.c_str());
	if (!overlapsFile)
	{
		cerr << "Error: unable to read from overlaps file " << overlapsFilename << endl;		
		exit(1);
	}
	
	// Parse file contents
	string line;
	int lineNumber = 0;
	while (getline(overlapsFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty overlaps line " << lineNumber << " of " << overlapsFilename << endl;
			exit(1);
		}
		
		vector<string> overlapsFields;
		split(overlapsFields, line, is_any_of("\t"));
		
		if (overlapsFields.size() < 2)
		{
			cerr << "Error: Format error for overlaps line " << lineNumber << " of " << overlapsFilename << endl;
			exit(1);
		}
		
		int clusterID1 = lexical_cast<int>(overlapsFields[0]);
		int clusterID2 = lexical_cast<int>(overlapsFields[1]);
		
		overlaps.push_back(IntegerPair(clusterID1, clusterID2));
	}
	
	overlapsFile.close();
}

void ReadAlignments(const string& bamFilename, StringVec& referenceNames, CompAlignVec& alignments)
{
	cout << "Reading alignments and creating fragment name index" << endl;
	
	samfile_t* inBamFile = samopen(bamFilename.c_str(), "rb", 0);
	
	referenceNames.clear();
	for (int targetIndex = 0; targetIndex < inBamFile->header->n_targets; targetIndex++)
	{
		referenceNames.push_back(string(inBamFile->header->target_name[targetIndex]));
	}

	bam1_t bamEntry;
	memset(&bamEntry, 0, sizeof(bam1_t));
	while (bam_read1(inBamFile->x.bam, &bamEntry))
	{
		// Split qname into id and end
		// Fragment index encoded in fragment name
		string qname = string((char*)bamEntry.data);
		vector<string> qnameFields;
		split(qnameFields, qname, is_any_of("/"));
		
		CompactAlignment compactAlignment;
		compactAlignment.readID.fragmentIndex = lexical_cast<int>(qnameFields[0]);
		compactAlignment.readID.readEnd = (qnameFields[1] == "1") ? 0 : 1;
		compactAlignment.refStrand.referenceIndex = bamEntry.core.tid;
		compactAlignment.refStrand.strand = (bamEntry.core.flag & 0x10) ? MinusStrand : PlusStrand;
		compactAlignment.region.start = bamEntry.core.pos + 1;
		compactAlignment.region.end = bamEntry.core.pos + bamEntry.core.l_qseq;
		
		alignments.push_back(compactAlignment);
		
		memset(&bamEntry, 0, sizeof(bam1_t));
	}
	
	samclose(inBamFile);
}

void IntepretAlignString(const string& alignString, Location& alignRegion)
{
	string::size_type colonDividerPos = alignString.find_first_of(":");
	if (colonDividerPos == string::npos || colonDividerPos == 0)
	{
		cerr << "Error: Unable to interpret strand for " << alignString << endl;
		exit(1);
	}
	
	char strand = alignString[colonDividerPos - 1];
	
	if (strand == '+')
	{
		alignRegion.strand = PlusStrand;
	}
	else if (strand == '-')
	{
		alignRegion.strand = MinusStrand;
	}
	else
	{
		cerr << "Error: Unable to interpret strand for " << alignString << endl;
		exit(1);
	}
	
	vector<string> alignFields;
	split(alignFields, alignString, is_any_of("+-:"));
	
	if (alignFields.size() != 4)
	{
		cerr << "Error: Unable to interpret alignment string " << alignString << endl;
		exit(1);
	}
	
	alignRegion.refName = alignFields[0];
	alignRegion.start = lexical_cast<int>(alignFields[2]);
	alignRegion.end = lexical_cast<int>(alignFields[3]);
}

void ReadAlignRegionPairs(const string& filename, LocationVecMap& alignRegionPairs)
{
	ifstream alignRegionPairsFile(filename.c_str());
	if (!alignRegionPairsFile.good())
	{
		cerr << "Error: Unable to open align region pairs file " << filename << endl;
		exit(1);
	}
	
	string line;
	int lineNumber = 0;
	
	while (getline(alignRegionPairsFile, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> alignRegionFields;
		split(alignRegionFields, line, is_any_of("\t"));
		
		if (alignRegionFields.size() < 5)
		{
			continue;
		}
		
		int pairID = lexical_cast<int>(alignRegionFields[0]);
		int pairEnd = lexical_cast<int>(alignRegionFields[1]);
		
		DebugCheck(pairEnd == 0 || pairEnd == 1);
		
		Location alignRegion;
		alignRegion.refName = alignRegionFields[2];
		alignRegion.strand = InterpretStrand(alignRegionFields[3]);
		alignRegion.start = lexical_cast<int>(alignRegionFields[4]);
		alignRegion.end = lexical_cast<int>(alignRegionFields[5]);
		
		alignRegionPairs[pairID].resize(2);
		alignRegionPairs[pairID][pairEnd] = alignRegion;
	}
	
	alignRegionPairsFile.close();
}
