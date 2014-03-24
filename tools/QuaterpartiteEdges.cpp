/*
 *  QuaterpartiteEdges.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-09-02.
 *
 */

#include "QuaterpartiteEdges.h"
#include "Common.h"
#include "DebugCheck.h"

#include <iostream>

using namespace std;


QuaterpartiteEdges::QuaterpartiteEdges(const IntegerTable& dnaClusters, const IntegerPairVec& overlaps, const IntegerTable& rnaClusters)
{
	int maxVertex = CreateMapping(dnaClusters, overlaps, rnaClusters);
	
	mEdges.resize(maxVertex);
	
	for (int dnaClusterIndex = 0; dnaClusterIndex < dnaClusters.size(); dnaClusterIndex++)
	{
		int dnaClusterVertex = EncodeDNACluster(dnaClusterIndex);

		for (IntegerVecConstIter dnaFragmentIndexIter = dnaClusters[dnaClusterIndex].begin(); dnaFragmentIndexIter != dnaClusters[dnaClusterIndex].end(); dnaFragmentIndexIter++)
		{
			int dnaReadVertex = EncodeDNAFragment(*dnaFragmentIndexIter);
			
			mEdges[dnaClusterVertex].push_back(dnaReadVertex);
			mEdges[dnaReadVertex].push_back(dnaClusterVertex);
		}
	}
	
	for (int rnaClusterIndex = 0; rnaClusterIndex < rnaClusters.size(); rnaClusterIndex++)
	{
		int rnaClusterVertex = EncodeRNACluster(rnaClusterIndex);

		for (IntegerVecConstIter rnaFragmentIndexIter = rnaClusters[rnaClusterIndex].begin(); rnaFragmentIndexIter != rnaClusters[rnaClusterIndex].end(); rnaFragmentIndexIter++)
		{
			int rnaReadVertex = EncodeRNAFragment(*rnaFragmentIndexIter);
			
			mEdges[rnaClusterVertex].push_back(rnaReadVertex);
			mEdges[rnaReadVertex].push_back(rnaClusterVertex);
		}
	}
	
	for (IntegerPairVecConstIter overlapIter = overlaps.begin(); overlapIter != overlaps.end(); overlapIter++)
	{
		int dnaClusterVertex = EncodeDNACluster(overlapIter->first);
		int rnaClusterVertex = EncodeRNACluster(overlapIter->second);
		
		mEdges[dnaClusterVertex].push_back(rnaClusterVertex);
		mEdges[rnaClusterVertex].push_back(dnaClusterVertex);
	}
}

int QuaterpartiteEdges::CreateMapping(const IntegerTable& dnaClusters, const IntegerPairVec& overlaps, const IntegerTable& rnaClusters)
{
	int maxDNAClusterIndex = dnaClusters.size();
	int maxRNAClusterIndex = rnaClusters.size();
	int maxDNAFragmentIndex = FindMaxElement(dnaClusters);
	int maxRNAFragmentIndex = FindMaxElement(rnaClusters);
	
	mDNAClusterOffset = 0;
	mDNAFragmentOffset = mDNAClusterOffset + maxDNAClusterIndex;
	mRNAClusterOffset = mDNAFragmentOffset + maxDNAFragmentIndex;
	mRNAFragmentOffset = mRNAClusterOffset + maxRNAClusterIndex;
	
	int maxVertex = mRNAFragmentOffset + maxRNAFragmentIndex;
	
	return maxVertex;
}

int QuaterpartiteEdges::EncodeDNACluster(int dnaClusterIndex)
{
	return dnaClusterIndex + mDNAClusterOffset;
}

int QuaterpartiteEdges::EncodeDNAFragment(int dnaFragmentIndex)
{
	return dnaFragmentIndex + mDNAFragmentOffset;
}

int QuaterpartiteEdges::EncodeRNACluster(int rnaClusterIndex)
{
	return rnaClusterIndex + mRNAClusterOffset;
}

int QuaterpartiteEdges::EncodeRNAFragment(int rnaFragmentIndex)
{
	return rnaFragmentIndex + mRNAFragmentOffset;
}

int QuaterpartiteEdges::DecodeDNACluster(int dnaClusterVertex)
{
	return dnaClusterVertex - mDNAClusterOffset;
}

int QuaterpartiteEdges::DecodeDNAFragment(int dnaFragmentVertex)
{
	return dnaFragmentVertex - mDNAFragmentOffset;
}

int QuaterpartiteEdges::DecodeRNACluster(int rnaClusterVertex)
{
	return rnaClusterVertex - mRNAClusterOffset;
}

int QuaterpartiteEdges::DecodeRNAFragment(int rnaFragmentVertex)
{
	return rnaFragmentVertex - mRNAFragmentOffset;
}

bool QuaterpartiteEdges::IsDNACluster(int vertex)
{
	return (vertex >= mDNAClusterOffset && vertex < mDNAFragmentOffset);
}

bool QuaterpartiteEdges::IsRNACluster(int vertex)
{
	return (vertex >= mRNAClusterOffset && vertex < mRNAFragmentOffset);
}

bool QuaterpartiteEdges::CalculateComponents(IntegerTable& dnaClusterComponents, IntegerTable& rnaClusterComponents)
{
	IntegerTable components;
	
	if (mEdges.empty())
	{
		return false;
	}
	
	cerr << "Finding components, '.'=10,000 components" << endl;
	
	// Visited flag for each vertex	
	IntegerVec visited(mEdges.size(), 0);
	
	// Next unvisited vertex in search
	int nextUnvisitedVertex = 0;
	
	// Stack of vertices to visit in search
	IntegerVec vertexStack;
	
	// Create the vertex sets for the first component
	components.push_back(IntegerVec());
	
	// Select the first vertex as the start
	vertexStack.push_back(0);
	
	while (!vertexStack.empty())
	{
		int nextVertex = vertexStack.back();
		vertexStack.pop_back();
		
		components.back().push_back(nextVertex);
		visited[nextVertex] = 1;
		
		for (IntegerVecConstIter vertexIter = mEdges[nextVertex].begin(); vertexIter != mEdges[nextVertex].end(); vertexIter++)
		{
			if (!visited[*vertexIter])
			{
				vertexStack.push_back(*vertexIter);
			}
		}
		
		if (vertexStack.empty())
		{
			while (nextUnvisitedVertex < visited.size() && visited[nextUnvisitedVertex] != 0)
			{
				nextUnvisitedVertex++;
			}
			
			if (nextUnvisitedVertex < visited.size())
			{
				vertexStack.push_back(nextUnvisitedVertex);				
			}
			
			if (!vertexStack.empty())
			{
				components.push_back(IntegerVec());
				
				if (components.size() % 10000 == 0)
				{
					cerr << ".";
					cerr.flush();
				}
			}
		}
	}
	
	cerr << endl;
	cerr << "Found " << components.size() << " components" << endl;
	
	dnaClusterComponents.clear();
	rnaClusterComponents.clear();
	
	for (IntegerTableIter componentIter = components.begin(); componentIter != components.end(); componentIter++)
	{
		dnaClusterComponents.push_back(IntegerVec());
		rnaClusterComponents.push_back(IntegerVec());
		
		for (IntegerVecIter vertexIter = componentIter->begin(); vertexIter != componentIter->end(); vertexIter++)
		{
			if (IsDNACluster(*vertexIter))
			{
				dnaClusterComponents.back().push_back(DecodeDNACluster(*vertexIter));
			}
			
			if (IsRNACluster(*vertexIter))
			{
				rnaClusterComponents.back().push_back(DecodeRNACluster(*vertexIter));
			}
		}
	}
	
	return true;
}


