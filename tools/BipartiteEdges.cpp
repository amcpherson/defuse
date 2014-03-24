/*
 *  BipartiteEdges.cpp
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-04.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "BipartiteEdges.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>

using namespace std;


bool BipartiteEdges::Read(istream& in)
{
	string line;
	long int lineNumber = 0;

	cerr << "Reading edge file, '.'=1,000,000 edges" << endl;
	
	while (getline(in, line))
	{
		lineNumber++;

		if (lineNumber % 1000000 == 0)
		{
			cerr << ".";
			cerr.flush();
		}
		
		if (line.length() == 0)
		{
			continue;
		}
		
		vector<string> edgeFields;
		split(edgeFields, line, is_any_of("\t"));
		
		if (edgeFields.size() < 2)
		{
			continue;
		}
		
		int leftVertex = lexical_cast<int>(edgeFields[0]);
		int rightVertex = lexical_cast<int>(edgeFields[1]);

		if (leftVertex + 1 > mLeftToRight.size())
		{
			mLeftToRight.resize(leftVertex + 1);
		}

		if (rightVertex + 1 > mRightToLeft.size())
		{
			mRightToLeft.resize(rightVertex + 1);
		}

		mLeftToRight[leftVertex].push_back(rightVertex);
		mRightToLeft[rightVertex].push_back(leftVertex);
	}

	cerr << endl;
	cerr << "Read " << mLeftToRight.size() << " left vertices" << endl;
	cerr << "Read " << mRightToLeft.size() << " right vertices" << endl;
	cerr << "Read " << lineNumber << " edges" << endl;
	
	return true;
}

bool BipartiteEdges::CalculateComponents(IntegerTable& leftVertices, IntegerTable& rightVertices)
{
	if (mLeftToRight.empty() || mRightToLeft.empty())
	{
		return false;
	}
	
	cerr << "Finding components, '.'=10,000 components" << endl;

	IntegerVec leftVisited(mLeftToRight.size(), 0);
	IntegerVec rightVisited(mRightToLeft.size(), 0);
	
	int leftVisitedVertex = 0;
	int rightVisitedVertex = 0;
	
	enum VisitEnum
	{
		LeftVisit = 0,
		RightVisit = 1,
	};
	
	IntegerPairVec vertexStack;

	// Create the first pair of vertex sets for the first component
	leftVertices.push_back(IntegerVec());
	rightVertices.push_back(IntegerVec());

	// Select the first left vertex as the start
	vertexStack.push_back(IntegerPair(0,LeftVisit));
	
	while (!vertexStack.empty())
	{
		IntegerPair nextVertex = vertexStack.back();
		vertexStack.pop_back();
		
		if (nextVertex.second == LeftVisit)
		{
			leftVertices.back().push_back(nextVertex.first);
			leftVisited[nextVertex.first] = 1;

			for (IntegerVecConstIter rightVertexIter = mLeftToRight[nextVertex.first].begin(); rightVertexIter != mLeftToRight[nextVertex.first].end(); rightVertexIter++)
			{
				if (!rightVisited[*rightVertexIter])
				{
					vertexStack.push_back(IntegerPair(*rightVertexIter,RightVisit));
				}
			}
		}
		else
		{
			rightVertices.back().push_back(nextVertex.first);
			rightVisited[nextVertex.first] = 1;
			
			for (IntegerVecConstIter leftVertexIter = mRightToLeft[nextVertex.first].begin(); leftVertexIter != mRightToLeft[nextVertex.first].end(); leftVertexIter++)
			{
				if (!leftVisited[*leftVertexIter])
				{
					vertexStack.push_back(IntegerPair(*leftVertexIter,LeftVisit));
				}
			}
		}
		
		if (vertexStack.empty())
		{
			while (leftVisitedVertex < leftVisited.size() && leftVisited[leftVisitedVertex] != 0)
			{
				leftVisitedVertex++;
			}
			
			while (rightVisitedVertex < rightVisited.size() && rightVisited[rightVisitedVertex] != 0)
			{
				rightVisitedVertex++;
			}
			
			if (leftVisitedVertex < leftVisited.size())
			{
				vertexStack.push_back(IntegerPair(leftVisitedVertex,LeftVisit));				
			}
			else if (rightVisitedVertex < rightVisited.size())
			{
				vertexStack.push_back(IntegerPair(rightVisitedVertex,RightVisit));				
			}
			
			if (!vertexStack.empty())
			{
				leftVertices.push_back(IntegerVec());
				rightVertices.push_back(IntegerVec());

				if (leftVertices.size() % 10000 == 0)
				{
					cerr << ".";
					cerr.flush();
				}
			}
		}
	}

	cerr << endl;
	cerr << "Found " << leftVertices.size() << " components" << endl;
	
	return true;
}


