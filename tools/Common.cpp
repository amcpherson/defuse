/*
 *  Common.cpp
 *  linkexons
 *
 *  Created by Andrew McPherson on 15/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Common.h"

#include <iostream>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;
using namespace boost;


int OtherStrand(int strand)
{
	return (1 - strand);
}

int OtherReadEnd(int readEnd)
{
	return (1 - readEnd);
}

void ReverseComplement(string& sequence)
{
	reverse(sequence.begin(), sequence.end());
	
	for (int seqIndex = 0; seqIndex < sequence.size(); seqIndex++)
	{
		char nucleotide = sequence[seqIndex];
		
		switch (nucleotide)
		{
			case 'A': nucleotide = 'T'; break;
			case 'C': nucleotide = 'G'; break;
			case 'T': nucleotide = 'A'; break;
			case 'G': nucleotide = 'C'; break;
			case 'a': nucleotide = 't'; break;
			case 'c': nucleotide = 'g'; break;
			case 't': nucleotide = 'a'; break;
			case 'g': nucleotide = 'c'; break;
		}
		
		sequence[seqIndex] = nucleotide;
	}
}

bool FragmentLessThan(const CompactAlignment& a1, const CompactAlignment& a2)
{
	return a1.readID.fragmentIndex < a2.readID.fragmentIndex;
}

double normalpdf(double x, double mu, double sigma)
{
	double coeff = 1.0 / ((double)sigma * sqrt(2 * M_PI));
	
	double dist = (((double)x - (double)mu) / (double)sigma);
	double exponent = -0.5 * dist * dist;
	
	return coeff * exp(exponent);
}

int FindMaxElement(const IntegerTable& clusters)
{
	int maxElement = -1;
	for (IntegerTableConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		for (IntegerVecConstIter elementIter = clusterIter->begin(); elementIter != clusterIter->end(); elementIter++)
		{
			if (*elementIter < 0)
			{
				cerr << "Error: negative elements not permitted" << endl;
				exit(1);
			}
			
			maxElement = max(maxElement, *elementIter);
		}
	}
	
	return maxElement;
}

int InterpretStrand(const string& strand)
{
	if (strand == "+")
	{
		return PlusStrand;
	}
	else if (strand == "-")
	{
		return MinusStrand;
	}
	else
	{
		cerr << "Error: Unable to intepret strand " << strand << endl;
		exit(1);			
	}
}

void CheckFile(const ios& file, const string& filename)
{
	if (!file.good())
	{
		cerr << "Error: Unable to open " << filename << endl;
		exit(1);
	}	
}

bool ParseTranscriptID(const string& transcriptID, string& gene, string& transcript)
{
	StringVec idFields;
	split(idFields, transcriptID, is_any_of("|"));

	if (idFields.size() < 2)
	{
		return false;
	}

	gene = idFields[0];
	transcript = idFields[1];

	return true;
}

int FindMaxElement(const IntegerVecMap& clusters)
{
	int maxElement = -1;
	for (IntegerVecMapConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		for (IntegerVecConstIter elementIter = clusterIter->second.begin(); elementIter != clusterIter->second.end(); elementIter++)
		{
			if (*elementIter < 0)
			{
				cerr << "Error: negative elements not permitted" << endl;
				exit(1);
			}
			
			maxElement = max(maxElement, *elementIter);
		}
	}
	
	return maxElement;
}

void Transpose(const IntegerVecMap& clusters, IntegerVecMap& elements)
{
	for (IntegerVecMapConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		for (IntegerVecConstIter elementIter = clusterIter->second.begin(); elementIter != clusterIter->second.end(); elementIter++)
		{
			elements[*elementIter].push_back(clusterIter->first);
		}
	}
}

void Print(const IntegerVecMap& clusters)
{
	for (IntegerVecMapConstIter clusterIter = clusters.begin(); clusterIter != clusters.end(); clusterIter++)
	{
		for (IntegerVecConstIter elementIter = clusterIter->second.begin(); elementIter != clusterIter->second.end(); elementIter++)
		{
			cout << clusterIter->first << "\t" << *elementIter << endl;
		}
	}
}

void Print(const IntegerVec& cluster)
{
	for (IntegerVecConstIter elementIter = cluster.begin(); elementIter != cluster.end(); elementIter++)
	{
		cout << *elementIter << endl;
	}
}


