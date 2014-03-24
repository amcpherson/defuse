/*
 *  QuaterpartiteEdges.h
 *  tools
 *
 *  Created by Andrew McPherson on 10-09-02.
 *
 */

#ifndef QUATERPARTITEEDGES_H_
#define QUATERPARTITEEDGES_H_

#include "Common.h"

#include <vector>

using namespace std;

class QuaterpartiteEdges
{
public:
	QuaterpartiteEdges(const IntegerTable& dnaClusters, const IntegerPairVec& overlaps, const IntegerTable& rnaClusters);
	bool CalculateComponents(IntegerTable& dnaClusterComponents, IntegerTable& rnaClusterComponents);
	
private:
	int CreateMapping(const IntegerTable& dnaClusters, const IntegerPairVec& overlaps, const IntegerTable& rnaClusters);
	int EncodeDNACluster(int dnaClusterIndex);
	int EncodeDNAFragment(int dnaFragmentIndex);
	int EncodeRNACluster(int rnaClusterIndex);
	int EncodeRNAFragment(int rnaFragmentIndex);
	int DecodeDNACluster(int dnaClusterVertex);
	int DecodeDNAFragment(int dnaFragmentVertex);
	int DecodeRNACluster(int rnaClusterVertex);
	int DecodeRNAFragment(int rnaFragmentVertex);
	bool IsDNACluster(int vertex);
	bool IsRNACluster(int vertex);
	
	int mDNAClusterOffset;
	int mDNAFragmentOffset;
	int mRNAClusterOffset;
	int mRNAFragmentOffset;
	IntegerTable mEdges;
};

#endif


