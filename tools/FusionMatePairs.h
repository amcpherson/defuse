/*
 *  FusionMatePairs.h
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#ifndef FUSIONMATEPAIRS_H_
#define FUSIONMATEPAIRS_H_

#include "Alignments.h"
#include "IMatePairClusterer.h"

#include <vector>

using namespace std;

class IAlignmentStream;
 
class FusionMatePairs
{
public:
	FusionMatePairs(const Alignments& alignments, const IMatePairClusterer* clusterer, int minClusterSize);
	
	void CreateMatePairClusters(IntegerPairTable& clusters);
	
private:
	const Alignments& mAlignments;
	const IMatePairClusterer* mClusterer;
	int mMinClusterSize;
};

#endif
