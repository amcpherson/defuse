/*
 *  HierarchicalClusterer.h
 *
 *  Created by Andrew McPherson on 1/27/10.
 *
 */

#ifndef HIERARCHICALCLUSTERER_H_
#define HIERARCHICALCLUSTERER_H_

#include "Common.h"

class HierarchicalClusterer
{
public:
	void DoClustering(const DoubleTable& distances, double threshold);
	const IntegerTable& GetClusters();

private:
	IntegerTable mClusters;
};

#endif

