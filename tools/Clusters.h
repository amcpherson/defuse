/*
 *  Clusters.h
 *
 *  Created by Andrew McPherson on 09-10-08.
 *
 */

#ifndef CLUSTERS_H_
#define CLUSTERS_H_

#include "Common.h"

#include <vector>
#include <boost/unordered_map.hpp>

using namespace std;
using namespace boost;

typedef IntegerTable ElementClusterMap;
typedef IntegerTableIter ElementClusterMapIter;
typedef IntegerTableConstIter ElementClusterMapConstIter;

class Clusters
{
public:
	Clusters();

	void OpenCluster();
	void AddToCluster(int element);
	int OpenClusterSize() const;
	int CloseCluster();
	
	int AddCluster(const IntegerVec& cluster);
	
	void Merge(Clusters& other);
	
	void RemoveRedundant();

	int GetNumClusters() const;
	const IntegerTable& GetClusters() const;
	const IntegerVec& GetCluster(int clusterIndex) const;

	void CreateElementClusterMap(ElementClusterMap& elementClusters) const;
	
	void Clear();

private:
	IntegerVec mCurrentCluster;
	IntegerTable mClusters;
};

#endif

