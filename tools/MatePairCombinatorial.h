/*
 *  MatePairCombinatorial.h
 *
 *  Created by Andrew McPherson.
 *
 */

#ifndef MATEPAIRCOMBINATORIAL_H_
#define MATEPAIRCOMBINATORIAL_H_

#include "IMatePairClusterer.h"

#include <vector>
#include <string.h>

using namespace std;

class MatePairCombinatorial : public IMatePairClusterer
{
public:
	void Initialize(const Parameters& parameters);
	void CreateMatePair(int id, const Region& alignment1, int strand1, const Region& alignment2, int strand2, MatePair& matePair) const;
	void CalcBreakRegion(const Region& alignment, int strand, Region& breakRegion) const;
	void DoClustering(const MatePairVec& matePairs, IntegerTable& clusters) const;
	
private:
	struct FeatureRegion
	{
		int matePairIndex;
		int start;
		int end;
	};
	
	typedef vector<FeatureRegion> FeatureRegionVec;
	typedef vector<FeatureRegion>::iterator FeatureRegionVecIter;
	typedef vector<FeatureRegion>::const_iterator FeatureRegionVecConstIter;
	
	struct FeatureRegionStartLessThan
	{
		bool operator()(const FeatureRegion& region1, const FeatureRegion& region2)
		{
			return region1.start < region2.start;
		}
	};

	struct FeatureRegionEndLessThan
	{
		bool operator()(const FeatureRegion& region1, const FeatureRegion& region2)
		{
			return region1.end < region2.end;
		}
	};
	
	struct ClusterTuple
	{
		int cluster1;
		int cluster2;
		int cluster3;
		
		bool operator==(const ClusterTuple& other) const
		{
			return (memcmp(this, &other, sizeof(ClusterTuple)) == 0);
		}
	};
	
	struct ClusterTupleHash : std::unary_function<ClusterTuple,std::size_t>
	{
		inline std::size_t operator()(ClusterTuple const& clusterTuple) const
		{
			std::size_t seed = 0;
			boost::hash_combine(seed, clusterTuple.cluster1);
			boost::hash_combine(seed, clusterTuple.cluster2);
			boost::hash_combine(seed, clusterTuple.cluster3);
			return seed;
		}
	};	
	
	typedef unordered_map<ClusterTuple,IntegerVec,ClusterTupleHash> ClusterElementMap;
	typedef unordered_map<ClusterTuple,IntegerVec,ClusterTupleHash>::const_iterator ClusterElementMapConstIter;
	
	void StrandRemap(const Region& region, int strand, Region& remapped) const;
	void CalcBreakRegion(const Region& alignment, Region& breakRegion) const;
	void CreateBreakRegion(int matePairIndex, const Region& alignment, FeatureRegion& breakRegion) const;
	void CreateLengthRange(int matePairIndex, const Region& alignment1, const Region& alignment2, FeatureRegion& lengthRange) const;
	void CreateOverlapClusters(FeatureRegionVec& regions, IntegerTable& overlapClusters) const;
	void RemoveRedundant(IntegerTable& clusters) const;
	void SelectFeatureRegions(const IntegerVec& indices, const FeatureRegionVec& features, FeatureRegionVec& selected) const;
	
	int mMaxFragmentLength;
	int mMaxFragmentLengthDiff;
	int mMinClusterSize;
	int mResolution;
};

#endif
