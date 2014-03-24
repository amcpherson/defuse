/*
 *  Parsers.h
 *  tools
 *
 *  Created by Andrew McPherson on 10-09-02.
 *
 */

#include "Common.h"

void ReadClusters(const string& clustersFilename, IntegerTable& clusters);
void WriteClusters(const string& inClustersFilename, const string& outClustersFilename, const IntegerTable& clusters, int minClusterSize);
void IntepretAlignString(const string& alignString, Location& alignRegion);
void ReadAlignRegionPairs(const string& filename, LocationVecMap& alignRegionPairs);

