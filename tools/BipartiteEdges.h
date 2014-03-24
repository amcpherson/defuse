/*
 *  BipartiteEdges.h
 *  tools
 *
 *  Created by Andrew McPherson on 10-08-04.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef BIPARTITEEDGES_H_
#define BIPARTITEEDGES_H_

#include "Common.h"

#include <vector>

using namespace std;

class BipartiteEdges
{
public:
	bool Read(istream& in);
	
	bool CalculateComponents(IntegerTable& leftVertices, IntegerTable& rightVertices);
	
private:
	IntegerTable mLeftToRight;
	IntegerTable mRightToLeft;	
};

#endif

