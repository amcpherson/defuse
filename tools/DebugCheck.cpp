/*
 *  DebugCheck.cpp
 *  findseqpath
 *
 *  Created by Andrew McPherson on 04/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "DebugCheck.h"

#include <iostream>

void DebugCheckFailure(const char* expr, const char* file, int line)
{
	std::cout << "Error: " << expr << " failed on line: " << line << " of " << file << std::endl;
}
