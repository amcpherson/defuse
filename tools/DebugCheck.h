/*
 *  DebugCheck.h
 *  findseqpath
 *
 *  Created by Andrew McPherson on 04/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DEBUGCHECK_H_
#define DEBUGCHECK_H_

void DebugCheckFailure(const char* expr, const char* file, int line);

#ifdef NO_DEBUG_CHECKS
#define DebugCheck(expr) (void)(0)
#else
#define DebugCheck(expr) ((expr) ? (void)(0) : DebugCheckFailure(#expr, __FILE__, __LINE__))
#endif

#endif

