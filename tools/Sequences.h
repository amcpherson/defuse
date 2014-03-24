/*
 *  Sequences.h
 *
 *  Created by Andrew McPherson on 07/16/11.
 *
 */

#ifndef SEQUENCES_H_
#define SEQUENCES_H_

#include "Common.h"

#include <map>
#include <string>
#include <iostream>

using namespace std;
using namespace boost;

class Sequences
{
public:
	void Read(const string& fastaFilename);
	void Get(const string& id, int start, int end, string& sequence) const;
	
private:
	unordered_map<string,string> mSequences;
};

#endif
