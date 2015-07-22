/*
 *  Common.h
 *  linkexons
 *
 *  Created by Andrew McPherson on 15/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

using namespace std;
using namespace boost;

enum Strand
{
	PlusStrand = 0,
	MinusStrand = 1,
};

int OtherStrand(int strand);
int InterpretStrand(const string& strand);

struct RefStrand
{
	union
	{
		struct
		{
			unsigned referenceIndex : 31;
			unsigned strand : 1;
		};
		
		unsigned int id;
	};
};

int OtherReadEnd(int readEnd);

struct Region
{
	int start;
	int end;	
};

inline bool operator==(const Region& r1, const Region& r2)
{
	return r1.start == r2.start && r1.end == r2.end;
}

inline int Length(const Region& r)
{
	return r.end - r.start + 1;
}

inline bool Overlap(const Region& r1, const Region& r2)
{
	if (r1.end < r2.start || r1.start > r2.end)
	{
		return false;
	}
	else
	{
		return true;
	}
}

typedef vector<Region> RegionVec;
typedef vector<Region>::iterator RegionVecIter;
typedef vector<Region>::const_iterator RegionVecConstIter;

typedef vector<RegionVec> RegionTable;
typedef vector<RegionVec>::iterator RegionTableIter;
typedef vector<RegionVec>::const_iterator RegionTableConstIter;

struct Location
{
	string refName;
	int strand;
	int start;
	int end;
};

typedef vector<Location> LocationVec;
typedef vector<Location>::iterator LocationVecIter;
typedef vector<Location>::const_iterator LocationVecConstIter;

typedef vector<LocationVec> LocationTable;
typedef vector<LocationVec>::iterator LocationTableIter;
typedef vector<LocationVec>::const_iterator LocationTableConstIter;

typedef map<int,LocationVec> LocationVecMap;
typedef map<int,LocationVec>::iterator LocationVecMapIter;
typedef map<int,LocationVec>::const_iterator LocationVecMapConstIter;

struct Counter
{
	Counter() : value(0) {}
	
	void operator++() { value++; }
	void operator++(int) { value++; }
	void operator+=(int a) { value += a; }
	operator int() { return value; }
	
	int value;
};

typedef vector<string> StringVec;
typedef vector<string>::iterator StringVecIter;
typedef vector<string>::const_iterator StringVecConstIter;

typedef vector<string> SequenceVec;
typedef vector<string>::iterator SequenceVecIter;
typedef vector<string>::const_iterator SequenceVecConstIter;

typedef vector<SequenceVec> SequenceTable;
typedef vector<SequenceVec>::iterator SequenceTableIter;
typedef vector<SequenceVec>::const_iterator SequenceTableConstIter;

typedef vector<int> IntegerVec;
typedef vector<int>::iterator IntegerVecIter;
typedef vector<int>::const_iterator IntegerVecConstIter;

typedef vector<IntegerVec> IntegerTable;
typedef vector<IntegerVec>::iterator IntegerTableIter;
typedef vector<IntegerVec>::const_iterator IntegerTableConstIter;

typedef vector<float> FloatVec;
typedef vector<float>::iterator FloatVecIter;
typedef vector<float>::const_iterator FloatVecConstIter;

typedef vector<FloatVec> FloatTable;
typedef vector<FloatVec>::iterator FloatTableIter;
typedef vector<FloatVec>::const_iterator FloatTableConstIter;

typedef vector<double> DoubleVec;
typedef vector<double>::iterator DoubleVecIter;
typedef vector<double>::const_iterator DoubleVecConstIter;

typedef vector<DoubleVec> DoubleTable;
typedef vector<DoubleVec>::iterator DoubleTableIter;
typedef vector<DoubleVec>::const_iterator DoubleTableConstIter;

typedef pair<int,int> IntegerPair;

typedef vector<IntegerPair> IntegerPairVec;
typedef vector<IntegerPair>::iterator IntegerPairVecIter;
typedef vector<IntegerPair>::const_iterator IntegerPairVecConstIter;

typedef vector<IntegerPairVec> IntegerPairTable;
typedef vector<IntegerPairVec>::iterator IntegerPairTableIter;
typedef vector<IntegerPairVec>::const_iterator IntegerPairTableConstIter;

typedef unordered_set<int> IntegerSet;
typedef unordered_set<int>::iterator IntegerSetIter;
typedef unordered_set<int>::const_iterator IntegerSetConstIter;

typedef unordered_map<int,int> IntegerMap;
typedef unordered_map<int,int>::iterator IntegerMapIter;
typedef unordered_map<int,int>::const_iterator IntegerMapConstIter;

typedef unordered_map<int,IntegerVec> IntegerVecMap;
typedef unordered_map<int,IntegerVec>::iterator IntegerVecMapIter;
typedef unordered_map<int,IntegerVec>::const_iterator IntegerVecMapConstIter;

typedef unordered_map<IntegerPair,int> IntegerPairMap;
typedef unordered_map<IntegerPair,int>::iterator IntegerPairMapIter;
typedef unordered_map<IntegerPair,int>::const_iterator IntegerPairMapConstIter;

void ReverseComplement(string& sequence);

inline bool operator==(const IntegerVec& vec1, const IntegerVec& vec2)
{
	return vec1.size() == vec2.size() && equal(vec1.begin(), vec1.end(), vec2.begin());
}

inline size_t hash_value(const IntegerVec& vec)
{
	size_t seed = 0;
	for (IntegerVecConstIter elementIter = vec.begin(); elementIter != vec.end(); elementIter++)
	{
		hash_combine(seed, *elementIter);
	}
    return seed;
}

struct ReadID
{
	union
	{
		struct
		{
			unsigned fragmentIndex : 31;
			unsigned readEnd : 1;
		};

		int id;
	};
};

struct ClusterID
{
	union
	{
		struct
		{
			unsigned clusterID : 31;
			unsigned clusterEnd : 1;
		};

		int id;
	};
};

struct CompactAlignment
{
	ReadID readID;
	RefStrand refStrand;
	Region region;
};

typedef vector<CompactAlignment> CompAlignVec;
typedef vector<CompactAlignment>::iterator CompAlignVecIter;
typedef vector<CompactAlignment>::const_iterator CompAlignVecConstIter;

typedef vector<CompAlignVec> CompAlignTable;
typedef vector<CompAlignVec>::iterator CompAlignTableIter;
typedef vector<CompAlignVec>::const_iterator CompAlignTableConstIter;

bool FragmentLessThan(const CompactAlignment& a1, const CompactAlignment& a2);

double normalpdf(double x, double mu, double sigma);

int FindMaxElement(const IntegerTable& clusters);
int FindMaxElement(const IntegerVecMap& clusters);

void CheckFile(const ios& file, const string& filename);

bool ParseTranscriptID(const string& transcriptID, string& gene, string& transcript);

void Transpose(const IntegerVecMap& clusters, IntegerVecMap& elements);

void Print(const IntegerVecMap& clusters);
void Print(const IntegerVec& cluster);

#endif
