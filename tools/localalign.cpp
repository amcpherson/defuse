/*
 *  localalign.cpp
 *
 *  Created by Andrew McPherson on 05-09-11.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "SimpleAligner.h"
#include "DebugCheck.h"
#include "Common.h"

#include <fstream>
#include <iostream>
#include <string>
#include <tclap/CmdLine.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	double threshold;
	
	try
	{
		TCLAP::CmdLine cmd("Local realignment tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<double> thresholdArg("t","threshold","Percent Perfect Threshold",false,0.0,"float",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		threshold = thresholdArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	string line;
	int lineNumber = 0;
	while (getline(cin, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty line " << lineNumber << endl;
			exit(1);
		}
		
		vector<string> inputFields;
		split(inputFields, line, is_any_of("\t"));
		
		if (inputFields.size() < 3)
		{
			cerr << "Error: Format error for line " << lineNumber << endl;
			exit(1);
		}
		
		const string& id = inputFields[0];
		const string& reference = inputFields[1];
		const string& sequence = inputFields[2];
		
		int score = aligner.Align(reference, sequence);
		
		int maxScore = sequence.size() * matchScore;
		double percent = (double) score / (double) maxScore;
		
		if (percent < threshold)
		{
			continue;
		}

		cout << id << "\t" << score << "\t" << percent << endl;
	}
}
	
