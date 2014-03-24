/*
 *  clustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Common.h"
#include "DebugCheck.h"
#include "Indexer.h"
#include "ReadStream.h"
#include "Sequences.h"
#include "SimpleAligner.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <tclap/CmdLine.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
using namespace std;


struct CompactPosition
{
	RefStrand refStrand;
	int position;
};

typedef vector<CompactPosition> CompPosVec;
typedef vector<CompactPosition>::iterator CompPosVecIter;
typedef vector<CompactPosition>::const_iterator CompPosVecConstIter;

int main(int argc, char* argv[])
{
	int matchScore;
	int misMatchScore;
	int gapScore;
	double threshold;
	int searchLength;
	string referenceFasta;
	string readSeqs1Filename;
	string readSeqs2Filename;
	
	try
	{
		TCLAP::CmdLine cmd("Mate Realignment Tool");
		TCLAP::ValueArg<int> matchScoreArg("m","match","Match Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> misMatchScoreArg("x","mismatch","Mismatch Score",true,0,"int",cmd);
		TCLAP::ValueArg<int> gapScoreArg("g","gap","Gap Score",true,0,"int",cmd);
		TCLAP::ValueArg<double> thresholdArg("t","threshold","Percent Perfect Threshold",false,0.0,"float",cmd);
		TCLAP::ValueArg<int> searchLengthArg("s","searchlength","Search Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> referenceFastaArg("r","reference","Reference Sequences Fasta",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs1FilenameArg("1","seq1","End 1 Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs2FilenameArg("2","seq2","End 2 Sequences",true,"","string",cmd);
		cmd.parse(argc,argv);
		
		matchScore = matchScoreArg.getValue();
		misMatchScore = misMatchScoreArg.getValue();
		gapScore = gapScoreArg.getValue();
		threshold = thresholdArg.getValue();
		searchLength = searchLengthArg.getValue();
		referenceFasta = referenceFastaArg.getValue();
		readSeqs1Filename = readSeqs1FilenameArg.getValue();
		readSeqs2Filename = readSeqs2FilenameArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	unordered_map<int,CompPosVec> readAlignments;
	NameIndex referenceNames;
	
	string line;
	int lineNumber = 0;
	while (getline(cin, line))
	{
		lineNumber++;
		
		if (line.length() == 0)
		{
			cerr << "Error: Empty alignment line " << lineNumber << endl;
			exit(1);
		}
		
		if (line[0] == '@')
		{
			continue;
		}
		
		vector<string> alignmentFields;
		split(alignmentFields, line, is_any_of("\t"));
		
		if (alignmentFields.size() < 10)
		{
			cerr << "Error: Format error for alignment line " << lineNumber << endl;
			exit(1);
		}
		
		string qname = alignmentFields[0];
		int flag = lexical_cast<int>(alignmentFields[1]);
		string rname = alignmentFields[2];
		int pos = lexical_cast<int>(alignmentFields[3]);
		string seq = alignmentFields[9];
		
		if (rname == "*")
		{
			continue;
		}
		
		// Interpret strand
		int strand;
		if ((flag & 0x0010) == 0)
		{
			strand = PlusStrand;
		}
		else
		{
			strand = MinusStrand;
		}
		
		// Split qname into id and end
		vector<string> qnameFields;
		split(qnameFields, qname, is_any_of("/"));
		if (qnameFields.size() != 2 || (qnameFields[1] != "1" && qnameFields[1] != "2"))
		{
			cerr << "Error: Unable to interpret qname for alignment line " << lineNumber << endl;
			exit(1);
		}
		
		// Retreive fragment name and read end
		string fragmentName = qnameFields[0];
		int readEnd = (qnameFields[1] == "1") ? 0 : 1;
		
		// 1-based start and end position
		int start = lexical_cast<int>(pos);
		int end = start + seq.length() - 1;
		
		// Fragment index encoded in fragment name
		int fragmentIndex = lexical_cast<int>(qnameFields[0]);
		
		ReadID readID;
		readID.readEnd = readEnd;
		readID.fragmentIndex = fragmentIndex;

		CompactPosition alignment;
		alignment.refStrand.referenceIndex = referenceNames.Index(rname);
		alignment.refStrand.strand = strand;
		alignment.position = (strand == PlusStrand) ? start : end;
		
		readAlignments[readID.id].push_back(alignment);
	}
	
	cerr << "Read alignments" << endl;
	
	Sequences referenceSequences;
	referenceSequences.Read(referenceFasta);
	
	cerr << "Read reference fasta" << endl;
	
	IReadStream* readSeqsStream[2];
	readSeqsStream[0] = IReadStream::Create(readSeqs1Filename);
	readSeqsStream[1] = IReadStream::Create(readSeqs2Filename);
	
	if (!readSeqsStream[0] || !readSeqsStream[1])
	{
		cout << "Error: unable to read sequences" << endl;
		exit(1);
	}
	
	SimpleAligner aligner(matchScore, misMatchScore, gapScore);
	
	for (int readSeqsIndex = 0; readSeqsIndex <= 1; readSeqsIndex++)
	{
		RawRead rawRead;
		while (readSeqsStream[readSeqsIndex]->GetNextRead(rawRead))
		{
			ReadID readID;
			readID.fragmentIndex = lexical_cast<int>(rawRead.fragment);
			readID.readEnd = rawRead.readEnd;
			
			ReadID otherReadID;
			otherReadID.fragmentIndex = lexical_cast<int>(rawRead.fragment);
			otherReadID.readEnd = OtherReadEnd(rawRead.readEnd);
			
			if (readAlignments.find(otherReadID.id) != readAlignments.end())
			{
				for(CompPosVecConstIter posIter = readAlignments[otherReadID.id].begin(); posIter != readAlignments[otherReadID.id].end(); posIter++)
				{
					string referenceSequence;
					
					if (posIter->refStrand.strand == PlusStrand)
					{
						referenceSequences.Get(referenceNames.Get(posIter->refStrand.referenceIndex), posIter->position, posIter->position + searchLength, referenceSequence);
						
						ReverseComplement(referenceSequence);
					}
					else
					{
						referenceSequences.Get(referenceNames.Get(posIter->refStrand.referenceIndex), posIter->position - searchLength, posIter->position, referenceSequence);						
					}
					
					int score = aligner.Align(referenceSequence, rawRead.sequence);
					
					int maxScore = rawRead.sequence.size() * matchScore;
					double percent = (double) score / (double) maxScore;
					
					if (percent < threshold)
					{
						continue;
					}
					
					cout << readID.fragmentIndex << "\t" << score << "\t" << percent << endl;
				}
			}
		}
	}
}

