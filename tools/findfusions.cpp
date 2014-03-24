/*
 *  clustermatepairs.cpp
 *
 *  Created by Andrew McPherson on 28/09/09.
 *
 */

#include "Alignments.h"
#include "DebugCheck.h"
#include "IMatePairClusterer.h"
#include "MatePairCombinatorial.h"
#include "MatePairProbabilistic.h"
#include "ExonRegions.h"
#include "FusionMatePairs.h"
#include "FusionPredictions.h"

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <tclap/CmdLine.h>

using namespace std;

int main(int argc, char* argv[])
{
	string referenceFilename;
	string exonRegionsFilename;
	string discordantFilename;
	string partialFilename;
	string readSeqs1Filename;
	string readSeqs2Filename;
	string outputDirectory;
	double fragmentLengthMean;
	double fragmentLengthStdDev;
	int readLength;
	double precision;
	int minClusterSize;
	string algorithm;

	try
	{
		TCLAP::CmdLine cmd("Mate Pair Clustering Tool");
		TCLAP::ValueArg<string> referenceArg("r","reference","Reference Fasta Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> exonRegionsArg("e","exonregions","Exon Regions Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> discordantArg("d","discordant","Discordant Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> partialArg("p","partial","Partial Alignments Filename",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs1Arg("1","seq1","End 1 Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> readSeqs2Arg("2","seq2","End 2 Sequences",true,"","string",cmd);
		TCLAP::ValueArg<string> outputArg("o","output","Output Directory",true,"","string",cmd);
		TCLAP::ValueArg<double> fragmentLengthMeanArg("u","fragmentmean","Fragment Length Mean",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> fragmentLengthStdDevArg("s","fragmentstddev","Fragment Length Standard Deviation",true,-1,"integer",cmd);
		TCLAP::ValueArg<int> readLengthArg("l","readlength","Read Length",true,-1,"integer",cmd);
		TCLAP::ValueArg<double> precisionArg("a","precision","Precision",true,-1,"double",cmd);
		TCLAP::ValueArg<int> minClusterSizeArg("m","minclustersize","Minimum Cluster Size",true,-1,"integer",cmd);
		TCLAP::ValueArg<string> algorithmArg("z","algorithm","Algorithm, combinatorial or probabilistic",true,"","string",cmd);
		cmd.parse(argc,argv);

		referenceFilename = referenceArg.getValue();
		exonRegionsFilename = exonRegionsArg.getValue();
		discordantFilename = discordantArg.getValue();
		partialFilename = partialArg.getValue();
		readSeqs1Filename = readSeqs1Arg.getValue();
		readSeqs2Filename = readSeqs2Arg.getValue();
		outputDirectory = outputArg.getValue();
		fragmentLengthMean = fragmentLengthMeanArg.getValue();
		fragmentLengthStdDev = fragmentLengthStdDevArg.getValue();
		readLength = readLengthArg.getValue();
		precision = precisionArg.getValue();
		minClusterSize = minClusterSizeArg.getValue();
		algorithm = algorithmArg.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
		exit(1);
	}
	
	Alignments alignments;
	
	cout << "Reading alignment file" << endl;
	IAlignmentStream* discordantStream = IAlignmentStream::Create(discordantFilename);
	if (!discordantStream)
	{
		cout << "Error: unable to read discordant alignments" << endl;		
		exit(1);
	}
	alignments.Read(discordantStream);
	
	cout << "Removing redundant alignments" << endl;
	alignments.RemoveRedundant();

	cout << "Filtering canonical alignments" << endl;
	alignments.FilterEnsemblTranscripts();
	
	// Create clusterer according to requested algorithm
	IMatePairClusterer* clusterer = 0;
	if (algorithm == "combinatorial")
	{
		clusterer = new MatePairCombinatorial();
	}
	else if (algorithm == "probabilistic")
	{
		clusterer = new MatePairProbabilistic();
	}
	else
	{
		cout << "Error: invalid algorithm" << endl;		
		exit(1);
	}
	
	// Initialize clusterer
	cout << "Initializing clusterer" << endl;
	IMatePairClusterer::Parameters clustererParams;
	clustererParams.fragmentMean = fragmentLengthMean;
	clustererParams.fragmentStdDev = fragmentLengthStdDev;
	clustererParams.readLength = readLength;
	clustererParams.precision = precision;
	clustererParams.minClusterSize = minClusterSize;
	clusterer->Initialize(clustererParams);
	
	// Do the clustering using the clusterer
	cout << "Doing Clustering" << endl;
	IntegerPairTable matePairClusters;
	FusionMatePairs fusionMatePairs(alignments, clusterer, minClusterSize);
	fusionMatePairs.CreateMatePairClusters(matePairClusters);
	
	// Read exon regions
	ifstream exonRegionsFile(exonRegionsFilename.c_str());
	if (!exonRegionsFile)
	{
		cout << "Error: unable to read exon regions" << endl;		
		exit(1);
	}
	
	ExonRegions exonRegions;
	exonRegions.Read(exonRegionsFile);
		
	// Create the fusion predictions
	int minFragmentSize = fragmentLengthMean - 3 * fragmentLengthStdDev;
	int maxFragmentSize = fragmentLengthMean + 3 * fragmentLengthStdDev;
			
	Parameters parameters(minFragmentSize, maxFragmentSize, readLength, 28, minClusterSize, outputDirectory);	
	FusionPredictions fusionPredictions(alignments, parameters, matePairClusters, exonRegions);
	
	cout << "Filtering low coverage clusters" << endl;	
	fusionPredictions.FilterLowCoverage();
	
	cout << "Finding cluster reads" << endl;
	IAlignmentStream* partialStream = IAlignmentStream::Create(partialFilename);
	if (!partialStream)
	{
		cout << "Error: unable to read partial alignments" << endl;		
		exit(1);
	}
	
	discordantStream->Reset();
	
	discordantStream = IAlignmentStream::Create(discordantFilename);
	
	fusionPredictions.ClearClusterReads();
	fusionPredictions.FindSpanningClusterReads();
	fusionPredictions.FindHangingClusterReads(partialStream);
	fusionPredictions.FindHangingClusterReads(discordantStream);
	
	Reference reference(referenceFilename);
	
	cout << "Retrieving cluster sequences" << endl;
	fusionPredictions.RetrieveClusterSequences(reference);
	
	cout << "Retrieving partner sequences" << endl;
	fusionPredictions.RetrievePartnerSequences(reference);
	
	cout << "Retrieving read sequences" << endl;
	IReadStream* readSeqs1Stream = IReadStream::Create(readSeqs1Filename);
	IReadStream* readSeqs2Stream = IReadStream::Create(readSeqs2Filename);
	
	if (!readSeqs1Stream || !readSeqs2Stream)
	{
		cout << "Error: unable to read sequences" << endl;
		exit(1);
	}
	
	fusionPredictions.RetrieveReadSequences(readSeqs1Stream);
	fusionPredictions.RetrieveReadSequences(readSeqs2Stream);
	
	cout << "Aligning cluster reads" << endl;
	fusionPredictions.FindBreaks();
	
	cout << "Printing Results" << endl;		
	fusionPredictions.PrintClusters();
}

