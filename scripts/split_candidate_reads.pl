#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 split_catalog candidate_reads > candidate_catalog\n";

my $split_catalog_filename = shift;
my $candidate_reads_filename = shift;

die $usage if not defined $candidate_reads_filename;

# Read in split catalog
my @split_catalog = read_split_catalog($split_catalog_filename);

# Iterate over splits
foreach my $split_info (@split_catalog)
{
	my $split_prefix = $split_info->[0];
	my $start_fragment_index = $split_info->[1];
	my $end_fragment_index = $split_info->[2];
	
	my $split_candidates_filename = $split_prefix.".candidate.reads";
	
	open SPL, ">".$split_candidates_filename or die "Error: Unable to open $split_candidates_filename: $!\n";
	open CAN, $candidate_reads_filename or die "Error: Unable to open $candidate_reads_filename: $!\n";
	while (<CAN>)
	{
		my $line = $_;

		chomp;
		my @fields = split /\t/;

		my $fragment_index = $fields[1];

		if ($fragment_index >= $start_fragment_index and $fragment_index <= $end_fragment_index)
		{
			print SPL $line;
		}
	}
	close CAN;
	close SPL;
	
	print $split_candidates_filename."\n";
}

sub read_split_catalog
{
	my $filename = shift;
	my @catalog;

	open RSC, $filename or die "Error: Unable to open $filename: $!\n";
	while (<RSC>)
	{
		chomp;
		my @fields = split /\t/;
		
		push @catalog, [$fields[0],$fields[1],$fields[2]];
	}
	close RSC;
	
	return @catalog;
}

