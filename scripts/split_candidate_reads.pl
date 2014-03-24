#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use FileCache;
no strict 'refs';

my $usage = "Usage: $0 split_catalog candidate_reads > candidate_catalog\n";

my $split_catalog_filename = shift;
my $candidate_reads_filename = shift;

die $usage if not defined $candidate_reads_filename;

# Read in split catalog
my @split_catalog = read_split_catalog($split_catalog_filename);

# Open files for writing
# Create binning for fragment indices
my $bin_size = 100000;
my %split_bins;
my @split_candidates_filenames;
foreach my $split_info (@split_catalog)
{
        my $split_prefix = $split_info->[0];
        my $start_fragment_index = $split_info->[1];
        my $end_fragment_index = $split_info->[2];

        my $split_candidates_filename = $split_prefix.".candidate.reads";
	
	cacheout $split_candidates_filename;
	
	my $start_bin = int($start_fragment_index / $bin_size);
	my $end_bin = int($end_fragment_index / $bin_size);

	foreach my $bin ($start_bin..$end_bin)
	{
		push @{$split_bins{$bin}}, [$start_fragment_index, $end_fragment_index, $split_candidates_filename];
	}

	push @split_candidates_filenames, $split_candidates_filename."\n";
}

# Output each line to the appropriate file
open CAN, $candidate_reads_filename or die "Error: Unable to open $candidate_reads_filename: $!\n";
while (<CAN>)
{
	my $line = $_;

	chomp;
	my @fields = split /\t/;

	my $fragment_index = $fields[1];
	my $bin = int($fragment_index / $bin_size);

	foreach my $split_info (@{$split_bins{$bin}})
	{
		if ($fragment_index >= $split_info->[0] and $fragment_index <= $split_info->[1])
		{
			my $filename = $split_info->[2];
			print $filename $line;
		}
	}
}
close CAN;

# Output the list of split canadidate filenames
print @split_candidates_filenames;

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

