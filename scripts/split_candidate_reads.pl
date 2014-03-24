#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
no strict 'refs';

my $usage = "Usage: $0 split_catalog candidate_reads > candidate_catalog\n";

my $split_catalog_filename = shift;
my $candidate_reads_filename = shift;

die $usage if not defined $candidate_reads_filename;

# Read in split catalog
my @split_catalog = read_split_catalog($split_catalog_filename);

# Maximum cached entries
my $max_cached = 10000;

# Open files for writing
# Create binning for fragment indices
my $bin_size = 100000;
my %split_bins;
my @split_candidates_filenames;
my %split_candidates;
foreach my $split_info (@split_catalog)
{
	my $split_prefix = $split_info->[0];
	my $start_fragment_index = $split_info->[1];
	my $end_fragment_index = $split_info->[2];

	my $split_candidates_filename = $split_prefix.".candidate.reads";
	
	my $start_bin = int($start_fragment_index / $bin_size);
	my $end_bin = int($end_fragment_index / $bin_size);

	foreach my $bin ($start_bin..$end_bin)
	{
		push @{$split_bins{$bin}}, [$start_fragment_index, $end_fragment_index, $split_candidates_filename];
	}

	push @split_candidates_filenames, $split_candidates_filename;

	unlink $split_candidates_filename;
}

sub flush
{
	my $filename = shift;
	my $cached_ref = shift;
	
	open OUT, ">>".$filename or die "Error: Unable to write to $filename\n $!\n";
	foreach my $line (@{$cached_ref})
	{
		print OUT $line;
	}
	close OUT;
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
			push @{$split_candidates{$filename}}, $line;
			
			if (scalar @{$split_candidates{$filename}} >= $max_cached)
			{
				flush($filename, \@{$split_candidates{$filename}});
				@{$split_candidates{$filename}} = ();
			}

			last;
		}
	}
}
close CAN;

foreach my $filename (@split_candidates_filenames)
{
	flush($filename, \@{$split_candidates{$filename}});
}

# Output the list of split canadidate filenames
foreach my $filename (@split_candidates_filenames)
{
	print $filename."\n";
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

