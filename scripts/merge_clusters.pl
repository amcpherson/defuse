#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use List::Util qw(min max);

die "Usage $0 clusters1 clusters2 ...\n" if scalar @ARGV == 0;

my $cluster_id = 0;
foreach my $cluster_filename (@ARGV)
{
	my $previous_cluster_id;
	open CLS, $cluster_filename or die "Error: Unable to open $cluster_filename\n";
	while (<CLS>)
	{
		my @fields = split /\t/;
		
		if (defined $previous_cluster_id and $previous_cluster_id != $fields[0])
		{
			$cluster_id++;
		}
		$previous_cluster_id = $fields[0];
		
		$fields[0] = $cluster_id;
		
		print join "\t", @fields;
	}
	close CLS;

	if (defined $previous_cluster_id)
	{
		$cluster_id++;
	}
}

