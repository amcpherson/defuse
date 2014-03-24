#!/usr/bin/perl

use strict;
use warnings;

my $clusters_filename = shift;

die "Usage: $0 clusters_filename < in.sam > out.clusters.sam\n" if not defined $clusters_filename;

my %align_index_clusters;
open CLU, $clusters_filename or die "Error: Unable to open $clusters_filename: $!\n";
while (<CLU>)
{
	chomp;
	my ($cluster_id, $fragment_index, $align_index_1, $align_index_2) = split /\t/;
	
	push @{$align_index_clusters{$align_index_1}}, $cluster_id;
	push @{$align_index_clusters{$align_index_2}}, $cluster_id;
}
close CLU;

my $current_align_index = 0;
while (<>)
{
	if (defined $align_index_clusters{$current_align_index})
	{
		foreach my $cluster_id (@{$align_index_clusters{$current_align_index}})
		{
			print $cluster_id."\t".$_;
		}
	}
	
	$current_align_index++;
}
