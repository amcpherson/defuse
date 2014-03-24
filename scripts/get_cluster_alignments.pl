#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $clusters_filename = shift;

die "Usage: $0 clusters_filename < in.sam > out.clusters.sam\n" if not defined $clusters_filename;

my %align_index_clusters;
open CLU, $clusters_filename or die "Error: Unable to open $clusters_filename: $!\n";
while (<CLU>)
{
	chomp;
	my ($cluster_id, $cluster_end, $fragment_index, $align_index) = split /\t/;
	
	push @{$align_index_clusters{$align_index}}, [$cluster_id,$cluster_end];
}
close CLU;

my $current_align_index = 0;
while (<>)
{
	if (defined $align_index_clusters{$current_align_index})
	{
		foreach my $cluster (@{$align_index_clusters{$current_align_index}})
		{
			print $cluster->[0]."\t".$cluster->[1]."\t".$_;
		}
	}
	
	$current_align_index++;
}
