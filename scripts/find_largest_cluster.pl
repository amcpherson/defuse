#!/usr/bin/perl

use strict;
use warnings;

my $clusters_filename = shift;

defined $clusters_filename or die "Usage: samtools view discordant.aligned.bam | $0 clusters.txt\n";

my $current_cluster_id;
my $current_cluster_size = 0;
my $current_cluster_align;
my $largest_cluster_id;
my $largest_cluster_size = -1;
my $largest_cluster_align;
open CLU, $clusters_filename or die "Error: Unable to open $clusters_filename: $!\n";
while (<CLU>)
{
	my ($cluster_id,$fragment_id,$align1_index,$align2_index) = split /\t/;
	
	if (defined $current_cluster_id and $current_cluster_id != $cluster_id)
	{
		if ($current_cluster_size > $largest_cluster_size)
		{
			$largest_cluster_id = $current_cluster_id;
			$largest_cluster_size = $current_cluster_size;
			$largest_cluster_align = $current_cluster_align;
		}

		$current_cluster_size = 0;
	}

	$current_cluster_id = $cluster_id;
	$current_cluster_size++;
	$current_cluster_align = [$align1_index,$align2_index];
}
close CLU;

if (defined $current_cluster_id)
{
	if ($current_cluster_size > $largest_cluster_size)
	{
		$largest_cluster_id = $current_cluster_id;
		$largest_cluster_size = $current_cluster_size;
		$largest_cluster_align = $current_cluster_align;
	}
}

print "Cluster: ".$largest_cluster_id."\n";
print "Size: ".$largest_cluster_size."\n";

my $current_align_index = 0;
while (<>)
{
	if ($current_align_index == $largest_cluster_align->[0] or $current_align_index == $largest_cluster_align->[1])
	{
		print $_;
	}
	
	$current_align_index++;
}


