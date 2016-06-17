#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

my %cluster_ref_name;
my %cluster_strand;
my %cluster_align_start;
my %cluster_align_end;
while (<>)
{
	chomp;
	my @fields = split /\t/;
	
	my $cluster_id = $fields[0];
	my $cluster_end = $fields[1];
	my $fragment_id = $fields[2];
	my $read_end = $fields[3];
	my $ref_name = $fields[4];
	my $strand = $fields[5];
	my $start = $fields[6];
	my $end = $fields[7];
	
	$cluster_ref_name{$cluster_id}{$cluster_end} = $ref_name;
	$cluster_strand{$cluster_id}{$cluster_end} = $strand;
	
	$cluster_align_start{$cluster_id}{$cluster_end} = $start if not defined $cluster_align_start{$cluster_id}{$cluster_end};
	$cluster_align_end{$cluster_id}{$cluster_end} = $end if not defined $cluster_align_end{$cluster_id}{$cluster_end};
	
	$cluster_align_start{$cluster_id}{$cluster_end} = min($cluster_align_start{$cluster_id}{$cluster_end}, $start);
	$cluster_align_end{$cluster_id}{$cluster_end} = max($cluster_align_end{$cluster_id}{$cluster_end}, $end);
}

foreach my $cluster_id (keys %cluster_strand)
{
	die "Error: Did not find 2 ends for cluster $cluster_id\n" if scalar keys %{$cluster_strand{$cluster_id}} != 2;

	foreach my $cluster_end (keys %{$cluster_strand{$cluster_id}})
	{
		print $cluster_id."\t";
		print $cluster_end."\t";
		print $cluster_ref_name{$cluster_id}{$cluster_end}."\t";
		print $cluster_strand{$cluster_id}{$cluster_end}."\t";
		print $cluster_align_start{$cluster_id}{$cluster_end}."\t";
		print $cluster_align_end{$cluster_id}{$cluster_end}."\n";
	} 
}


