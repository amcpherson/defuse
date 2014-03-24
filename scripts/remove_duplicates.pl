#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];


my $min_cluster_size = shift;

defined $min_cluster_size or die "Usage: $0 min_cluster_size < in_clusters > out_clusters\n";

sub remove_duplicates
{
	my %positions;
	my %lines;
	foreach my $line (@_)
	{
		chomp $line;
		my @fields = split /\t/, $line;
		
		my $cluster_id = $fields[0];
		my $cluster_end = $fields[1];
		my $fragment_id = $fields[2];
		my $read_end = $fields[3];
		my $ref_name = $fields[4];
		my $strand = $fields[5];
		my $start = $fields[6];
		my $end = $fields[7];

		my $position = ($strand eq "+") ? $start : $end;

		$positions{$fragment_id}{$cluster_end} = $position;
		$lines{$fragment_id}{$cluster_end} = $line;
	}

	my @undup;
	my %position_pairs;
	foreach my $fragment_id (keys %positions)
	{
		my $position1 = $positions{$fragment_id}{"0"};
		my $position2 = $positions{$fragment_id}{"1"};

		my $line1 = $lines{$fragment_id}{"0"};
		my $line2 = $lines{$fragment_id}{"1"};

		next if defined $position_pairs{$position1."-".$position2};

		$position_pairs{$position1."-".$position2} = 1;
		push @undup, $line1."\n";
		push @undup, $line2."\n";
	}

	return @undup;
}

my $current_id;
my @lines;
while (<>)
{
	my $line = $_;
	
	chomp;
	my @fields = split /\t/, $line;
	
	my $cluster_id = $fields[0];
	my $cluster_end = $fields[1];
	my $fragment_id = $fields[2];
	my $read_end = $fields[3];
	my $ref_name = $fields[4];
	my $strand = $fields[5];
	my $start = $fields[6];
	my $end = $fields[7];
	
	if (defined $current_id and $current_id != $cluster_id)
	{
		my @undup_lines = remove_duplicates(@lines);
		
		if (scalar @undup_lines >= 2 * $min_cluster_size)
		{
			print @undup_lines;
		}

		@lines = ();
	}
	
	$current_id = $cluster_id;
	
	push @lines, $line;
}

if (defined $current_id)
{
	my @undup_lines = remove_duplicates(@lines);
	
	if (scalar @undup_lines >= 2 * $min_cluster_size)
	{
		print @undup_lines;
	}
}


