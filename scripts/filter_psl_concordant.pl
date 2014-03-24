#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
	

my $min_percent_identity = shift;
my $max_align_region = shift;

die "Usage: $0 min_percent_identity max_align_region\n" if not defined $max_align_region;

while (<>)
{
	my $line = $_;

	chomp;
	my @psl_fields = split /\t/;

	next if scalar @psl_fields < 1 or not $psl_fields[0] =~ /^\d+$/;

	my $num_matches = $psl_fields[0];
	my $num_target_bases_inserted = $psl_fields[7];
	my $strand = $psl_fields[8];
	my $cluster_id = $psl_fields[9];
	my $query_size = $psl_fields[10];
	my $query_start = $psl_fields[11] + 1;
	my $query_end = $psl_fields[12];
	my $target_seq_name = $psl_fields[13];
	my $target_size = $psl_fields[14];
	my $target_start = $psl_fields[15] + 1;
	my $target_end = $psl_fields[16];
	my @block_sizes = split /,/, $psl_fields[18];
	my @block_query_starts = split /,/, $psl_fields[19];
	my @block_target_starts = split /,/, $psl_fields[20];

	my $percent_identity = $num_matches / $query_size;
	my $align_region = $target_end - $target_start + 1;

	next if $percent_identity < $min_percent_identity;
	next if $align_region > $max_align_region;

	print $line;
}


