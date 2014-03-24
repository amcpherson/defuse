#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use List::Util qw(min max);

die "Usage $0 read_stats1 read_stats2 ...\n" if scalar @ARGV == 0;

sub get_stats
{
	my $stats_filename = shift;
	my $stats_outref = shift;
	
	open STATS, $stats_filename or die "Error: Unable to open $stats_filename\n";
	my @stats = <STATS>;
	chomp(@stats);
	close STATS;

	scalar @stats == 2 or die "Error: Stats file $stats_filename does not have 2 lines\n";

	my @keys = split /\t/, $stats[0];
	my @values = split /\t/, $stats[1];
	
	while (scalar @keys > scalar @values)
	{
		push @values, "";
	}

	scalar @keys == scalar @values or die "Error: Stats file $stats_filename with column mismatch\n";

	foreach my $stat_index (0..$#keys)
	{
		my $key = $keys[$stat_index];
		my $value = $values[$stat_index];

		$stats_outref->{$key} = $value;
	}
}

my @all_read_stats;

foreach my $read_stats_filename (@ARGV)
{
	my %read_stats;
	get_stats($read_stats_filename, \%read_stats);

	die "Error: $read_stats_filename is incomplete" if not defined $read_stats{"frag_count"};
	die "Error: $read_stats_filename is incomplete" if not defined $read_stats{"fraglength_mean"};
	die "Error: $read_stats_filename is incomplete" if not defined $read_stats{"fraglength_stddev"};
	die "Error: $read_stats_filename is incomplete" if not defined $read_stats{"readlength_min"};
	die "Error: $read_stats_filename is incomplete" if not defined $read_stats{"readlength_max"};

	push @all_read_stats, \%read_stats;
}

my $fraglength_num = 0;
my $fraglength_sum = 0;
my $fraglength_sum_sq = 0;
my @read_lengths;

foreach my $read_stats (@all_read_stats)
{
	my $current_num_fragments = $read_stats->{"frag_count"};
	my $current_fraglength_sum = $read_stats->{"fraglength_mean"} * $read_stats->{"frag_count"};
	my $current_fraglength_variance = $read_stats->{"fraglength_stddev"} ** 2;
	my $current_fraglength_sum_sq = ($current_fraglength_variance + $read_stats->{"fraglength_mean"} ** 2) * $current_num_fragments;
	
	next unless $current_num_fragments > 0;
	
	$fraglength_num += $current_num_fragments;
	$fraglength_sum += $current_fraglength_sum;
	$fraglength_sum_sq += $current_fraglength_sum_sq;
	
	push @read_lengths, $read_stats->{"readlength_min"};
	push @read_lengths, $read_stats->{"readlength_max"};
}

my $fraglength_mean = $fraglength_sum / $fraglength_num;
my $fraglength_variance = $fraglength_sum_sq / $fraglength_num - $fraglength_mean ** 2;
my $fraglength_stddev = $fraglength_variance ** 0.5;
my $readlength_min = min(@read_lengths);
my $readlength_max = max(@read_lengths);

print "frag_count\tfraglength_mean\tfraglength_stddev\treadlength_min\treadlength_max\n";
print "$fraglength_num\t$fraglength_mean\t$fraglength_stddev\t$readlength_min\t$readlength_max\n";


