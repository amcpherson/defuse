#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Temp qw[tempfile tempdir];

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Calculate spanning fragment length pvalue\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --clusters  Spanning Clusters Sam Filename\n";
push @usage, "  -b, --breaks    Breakpoint Predictions Filename\n";
push @usage, "  -s, --seqs      Sequence Predictions Filename\n";

my $help;
my $clusters_filename;
my $breaks_filename;
my $seqs_filename;

GetOptions
(
 	'help'        => \$help,
	'clusters=s'  => \$clusters_filename,
 	'breaks=s'    => \$breaks_filename,
	'seqs=s'      => \$seqs_filename,
);

not defined $help or usage() and exit;

defined $clusters_filename or die @usage;
defined $breaks_filename or die @usage;
defined $seqs_filename or die @usage;

my %cluster_break_pos;
open BR, $breaks_filename or die "Error: Unable to find $breaks_filename: $!\n";
while (<BR>)
{
	chomp;
	my @fields = split /\t/;
	
	my $cluster_id = $fields[0];
	my $cluster_end = $fields[1];
	my $ref_name = $fields[2];
	my $strand = $fields[3];
	my $break_pos = $fields[4];
	
	$cluster_break_pos{$cluster_id}{$cluster_end} = $break_pos;
}
close BR;

my %cluster_inter_length;
open SEQ, $seqs_filename or die "Error: Unable to find $seqs_filename: $!\n";
while (<SEQ>)
{
	chomp;
	my @fields = split /\t/;
	
	my $cluster_id = $fields[0];
	my $inter_length = $fields[2];
	
	$cluster_inter_length{$cluster_id} = $inter_length;
}
close SEQ;

my %cluster_strand;
my %cluster_align_start;
my %cluster_align_end;
open SP, $clusters_filename or die "Error: Unable to find $clusters_filename: $!\n";
while (<SP>)
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
	
	$cluster_strand{$cluster_id}{$cluster_end} = $strand;
	
	$cluster_align_start{$cluster_id}{$fragment_id}{$cluster_end} = $start;
	$cluster_align_end{$cluster_id}{$fragment_id}{$cluster_end} = $end;
}
close SP;

foreach my $cluster_id (keys %cluster_strand)
{
	next if not defined $cluster_break_pos{$cluster_id};
	
	die "Error: Did not find 2 reference names for cluster $cluster_id\n" if scalar keys %{$cluster_strand{$cluster_id}} != 2;

	my $spanning_fragment_sum = 0;
	my $spanning_fragment_count = 0;

	foreach my $fragment_id (keys %{$cluster_align_start{$cluster_id}})
	{
		my $fragment_length = 0;
		foreach my $cluster_end (keys %{$cluster_align_start{$cluster_id}{$fragment_id}})
		{
			my $strand = $cluster_strand{$cluster_id}{$cluster_end};
			if ($strand eq "+")
			{
				$fragment_length += $cluster_break_pos{$cluster_id}{$cluster_end} - $cluster_align_start{$cluster_id}{$fragment_id}{$cluster_end} + 1;
			}
			else
			{
				$fragment_length += $cluster_align_end{$cluster_id}{$fragment_id}{$cluster_end} - $cluster_break_pos{$cluster_id}{$cluster_end} + 1;
			}
		}
		$fragment_length += $cluster_inter_length{$cluster_id};

		$spanning_fragment_sum += $fragment_length;
		$spanning_fragment_count++;
	}

	my $spanning_fragment_mean = $spanning_fragment_sum / $spanning_fragment_count;
	
	print $cluster_id."\t".$spanning_fragment_mean."\t".$spanning_fragment_count."\n";
}


