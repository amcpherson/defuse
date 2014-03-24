#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

use lib dirname($0)."/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my %cluster_strand;
my %cluster_align_start;
my %cluster_align_end;
while (<>)
{
	chomp;
	my @fields = split /\t/;
	
	my $cluster_id = $fields[0];
	my $read_id = $fields[1];
	my $flag = $fields[2];
	my $ref_name = $fields[3];
	my $start = $fields[4];
	my $end = $start + length($fields[10]) - 1;

	$read_id =~ /(.*)\/([12])/;
	my $fragment_id = $1;
	my $read_end = $2;

	my $strand;
	if ($flag & hex('0x0010'))
	{
		$strand = "-";
	}
	else
	{
		$strand = "+";
	}

	$cluster_strand{$cluster_id}{$ref_name} = $strand;
	
	$cluster_align_start{$cluster_id}{$ref_name} = $start if not defined $cluster_align_start{$cluster_id}{$ref_name};
	$cluster_align_end{$cluster_id}{$ref_name} = $end if not defined $cluster_align_end{$cluster_id}{$ref_name};
	
	$cluster_align_start{$cluster_id}{$ref_name} = min($cluster_align_start{$cluster_id}{$ref_name}, $start);
	$cluster_align_end{$cluster_id}{$ref_name} = max($cluster_align_end{$cluster_id}{$ref_name}, $end);
}

foreach my $cluster_id (keys %cluster_strand)
{
	die "Error: Did not find 2 reference names for cluster $cluster_id\n" if scalar keys %{$cluster_strand{$cluster_id}} != 2;

	foreach my $ref_name (keys %{$cluster_strand{$cluster_id}})
	{
		print $cluster_id."\t";
		print $ref_name."\t";
		print $cluster_strand{$cluster_id}{$ref_name}."\t";
		print $cluster_align_start{$cluster_id}{$ref_name}."\t";
		print $cluster_align_end{$cluster_id}{$ref_name}."\n";
	} 
}


