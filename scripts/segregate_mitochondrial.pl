#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin";
use gene_models;

my $gene_models_filename = shift;
my $mt_chromosome = shift;

defined $mt_chromosome or die "Usage: $0 gene_models_filename mt_chromosome\n";

my $gene_models = gene_models->new($gene_models_filename);

my $current_id;
my %chromosomes;
my @lines;
while (<>)
{
	my $line = $_;
	
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
	
	if (defined $current_id and $current_id != $cluster_id)
	{
		die "Error: Unmatched cluster $current_id\n" if not defined $chromosomes{"0"} or not defined $chromosomes{"1"};
		
		if (($chromosomes{"0"} ne $mt_chromosome and $chromosomes{"1"} ne $mt_chromosome) or $chromosomes{"0"} eq $chromosomes{"1"})
		{
			print @lines;
		}
		
		%chromosomes = ();
		@lines = ();
	}
	
	my $chromosome = $gene_models->calc_genomic_chromosome($ref_name);
	
	$current_id = $cluster_id;
	$chromosomes{$cluster_end} = $chromosome;
	push @lines, $line;
}

if (defined $current_id)
{
	die "Error: Unmatched cluster $current_id\n" if not defined $chromosomes{"0"} or not defined $chromosomes{"1"};
	
	if (($chromosomes{"0"} ne $mt_chromosome and $chromosomes{"1"} ne $mt_chromosome) or $chromosomes{"0"} eq $chromosomes{"1"})
	{
		print @lines;
	}
}
