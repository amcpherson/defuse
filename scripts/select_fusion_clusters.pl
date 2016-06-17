#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin";
use gene_models;

my $gene_models_filename = shift;

defined $gene_models_filename or die "Usage: $0 gene_models_filename\n";

my $gene_models = gene_models->new($gene_models_filename);

sub is_fusion
{
	my $cluster_ref = shift;
	
	my $pos1 = ($cluster_ref->{"0"}{start} + $cluster_ref->{"0"}{end}) / 2;
	my $pos2 = ($cluster_ref->{"1"}{start} + $cluster_ref->{"1"}{end}) / 2;
	
	my $gene1 = $gene_models->calc_gene($cluster_ref->{"0"}{ref_name}, $pos1);
	my $gene2 = $gene_models->calc_gene($cluster_ref->{"1"}{ref_name}, $pos2);
	
	my $genomic_pos1 = $gene_models->calc_genomic_position($cluster_ref->{"0"}{ref_name}, $pos1);
	my $genomic_pos2 = $gene_models->calc_genomic_position($cluster_ref->{"1"}{ref_name}, $pos2);
	
	my $gene_location1 = $gene_models->calc_gene_location($gene1, $genomic_pos1);
	my $gene_location2 = $gene_models->calc_gene_location($gene2, $genomic_pos2);
	
	# Filter clusters involving the same gene
	return 0 if $gene1 eq $gene2;
	
	# Filter clusters that are intergenic for both sides
	my $intergenic1 = ($gene_location1 eq "upstream" or $gene_location1 eq "downstream");
	my $intergenic2 = ($gene_location2 eq "upstream" or $gene_location2 eq "downstream");
	return 0 if $intergenic1 and $intergenic2;
	
	return 1;
}


my $current_id;
my %cluster;
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
		die "Error: Unmatched cluster $current_id\n" if not defined $cluster{"0"} or not defined $cluster{"1"};
		
		if (is_fusion(\%cluster))
		{
			print @lines;
		}
		
		%cluster = ();
		@lines = ();
	}
	
	my $chromosome = $gene_models->calc_genomic_chromosome($ref_name);
	
	$current_id = $cluster_id;
	$cluster{$cluster_end}{ref_name} = $ref_name;
	
	$cluster{$cluster_end}{start} = $start if not defined $cluster{$cluster_end}{start};
	$cluster{$cluster_end}{end} = $end if not defined $cluster{$cluster_end}{end};
	
	$cluster{$cluster_end}{start} = min($cluster{$cluster_end}{start}, $start);
	$cluster{$cluster_end}{end} = max($cluster{$cluster_end}{end}, $end);
	
	push @lines, $line;
}

if (defined $current_id)
{
	die "Error: Unmatched cluster $current_id\n" if not defined $cluster{"0"} or not defined $cluster{"1"};
	
	if (is_fusion(\%cluster))
	{
		print @lines;
	}
}
