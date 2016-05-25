#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin";
use gene_models;

my $gene_models_filename = shift;
my $region_length = shift;

defined $region_length or die "Usage: $0 gene_models_filename region_length < sam_alignments > concordant_readids\n";

my $gene_models = gene_models->new($gene_models_filename);

my $extend_length = int($region_length / 2);
my $bin_spacing = int($region_length / 2);

my %align_bins;
while (<>)
{
	next if /^\@/;

	chomp;
	my @fields = split /\t/;

	my $flag = $fields[1];

	next if ($flag & hex('0x0004'));

	my $readid = $fields[0];
	my $reference = $fields[2];
	my $pos = $fields[3];
	my $seq = $fields[9];
	
	$readid =~ /^(.*)\/([12])$/ or die "Fastq error, unable to interpret readid $readid\n";
	my $fragment = $1;
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
	
	my $start = $pos - $extend_length;
	my $end = $pos + length($seq) - 1 + $extend_length;
	
	my $chromosome = $gene_models->calc_genomic_chromosome($reference);
	my @genomic_regions = $gene_models->calc_genomic_regions($reference, [$start,$end]);
	
	foreach my $bin (gene_models::get_bins([$genomic_regions[0]->[0],$genomic_regions[$#genomic_regions]->[1]], $bin_spacing))
	{
		$align_bins{$fragment}{$read_end}{$chromosome}{$bin} = 1;
	}
}

foreach my $fragment (keys %align_bins)
{
	my $concordant = 0;
	foreach my $chromosome (keys %{$align_bins{$fragment}{"1"}})
	{
		foreach my $bin (keys %{$align_bins{$fragment}{"1"}{$chromosome}})
		{
			if (defined $align_bins{$fragment}{"2"}{$chromosome}{$bin})
			{
				$concordant = 1;
				last;
			}
		}
		
		if ($concordant)
		{
			last;
		}
	}
	
	if ($concordant)
	{
		print $fragment."\n";
	}
}

