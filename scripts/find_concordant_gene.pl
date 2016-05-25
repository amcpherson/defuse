#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin";
use gene_models;

my $gene_models_filename = shift;

defined $gene_models_filename or die "Usage: $0 gene_models_filename < sam_alignments > concordant_readids\n";

my $gene_models = gene_models->new($gene_models_filename);

my %align_genes;
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
	
	my $start = $pos;
	my $end = $pos + length($seq) - 1;
	
	foreach my $gene ($gene_models->calc_overlapping_genes($reference, [$start,$end]))
	{
		$align_genes{$fragment}{$read_end}{$gene} = 1;
	}
}

foreach my $fragment (keys %align_genes)
{
	my $concordant = 0;
	foreach my $gene (keys %{$align_genes{$fragment}{"1"}})
	{
		if (defined $align_genes{$fragment}{"2"}{$gene})
		{
			$concordant = 1;
			last;
		}
	}
	
	if ($concordant)
	{
		print $fragment."\n";
	}
}

