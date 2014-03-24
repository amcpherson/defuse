#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';


# Count the number of paired end sequences aligning to each gene

my %gene_counts;

while (my $line = <>)
{
	chomp($line);
	
	next if $line =~ /^\@/;
	
	my @fields = split /\t/, $line;

	my $reference = $fields[2];
	$reference =~ /^([^|]*)/;
	my $gene = $1;

	$gene_counts{$gene} += 0.5;
}

foreach my $gene (keys %gene_counts)
{
	print $gene."\t".$gene_counts{$gene}."\n";
}

