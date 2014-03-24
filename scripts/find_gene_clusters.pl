#!/usr/bin/perl

use strict;
use warnings;

my %mappings;
while (<>)
{
	my $line = $_;

	chomp;
	next if /^\@/;
	
	my @fields = split /\t/;

	my $sam_readid = $fields[0];
	my $reference = $fields[2];
	
	$sam_readid =~ s/(.*)\/([12])//;
	my $fragment = $1;
	my $read_end = $2;

	$reference =~ /(ENSG\d+)/;
	my $gene = $1;
	
	$mappings{$fragment}{$read_end}{$gene} = 1;
}

my %fusions;
foreach my $fragment (keys %mappings)
{
	foreach my $gene1 (keys %{$mappings{$fragment}{"1"}})
	{
		foreach my $gene2 (keys %{$mappings{$fragment}{"2"}})
		{
			($gene1,$gene2) = sort ($gene1,$gene2);
			$fusions{$gene1."\t".$gene2}++;
		}
	}
}

my @fusions_sorted = sort { $fusions{$a} <=> $fusions{$b} } keys %fusions;

foreach my $fusion (@fusions_sorted)
{
	print $fusion."\t".$fusions{$fusion}."\n";
}
