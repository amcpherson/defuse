#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $genes_filename = shift;

die "Usage: $0 genes\n" if not defined $genes_filename;

my %genes;

open GENES, $genes_filename or die "Error: Unable to open $genes_filename: $!\n";
while (my $gene = <GENES>)
{
	chomp($gene);
	$genes{$gene} = 1;	
}
close GENES;

while (<>)
{
	my $line = $_;

	chomp;
	next if /^\@/;
	
	my @fields = split /\t/;

	my $reference = $fields[2];
	$reference =~ /^([^|]*)/;
	my $gene = $1;

	if ($genes{$gene})
	{
		print $line;
	}
}

