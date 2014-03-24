#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 references\n";

my $references_filename = shift;

die $usage if not defined $references_filename;

my %references;

open REFS, $references_filename or die "Error: Unable to open $references_filename: $!\n";
while (my $reference = <REFS>)
{
	chomp($reference);
	$references{$reference} = 1;	
}
close REFS;

while (my $line = <>)
{
	chomp($line);

	my @bowtie_fields = split /\t/, $line;

	my $reference = $bowtie_fields[2];

	if ($references{$reference})
	{
		print $line."\n";
	}
}

