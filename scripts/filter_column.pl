#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $values_filename = shift;
my $column_index = shift;
my $invert = shift;

die "Usage: $0 values_filename column_index invert\n" if not defined $column_index;

$invert = 0 if not defined $invert;

my %values;

open VLS, $values_filename or die "Error: Unable to open $values_filename: $!\n";
while (my $line = <VLS>)
{
	chomp($line);
	
	my @fields = split /\t/, $line;
	
	my $value = $fields[0];
	
	$values{$value} = 1;	
}
close VLS;

while (my $line = <>)
{
	chomp($line);

	my @fields = split /\t/, $line;
	
	die "Error: invalid column" if $column_index >= scalar @fields;

	my $value = $fields[$column_index];

	if ((not $invert and defined $values{$value}) or ($invert and not defined $values{$value}))
	{
		print $line."\n";
	}
}

