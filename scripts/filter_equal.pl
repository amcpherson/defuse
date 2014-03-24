#!/usr/bin/perl

use strict;
use warnings;

my $column = shift;
my $yn = shift;

die "Usage: $0 column [Y|N]\n" if not defined $yn;

my $found_header = 0;
my %ind;
while (<>)
{
	chomp;
	my @fields = split /\t/;

	if (not $found_header)
	{
		foreach my $field_index (0..$#fields)
		{
			$ind{$fields[$field_index]} = $field_index;
		}

		$found_header = 1;

		print $_."\n";
		next;
	}

	if ($found_header and $fields[0] eq "cluster_id")
	{
		next;
	}

	if ($ind{$column} >= scalar @fields)
	{
		die "Error: invalid field for line:\n$_";
	}
	
	my $value = $fields[$ind{$column}];
	
	next unless $value eq $yn;
	
	print $_."\n";
}
