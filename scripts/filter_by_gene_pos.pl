#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $distlimit = 50000;

while (<>)
{
	chomp;
	my @fields = split /\t/;

	next if ($fields[4] eq $fields[12] and (abs($fields[7] - $fields[14]) < $distlimit or abs($fields[6] - $fields[15]) < $distlimit));

	print $_."\n";
}


