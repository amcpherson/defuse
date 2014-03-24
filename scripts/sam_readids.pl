#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

while (<>)
{
	chomp;
	my @sam_fields = split /\t/;
	
	next if $sam_fields[0] =~ /^\@/g;
	
	$sam_fields[0] =~ s/\/[12]//;
	print $sam_fields[0]."\n";
}


