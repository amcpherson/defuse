#!/usr/bin/perl

use strict;
use warnings;

my $seqnum = 1;
while (<>)
{
	chomp;
	my @fields = split;
	die if scalar @fields < 1;
	if (scalar @fields == 1)
	{
		print ">".$seqnum."\n".$fields[0]."\n";
		$seqnum++;	
	}
	else
	{
		my $seq = pop @fields;
		my $seqname = join "_", @fields;
		print ">".$seqname."\n".$seq."\n";
	}
}
