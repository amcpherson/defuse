#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
	

my $percentIdentityThreshold = shift;

die "Usage: $0 percentIdentityThreshold\n" if not defined $percentIdentityThreshold;

while (<>)
{
	my $line = $_;

	chomp;
	my @pslFields = split /\t/;

	next unless scalar(@pslFields) >= 11;

	my $numMatches = $pslFields[0];
	my $breakpointSeqLength = $pslFields[10];

	next unless $numMatches =~ /^\d+$/;

	my $percentIdentity = $numMatches / $breakpointSeqLength;

	print $line if $percentIdentity > $percentIdentityThreshold;
}


