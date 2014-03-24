#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw[min max];

my $usage = "Usage: $0 trim_length\n";

my $trim_length = shift;

die $usage if not defined $trim_length;

while (1)
{
	my $readid = <>;
	my $sequence = <>;
	my $comment = <>;
	my $quality = <>;

	last if not defined $quality;

	chomp($readid);
	chomp($sequence);
	chomp($comment);
	chomp($quality);
	
	my $read_trim_length = min(length($sequence), $trim_length);

	$sequence = substr($sequence, 0, $read_trim_length);
	$quality = substr($quality, 0, $read_trim_length);

	print "$readid\n$sequence\n$comment\n$quality\n";
}

