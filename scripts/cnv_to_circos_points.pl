#!/usr/bin/perl

use strict;
use warnings;

my $level_list = shift;

defined $level_list or die "Usage: $0 levels\n";

my %levels = map { $_ => 1 } split /,/, $level_list;

while (<>)
{
	chomp;
	
	my ($chr,$start,$end,$level,$value) = split /\t/;
	
	$chr =~ s/23/X/;
	$chr = "hs".$chr;
	
	if (defined $levels{$level})
	{
		print $chr."\t".$start."\t".$end."\t".$value."\n";
	}
}