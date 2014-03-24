#!/usr/bin/perl

use strict;
use warnings;

my $list1_filename = shift;
my $list2_filename = shift;

defined $list2_filename or die "Usage: $0 list1_filename list2_filename\n";

my %list1;
open LIST1, $list1_filename or die "Error: Unable to open $list1_filename $!\n";
while (<LIST1>)
{
	chomp;
	$list1{$_} = 1;
}
close LIST1;

open LIST2, $list2_filename or die "Error: Unable to open $list2_filename $!\n";
while (<LIST2>)
{
	chomp;
	
	if (defined $list1{$_})
	{
		print $_."\n";
	}
}
close LIST2;

