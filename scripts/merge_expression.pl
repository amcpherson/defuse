#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use List::Util qw(min max);

die "Usage $0 expressio1 expression2 ...\n" if scalar @ARGV == 0;

my %expression;

foreach my $expression_filename (@ARGV)
{
	open EXP, $expression_filename or die "Error: Unable to open $expression_filename\n";
	my @expression_entries = <EXP>;
	chomp(@expression_entries);
	close EXP;
	
	foreach my $expression_entry (@expression_entries)
	{
		my ($gene,$count) = split /\t/, $expression_entry;
		next if not defined $gene or not defined $count;
		$expression{$gene} += $count;
	}
}

foreach my $gene (keys %expression)
{
	print $gene."\t".$expression{$gene}."\n";
}

