#!/usr/bin/perl

use strict;
use warnings;

my $usage = "Usage: $0 sam1 sam2\n";

my $sam1 = shift;
my $sam2 = shift;

die $usage if not defined $sam2;

sub read_align
{
	my $sam_filename = shift;
	my $align_ref = shift;

	open SAM, $sam_filename or die "Error: Unable to open $sam_filename\n";
	while (<SAM>)
	{
		next if /^\@/;

		chomp;
		my @fields = split /\t/;

		my $readid = $fields[0];
		my $flag = $fields[1];
		my $reference = $fields[2];

		$readid =~ s/\/[12]$//;

		next if ($flag & hex('0x0004'));

		$reference =~ /(ENSG\d+)/;
		my $gene = $1;

		die "Error: Unrecognized reference name format\n" if not defined $gene;

		$align_ref->{$readid}{$gene} = 1;
	}
	close SAM;
}

my %align1;
my %align2;

read_align($sam1,\%align1);
read_align($sam2,\%align2);

foreach my $readid1 (keys %align1)
{
	my $concordant = 0;

	foreach my $gene1 (keys %{$align1{$readid1}})
	{
		if (defined $align2{$readid1}{$gene1})
		{
			$concordant = 1;
			last;
		}
	}

	if ($concordant)
	{
		print $readid1."\n";
	}
}

