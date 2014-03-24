#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 read_ids_filename\n";

my $readids_filename = shift;

die $usage if not defined $readids_filename;

my %readids;

open RIDS, $readids_filename or die "Error: Unable to open $readids_filename: $!\n";
while (my $readid = <RIDS>)
{
	chomp($readid);
	$readids{$readid} = 1;	
}
close RIDS;

while (my $line = <>)
{
	chomp($line);

	my @bowtie_fields = split /\t/, $line;

	my $readid = $bowtie_fields[0];

	$readid =~ /^(.+)\/[12]/;

	if ($readids{$1})
	{
		print $line."\n";
	}
}

