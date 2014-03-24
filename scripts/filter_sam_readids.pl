#!/usr/bin/perl

use strict;
use warnings FATAL => qw( all );

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

	next if $line =~ /^\@/;

	my @sam_fields = split /\t/, $line;

	my $sam_readid = $sam_fields[0];
	$sam_readid =~ s/\/[12]//;

	if ($readids{$sam_readid})
	{
		print $line."\n";
	}
}

