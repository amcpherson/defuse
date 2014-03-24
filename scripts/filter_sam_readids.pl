#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;

my $invert;
GetOptions('invert' => \$invert);

my $readids_filename = shift;

my $usage = "Usage: $0 [-i|--invert] read_ids_filename\n";

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

	if ((defined $readids{$sam_readid} and not defined $invert) or (not defined $readids{$sam_readid} and defined $invert))
	{
		print $line."\n";
	}
}

