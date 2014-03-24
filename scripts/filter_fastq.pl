#!/usr/bin/perl

use strict;
use warnings;
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

	$readid =~ /^@(.+)\/[12]/;

	if ((defined $readids{$1} and not defined $invert) or (not defined $readids{$1} and defined $invert))
	{
		print "$readid\n$sequence\n$comment\n$quality\n";
	}
}

