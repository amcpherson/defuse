#!/usr/bin/perl -w

# Author: lh3
# Version: 0.1.7

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;

my $all_flag;
my $help_flag;

GetOptions("all" => \$all_flag, "help" => \$help_flag);

my $usage = "$0 [-a (use all reads)] < in.qseq";

if (defined $help_flag)
{
	print $usage."\n";
	exit;
}

# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64)
{
	$conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}

while (<STDIN>)
{
	chomp;
	my @fields = split /\t/;

	die "Error: not enough fields\n" if scalar @fields < 11;

	next if not defined $all_flag and $fields[10] == 0;

	my $name = "\@".$fields[0]."_".$fields[1].":".$fields[2].":".$fields[3].":".$fields[4].":".$fields[5]."/".$fields[7];

	my $sequence = $fields[8];
	$sequence =~ tr/\./N/;

	my $quality = "";
	foreach my $basequal (split '',$fields[9])
	{
		$quality .= $conv_table[ord($basequal)];
	}

	print $name."\n".$sequence."\n+\n".$quality."\n";
}

