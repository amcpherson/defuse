#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

my $prefix = shift;

die "Usage: $0 prefix\n" if not defined $prefix;

my $fastq1 = $prefix.".1.fastq";
my $fastq2 = $prefix.".2.fastq";

my $fastq1_tmp = $prefix.".1.fastq.tmp";
my $fastq2_tmp = $prefix.".2.fastq.tmp";

open FQ1, $fastq1 or die "Error: Unable to open $fastq1\n";
open FQ2, $fastq2 or die "Error: Unable to open $fastq2\n";

open TMP1, ">".$fastq1_tmp or die "Error: Unable to open $fastq1_tmp\n";
open TMP2, ">".$fastq2_tmp or die "Error: Unable to open $fastq2_tmp\n";

while (1)
{
        my $readid1 = <FQ1>;
        my $sequence1 = <FQ1>;
        my $comment1 = <FQ1>;
        my $quality1 = <FQ1>;

        last if not defined $quality1;

        chomp($readid1);
        chomp($sequence1);
        chomp($comment1);
        chomp($quality1);

        my $readid2 = <FQ2>;
        my $sequence2 = <FQ2>;
        my $comment2 = <FQ2>;
        my $quality2 = <FQ2>;

        last if not defined $quality2;

        chomp($readid2);
        chomp($sequence2);
        chomp($comment2);
        chomp($quality2);

        if (not is_garbage($sequence1) and not is_garbage($sequence2))
        {
                print TMP1 "$readid1\n$sequence1\n$comment1\n$quality1\n";
                print TMP2 "$readid2\n$sequence2\n$comment2\n$quality2\n";
        }
}

close FQ1;
close FQ2;

close TMP1;
close TMP2;

rename $fastq1_tmp, $fastq1;
rename $fastq2_tmp, $fastq2;

sub is_garbage
{
	my $sequence = $_[0];

	my $max_single_nt = 0.9 * length($sequence);
	my $max_run_nt = 0.5 * length($sequence);

	return 1 if @{[$sequence =~ /(A)/g]} > $max_single_nt;
	return 1 if @{[$sequence =~ /(C)/g]} > $max_single_nt;
	return 1 if @{[$sequence =~ /(T)/g]} > $max_single_nt;
	return 1 if @{[$sequence =~ /(G)/g]} > $max_single_nt;
	return 1 if @{[$sequence =~ /(N)/g]} > $max_single_nt;
	return 1 if $sequence =~ /(A{$max_run_nt,})/;
	return 1 if $sequence =~ /(C{$max_run_nt,})/;
	return 1 if $sequence =~ /(T{$max_run_nt,})/;
	return 1 if $sequence =~ /(G{$max_run_nt,})/;
	return 1 if $sequence =~ /(N{$max_run_nt,})/;

	return 0;
}

