#!/usr/bin/perl

use strict;
use warnings;
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

my $fq_index = $prefix.".fqi";
my $fragment_names = $prefix.".names";

open FQ1, $fastq1 or die "Error: Unable to open $fastq1\n";
open FQ2, $fastq2 or die "Error: Unable to open $fastq2\n";

open TMP1, ">".$fastq1_tmp or die "Error: Unable to open $fastq1_tmp\n";
open TMP2, ">".$fastq2_tmp or die "Error: Unable to open $fastq2_tmp\n";

open FQI, ">".$fq_index or die "Error: Unable to open $fq_index\n";
binmode(FQI);

open NAM, ">".$fragment_names or die "Error: Unable to open $fragment_names\n";

my $current_fragment_index = 0;
while (1)
{
	my $filepos1 = pack('q',tell(TMP1));
	my $filepos2 = pack('q',tell(TMP2));
	
	my $readid1 = <FQ1>;
	my $sequence1 = <FQ1>;
	my $comment1 = <FQ1>;
	my $quality1 = <FQ1>;

	last if not defined $quality1;

	chomp($readid1);
	chomp($sequence1);
	chomp($comment1);
	chomp($quality1);
	
	$readid1 =~ /^\@(.*)\/([12])$/ or die "Fastq error\n";
	my $fragment1 = $1;
	my $end1 = $2;

	my $readid2 = <FQ2>;
	my $sequence2 = <FQ2>;
	my $comment2 = <FQ2>;
	my $quality2 = <FQ2>;

	last if not defined $quality2;

	chomp($readid2);
	chomp($sequence2);
	chomp($comment2);
	chomp($quality2);

	$readid2 =~ /^\@(.*)\/([12])$/ or die "Fastq error\n";
	my $fragment2 = $1;
	my $end2 = $2;
	
	die "Fastq error\n" if $fragment1 ne $fragment2 or $end1 ne '1' or $end2 ne '2';

	print FQI $filepos1;
	print FQI $filepos2;
	print TMP1 "\@$current_fragment_index/1\n$sequence1\n$comment1\n$quality1\n";
	print TMP2 "\@$current_fragment_index/2\n$sequence2\n$comment2\n$quality2\n";
	print NAM "$current_fragment_index\t$fragment1\n";
	
	$current_fragment_index++;
}

close FQ1;
close FQ2;

close TMP1;
close TMP2;

close FQI;

close NAM;

rename $fastq1_tmp, $fastq1;
rename $fastq2_tmp, $fastq2;

