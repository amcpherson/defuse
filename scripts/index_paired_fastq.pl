#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

my $input_reads_end_1_fastq = shift;
my $input_reads_end_2_fastq = shift;
my $output_reads_end_1_fastq = shift;
my $output_reads_end_2_fastq = shift;
my $output_reads_index_filename = shift;
my $output_reads_names_filename = shift;

die "Usage: $0 infastq1 infastq2 outfastq1 outfastq2 outindex outnames\n" if not defined $output_reads_names_filename;

# Read from both input fastq files
open IN1, $input_reads_end_1_fastq or die "Error: Unable to write to $input_reads_end_1_fastq: $!\n";
open IN2, $input_reads_end_2_fastq or die "Error: Unable to write to $input_reads_end_2_fastq: $!\n";

# Open both fastq files for simultaneous writing
open FQ1, ">".$output_reads_end_1_fastq or die "Error: Unable to write to $output_reads_end_1_fastq: $!\n";
open FQ2, ">".$output_reads_end_2_fastq or die "Error: Unable to write to $output_reads_end_2_fastq: $!\n";
open FQI, ">".$output_reads_index_filename or die "Error: Unable to open $output_reads_index_filename\n"; binmode(FQI);
open NAM, ">".$output_reads_names_filename or die "Error: Unable to open $output_reads_names_filename\n";

my $current_fragment_index = 0;
while (1)
{
	my $readid1 = <IN1>;
	my $sequence1 = <IN1>;
	my $comment1 = <IN1>;
	my $quality1 = <IN1>;

	last if not defined $quality1;

	chomp($readid1);
	chomp($sequence1);
	chomp($comment1);
	chomp($quality1);
	
	my $readid2 = <IN2>;
	my $sequence2 = <IN2>;
	my $comment2 = <IN2>;
	my $quality2 = <IN2>;

	last if not defined $quality2;

	chomp($readid2);
	chomp($sequence2);
	chomp($comment2);
	chomp($quality2);

	my $filepos1 = pack('q',tell(FQ1));
	my $filepos2 = pack('q',tell(FQ2));
	
	print FQI $filepos1;
	print FQI $filepos2;
	print FQ1 "\@$current_fragment_index/1\n$sequence1\n$comment1\n$quality1\n";
	print FQ2 "\@$current_fragment_index/2\n$sequence2\n$comment2\n$quality2\n";
	print NAM "$current_fragment_index\t$readid1\t$readid2\n";
	
	$current_fragment_index++;
}

close IN1;
close IN2;

# Check if we found reads
die "Error: No reads found\n" if $current_fragment_index == 0;
print "Imported $current_fragment_index reads\n";

# Finished writing fastq sequences
close FQ1;
close FQ2;

close FQI;

close NAM;

print "Finished Import\n";

