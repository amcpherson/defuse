#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 fastq_filename reads_per_file split_filename\n\tA % in split_filename will be replaced by the file number\n";

my $fastq_filename = shift;
my $reads_per_file = shift;
my $split_filename = shift;

die $usage if not defined $split_filename;

my $file_num = 1;
my $file;
my $reads_this_file = $reads_per_file;

open FQ, $fastq_filename or die "Error: Unable to open $fastq_filename: $!\n";
while (1)
{
	my $readid = <FQ>;
	my $sequence = <FQ>;
	my $comment = <FQ>;
	my $quality = <FQ>;

	last if not defined $quality;

	chomp($readid);
	chomp($sequence);
	chomp($comment);
	chomp($quality);

	if ($reads_this_file == $reads_per_file)
	{
		my $file_num_str = $file_num;
		$file_num_str = "0".$file_num_str while length($file_num_str) < 3;

		my $filename = $split_filename;
		$filename =~ s/%/$file_num_str/;

		print $filename."\n";

		close $file if defined $file;

		open $file, ">".$filename or die "Error: Unable to open $filename: $!\n";

		$reads_this_file = 0;
		$file_num++;
	}

	$reads_this_file++;

	print $file "$readid\n$sequence\n$comment\n$quality\n";
}
close FQ;

close $file if defined $file;

