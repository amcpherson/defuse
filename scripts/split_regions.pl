#!/usr/bin/perl

use strict;
use warnings;

my $usage = "Usage: $0 regions_filename regions_per_file split_filename\n\tA % in split_filename will be replaced by the file number\n";

my $regions_filename = shift;
my $regions_per_file = shift;
my $split_filename = shift;

die $usage if not defined $split_filename;

my $file_num = 1;
my $file;
my $regions_this_file = $regions_per_file;

open REG, $regions_filename or die "Error: Unable to open $regions_filename: $!\n";
while (1)
{
	my $regions1 = <REG>;
	my $regions2 = <REG>;

	last if not defined $regions2;

	if ($regions_this_file == $regions_per_file)
	{
		my $file_num_str = $file_num;
		$file_num_str = "0".$file_num_str while length($file_num_str) < 3;

		my $filename = $split_filename;
		$filename =~ s/%/$file_num_str/;

		print $filename."\n";

		close $file if defined $file;

		open $file, ">".$filename or die "Error: Unable to open $filename: $!\n";

		$regions_this_file = 0;
		$file_num++;
	}

	$regions_this_file++;

	print $file $regions1;
	print $file $regions2;
}
close REG;

close $file if defined $file;


