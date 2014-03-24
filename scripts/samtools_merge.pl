#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $byread;
GetOptions ("n" => \$byread);

my $samtools_bin = shift;
my $output_filename = shift;

my @input_filenames;
while (my $input_filename = shift)
{
	push @input_filenames, $input_filename;
}

scalar @input_filenames > 0 or die "Usage: $0 samtools_bin output_filename input_filename1 input_filename2 ...\n";

my @nonempty_input_filenames;
foreach my $input_filename (@input_filenames)
{
	my $nonempty = 0;
	
	open SAM, "$samtools_bin view $input_filename |" or die "Error: Unable to samtools view $input_filename\n";
	while (<SAM>)
	{
		next if /^\@/;
		
		if (length($_) > 0)
		{
			$nonempty = 1;
			last;
		}
	}
	close SAM;
	
	if ($nonempty)
	{
		push @nonempty_input_filenames, $input_filename;
	}
}

if (scalar @nonempty_input_filenames <= 1)
{
	rename $input_filenames[0], $output_filename;
}
else
{
	my $command = join " ", ("$samtools_bin", "merge", (defined $byread) ? "-n" : "", $output_filename, @nonempty_input_filenames);
	my $result = system $command;
	
	$result == 0 or die "Error: Merge failed\n";
}

