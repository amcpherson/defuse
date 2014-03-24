#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 fasta_filename sequences_per_file split_prefix > split_catalog\n";

my $fasta_filename = shift;
my $sequences_per_file = shift;
my $split_prefix = shift;

die $usage if not defined $split_prefix;

my $file_num = 1;
my $split_fasta_file;
my $split_fasta_filename;
my $sequences_this_file = $sequences_per_file;

sub write_seq
{
	my $seq_id = shift;
	my $seq = shift;
	
	return if not defined $seq_id;
	
	if ($sequences_this_file == $sequences_per_file)
	{
		my $file_num_str = $file_num;
		$file_num_str = "0".$file_num_str while length($file_num_str) < 3;

		my $split_fasta_filename = $split_prefix.".split.".$file_num_str.".fa";

		close $split_fasta_file if defined $split_fasta_file;
		open $split_fasta_file, ">".$split_fasta_filename or die "Error: Unable to open $split_fasta_filename: $!\n";
		
		print $split_fasta_filename."\n";
		
		$file_num++;
		
		$sequences_this_file = 0;
	}
	
	print $split_fasta_file ">".$seq_id."\n".$seq."\n";
	$sequences_this_file++;
}

my $current_seq_id;
my $current_seq;

open FA, $fasta_filename or die "Error: Unable to open $fasta_filename: $!\n";
while (<FA>)
{
	chomp;
	
	if (/^>(.*)$/)
	{
		my $next_seq_id = $1;
		
		write_seq($current_seq_id, $current_seq);
		
		$current_seq_id = $next_seq_id;
		$current_seq = "";
	}
	else
	{
		$current_seq = $current_seq.$_;
	}
}

write_seq($current_seq_id, $current_seq);

close FA;
close $split_fasta_file if defined $split_fasta_file;

	
