#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use FileCache;
no strict 'refs';

my $usage = "Usage: $0 fasta_filename fasta_out_filename1 fasta_out_filename2 ...\n";

my @filenames = @ARGV;

my $input_fasta_filename = shift @filenames;
my @output_fasta_filenames = @filenames;

die $usage if not defined $input_fasta_filename or scalar @output_fasta_filenames == 0;

my @output_fasta_files;
foreach my $output_fasta_filename (@output_fasta_filenames)
{
	my $output_fasta_file = $output_fasta_filename;
	cacheout $output_fasta_file;
	push @output_fasta_files, $output_fasta_file;
}

my $out_file_index = 0;

sub write_seq
{
	my $seq_id = shift;
	my $seq = shift;
	
	return if not defined $seq_id;
	
	my $output_fasta_file = $output_fasta_files[$out_file_index];
	
	print $output_fasta_file ">".$seq_id."\n".$seq."\n";

	$out_file_index = ($out_file_index + 1) % scalar @output_fasta_filenames;
}

my $current_seq_id;
my $current_seq;

open FA, $input_fasta_filename or die "Error: Unable to open $input_fasta_filename: $!\n";
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

	
