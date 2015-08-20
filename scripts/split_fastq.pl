#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 fastq_prefix reads_per_file split_prefix > split_catalog\n";

my $fastq_prefix = shift;
my $reads_per_file = shift;
my $split_prefix = shift;

die $usage if not defined $split_prefix;

my @read_ends;
foreach my $read_end ("1","2")
{
	push @read_ends, $read_end;
}

my @fastq_files;
my @split_fastq_files;
foreach my $read_end_index (0..1)
{
	my $fastq_filename = $fastq_prefix.".".$read_ends[$read_end_index].".fastq";

	my $fastq_file;
	open $fastq_file, $fastq_filename or die "Error: Unable to open $fastq_filename: $!\n";
	push @fastq_files, $fastq_file;

	push @split_fastq_files, undef;
}

my @split_fastq_prefixes;
my %first_fragment_index;
my %last_fragment_index;
my $reads_this_file = $reads_per_file;

while (1)
{
	my @fastq_entries;
	foreach my $fastq_file (@fastq_files)
	{
		push @fastq_entries, read_fastq_entry($fastq_file);
	}
	
	last if not check_fastq_entries(@fastq_entries);
	
	my $fragment_index = get_entries_fragment_index(@fastq_entries);
	
	if ($reads_this_file == $reads_per_file)
	{
		my $file_num_str = scalar @split_fastq_prefixes;
		$file_num_str = "0".$file_num_str while length($file_num_str) < 3;
		
		my $split_fastq_prefix = $split_prefix.".split.".$file_num_str;

		foreach my $read_end_index (0..1)
		{
			my $split_fastq_filename = $split_fastq_prefix.".".$read_ends[$read_end_index].".fastq";
			
			close $split_fastq_files[$read_end_index] if defined $split_fastq_files[$read_end_index];
		
			open $split_fastq_files[$read_end_index], ">".$split_fastq_filename or die "Error: Unable to open $split_fastq_filename: $!\n";
			
			$first_fragment_index{$split_fastq_prefix} = $fragment_index;
		}
		
		$reads_this_file = 0;
		
		push @split_fastq_prefixes, $split_fastq_prefix;
	}
	
	$last_fragment_index{$split_fastq_prefixes[$#split_fastq_prefixes]} = $fragment_index;
	
	foreach my $read_end_index (0..1)
	{
		write_fastq_entry($split_fastq_files[$read_end_index], $fastq_entries[$read_end_index]);
	}

	$reads_this_file++;
}

foreach my $read_end_index (0..1)
{
	close $fastq_files[$read_end_index];
	close $split_fastq_files[$read_end_index] if defined $split_fastq_files[$read_end_index];
}

foreach my $split_fastq_prefix (@split_fastq_prefixes)
{
	print $split_fastq_prefix."\t".$first_fragment_index{$split_fastq_prefix}."\t".$last_fragment_index{$split_fastq_prefix}."\n";
}

sub read_fastq_entry
{
	my $file = $_[0];

	my $readid = <$file>;
	my $sequence = <$file>;
	my $comment = <$file>;
	my $quality = <$file>;
	
	return [$readid,$sequence,$comment,$quality];
}

sub write_fastq_entry
{
	my $file = $_[0];
	my $entry = $_[1];

	print $file $entry->[0];
	print $file $entry->[1];
	print $file $entry->[2];
	print $file $entry->[3];
}

sub check_fastq_entries
{
	return (defined $_[0]->[0] and defined $_[0]->[1] and defined $_[0]->[2] and defined $_[0]->[3] and
	        defined $_[1]->[0] and defined $_[1]->[1] and defined $_[1]->[2] and defined $_[1]->[3]);
}

sub interpret_read_id
{
	my $entry = $_[0];
	
	$entry->[0] =~ /(\d+)\/(\d)/ or die "Error: Unable to interpret read id $_[0]->[0]\n";
	
	return [$1,$2];
}

sub get_entries_fragment_index
{
	my $entry1 = $_[0];
	my $entry2 = $_[1];
	
	my $read_id1 = interpret_read_id($entry1);
	my $read_id2 = interpret_read_id($entry2);
	
	($read_id1->[0] == $read_id2->[0]) or die "Error: Inconsistent fragment index for $entry1->[0] and $entry2->[0]\n";
	($read_id1->[1] == "1") and ($read_id2->[1] == "2") or die "Error: Inconsistent read end for $entry1->[0] and $entry2->[0]\n";
	
	return $read_id1->[0];
}

