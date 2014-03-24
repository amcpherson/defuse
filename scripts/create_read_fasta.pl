#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

my $clusters_sc_filename = shift;
my $reads_prefix = shift;

defined $reads_prefix or die "Usage: $0 clusters_sc_filename reads_prefix\n";

# Read clusters
my %fragments;
read_cluster_fragments($clusters_sc_filename, \%fragments);

# Create a read sequence fasta
my $reads_fastq_index = $reads_prefix.".fqi";
my $reads_end_1_fastq = $reads_prefix.".1.fastq";
my $reads_end_2_fastq = $reads_prefix.".2.fastq";
output_reads_to_fasta(\%fragments, $reads_fastq_index, $reads_end_1_fastq, $reads_end_2_fastq);

sub read_cluster_fragments
{
	my $clusters_filename = shift;
	my $fragments_hash_ref = shift;
	
	open CLU, $clusters_filename or die "Error: Unable to open $clusters_filename: $!\n";
	while (<CLU>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $cluster_end = $fields[1];
		my $fragment_id = $fields[2];
		my $read_end = $fields[3];
		my $ref_name = $fields[4];
		my $strand = $fields[5];
		my $start = $fields[6];
		my $end = $fields[7];
		
		$fragments_hash_ref->{$fragment_id} = 1;
	}
	close CLU;
}

sub output_reads_to_fasta
{
	my $fragments_ref = shift;
	my $fastq_index_filename = shift;
	my $end_1_fastq_filename = shift;
	my $end_2_fastq_filename = shift;
	
	open FQI, $fastq_index_filename or die "Error: Unable to open $fastq_index_filename\n";
	binmode(FQI);
	
	foreach my $read_end ("1","2")
	{
		my $fastq_filename = ($read_end eq "1") ? $end_1_fastq_filename : $end_2_fastq_filename;
		open FQ, $fastq_filename or die "Error: Unable to open $fastq_filename\n";
		
		my $read_end_index = ($read_end eq "1") ? 0 : 1;
		foreach my $fragment_index (sort {$a <=> $b} keys %{$fragments_ref})
		{
			# Seek to position in the index
			my $fragment_index_file_pos = $fragment_index * 2 * 8 + $read_end_index * 8;
			seek FQI, $fragment_index_file_pos, 0;
			
			# Read from index
			my $fragment_file_pos_bytes;
			read FQI, $fragment_file_pos_bytes, 8;
			my $fragment_file_pos = unpack "q", $fragment_file_pos_bytes;
			
			# Seek to position in fastq
			seek FQ, $fragment_file_pos, 0;
			
			# Read fastq entry
			my $read_id = <FQ>;
			my $sequence = <FQ>;
			my $comment = <FQ>;
			my $quality = <FQ>;
			
			chomp($read_id);
			chomp($sequence);
			
			$read_id =~ s/\@// or die "Error: Unable to retrieve correct index for $fragment_index/$read_end\n";
			$read_id = $fragment_index."/".$read_end or die "Error: Unable to retrieve correct index for $fragment_index/$read_end\n";
			
			print ">".$read_id."\n";
			print $sequence."\n";
		}
		
		close FQ;
	}
	
	close FQI;
}

