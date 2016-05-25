#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Spec;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin";
use configdata;
use parsers;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Retrieve reads supporting a fusion.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -i, --id        Cluster ID (mutually exclusive with --list)\n";
push @usage, "  -l, --list      Filename of Cluster ID List\n";
push @usage, "      --fastq1    Filename of End 1 Fastq\n";
push @usage, "      --fastq2    Filename of End 2 Fastq\n";

my $help;
my $config_filename;
my $dataset_directory;
my $output_directory;
my $query_cluster_id;
my $query_filename;
my $fastq1_output_filename;
my $fastq2_output_filename;

GetOptions
(
	'help'        => \$help,
	'output=s'    => \$output_directory,
	'id=s'        => \$query_cluster_id,
	'list=s'      => \$query_filename,
	'fastq1=s'    => \$fastq1_output_filename,
	'fastq2=s'    => \$fastq2_output_filename,
);

not defined $help or die @usage;

defined $output_directory or die @usage;
defined $query_cluster_id or defined $query_filename or die @usage;
not (defined $query_cluster_id and defined $query_filename) or die @usage;
defined $fastq1_output_filename or die @usage;
defined $fastq2_output_filename or die @usage;

my %cluster_ids;

if (defined $query_cluster_id)
{
	$cluster_ids{$query_cluster_id} = 1;
}
elsif (defined $query_filename)
{
	open LST, $query_filename or die "Error: Unable to find $query_filename: $!\n";
	while (<LST>)
	{
		chomp;
		$cluster_ids{$_} = 1;
	}
	close LST;	
}
else
{
	die;
}

my $reads_prefix = $output_directory."/reads";
my $clusters_sc = $output_directory."/clusters.sc";
my $split_pred_align = $output_directory."/splitreads.predalign";

my $reads_fastq_index = $reads_prefix.".fqi";
my $reads_end_1_fastq = $reads_prefix.".1.fastq";
my $reads_end_2_fastq = $reads_prefix.".2.fastq";

my %fragment_ids;
my %fragment_names;
read_clusters($clusters_sc, \%cluster_ids, \%fragment_ids, \%fragment_names);
read_splitread_ids($split_pred_align, \%cluster_ids, \%fragment_ids, \%fragment_names);

output_reads_to_fastq(\%fragment_ids, \%fragment_names, $reads_fastq_index, 
	$reads_end_1_fastq, $reads_end_2_fastq, 
	$fastq1_output_filename, $fastq2_output_filename);

sub read_clusters
{
	my $clusters_filename = shift;
	my $cluster_ids = shift;
	my $fragment_ids = shift;
	my $fragment_names = shift;
	
	open CLU, $clusters_filename or die "Error: Unable to find $clusters_filename: $!\n";
	while (<CLU>)
	{
		chomp;
		my @fields = split /\t/;

		my $cluster_id = $fields[0];
		my $fragment_id = $fields[2];
		
		next unless defined $cluster_ids->{$cluster_id};

		$fragment_ids->{$fragment_id} = 1;
		$fragment_names->{$fragment_id} = $cluster_id."_".$fragment_id;
	}
	close CLU;
}

sub read_splitread_ids
{
	my $split_readsids_filename = shift;
	my $cluster_ids = shift;
	my $fragment_ids = shift;
	my $fragment_names = shift;

	open RID, $split_readsids_filename or die "Error: Unable to find $split_readsids_filename: $!\n";
	while (<RID>)
	{
		chomp;
		my @fields = split /\t/;

		my $cluster_id = $fields[0];
		my $fragment_id = $fields[1];
		
		next unless defined $cluster_ids->{$cluster_id};

		$fragment_ids->{$fragment_id} = 1;
		$fragment_names->{$fragment_id} = $cluster_id."_".$fragment_id;
	}
	close RID;
}

sub output_reads_to_fastq
{
	my $fragments_ref = shift;
	my $fragment_names_ref = shift;
	my $fastq_index_filename = shift;
	my $fastq1_filename = shift;
	my $fastq2_filename = shift;
	my $fastq1_output_filename = shift;
	my $fastq2_output_filename = shift;
	
	open FQI, $fastq_index_filename or die "Error: Unable to open $fastq_index_filename\n";
	binmode(FQI);
	
	foreach my $read_end ("1", "2")
	{
		my $fastq_filename = ($read_end eq "1") ? $fastq1_filename : $fastq2_filename;
		open FQ, $fastq_filename or die "Error: Unable to open $fastq_filename\n";
		
		my $fastq_output_filename = ($read_end eq "1") ? $fastq1_output_filename : $fastq2_output_filename;
		open FQO, ">".$fastq_output_filename or die "Error: Unable to open $fastq_output_filename\n";
		
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
			
			$read_id =~ /\@$fragment_index\/$read_end/ or die "Error: Unable to retrieve correct index for $fragment_index/$read_end\n";

			my $read_name = $fragment_names_ref->{$fragment_index}."/".$read_end;
			
			print FQO "@".$read_name."\n";
			print FQO $sequence;
			print FQO $comment;
			print FQO $quality;
		}
		
		close FQ;
		close FQO;
	}
	
	close FQI;
}


