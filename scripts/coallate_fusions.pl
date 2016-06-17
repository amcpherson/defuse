#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];

use FindBin;
use lib "$FindBin::RealBin";
use configdata;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Coallate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -l, --list      List of Cluster IDs\n";

my $help;
my $config_filename;
my $dataset_directory;
my $output_directory;
my $cluster_list;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'output=s'    => \$output_directory,
	'list=s'      => \$cluster_list,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $dataset_directory or die @usage;
defined $output_directory or die @usage;
defined $cluster_list or die @usage;

my $source_directory = abs_path("$FindBin::RealBin/../");

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

my $scripts_directory		= $config->get_value("scripts_directory");

my %span_count_failed;
my %span_align_failed;
my %split_count_failed;
my %split_pos_pvalue_failed;
my %split_min_pvalue_failed;
my %span_length_pvalue_failed;
my %alignment_failed;
my %estislands_failed;

# Read in list of clusters to coallate
my %cluster_ids;
open CID, $cluster_list or die "Error: Unable to open $cluster_list: $!\n";
while (<CID>)
{
	chomp;
	my ($cluster_id) = split /\t/;
	$cluster_ids{$cluster_id} = 1;
}
close CID;

# Read in cluster sam file
my %clusters;
my $clusters_sc = $output_directory."/clusters.sc";
read_clusters($clusters_sc, \%clusters);

# Find cluster spanning read count and strands
my %fusion_span_count;
my %fusion_strand;
foreach my $cluster_id (keys %clusters)
{
	foreach my $cluster_end (keys %{$clusters{$cluster_id}})
	{
		$fusion_span_count{$cluster_id} = scalar keys %{$clusters{$cluster_id}{$cluster_end}};
		
		my $fragment_id = (keys %{$clusters{$cluster_id}{$cluster_end}})[0];
		$fusion_strand{$cluster_id}{$cluster_end} = $clusters{$cluster_id}{$cluster_end}{$fragment_id}{strand};
	}
}

# Read in annotations file
my %annotations;
my $annotations_filename = $output_directory."/annotations";
read_annotations($annotations_filename, \%annotations);

# Read in split read break predictions
my %splitreads_break;
my $splitreads_break_filename = $output_directory."/splitreads.break";
read_breaks($splitreads_break_filename, \%splitreads_break);

# Read in split read seq predictions
my %splitreads_seq;
my $splitreads_seq_filename = $output_directory."/splitreads.seq";
read_seq($splitreads_seq_filename, \%splitreads_seq);

# Read in spanning pvalue file for split and denovo predictions
my %splitreads_span_pval;
my $splitreads_span_pval_filename = $output_directory."/splitreads.span.pval";
read_span_pval($splitreads_span_pval_filename, \%splitreads_span_pval);

# Read in split pvalue file for split predictions
my %splitreads_split_pval;
my $splitreads_split_pval_filename = $output_directory."/splitreads.split.pval";
read_split_pval($splitreads_split_pval_filename, \%splitreads_split_pval);

# Collect a list of all annotation types
my %anno_types_hash;
foreach my $cluster_id (keys %annotations)
{
	foreach my $anno_type (keys %{$annotations{$cluster_id}})
	{
		$anno_types_hash{$anno_type} = 1;
	}
}
my @anno_types = sort {$a cmp $b} (keys %anno_types_hash);

print "cluster_id\t";

print "splitr_sequence\t";
print "splitr_count\t";
print "splitr_span_pvalue\t";
print "splitr_pos_pvalue\t";
print "splitr_min_pvalue\t";

foreach my $anno_type (@anno_types)
{
	print $anno_type."\t";
}
print "\n";

foreach my $cluster_id (keys %cluster_ids)
{
	print $cluster_id."\t";

	print $splitreads_seq{$cluster_id}{sequence}."\t";
	print $splitreads_seq{$cluster_id}{split_count}."\t";
	print $splitreads_span_pval{$cluster_id}{pvalue}."\t";
	print $splitreads_split_pval{$cluster_id}{pos_pvalue}."\t";
	print $splitreads_split_pval{$cluster_id}{min_pvalue}."\t";
	
	foreach my $anno_type (@anno_types)
	{
		my $anno_value = $annotations{$cluster_id}{$anno_type};
		$anno_value = "" if not defined $anno_value;
		print $anno_value."\t";
	}
	
	print "\n";
}

sub read_clusters
{
	my $clusters_filename = shift;
	my $clusters_hash_ref = shift;
	
	open CLU, $clusters_filename or die "Error: Unable to find $clusters_filename: $!\n";
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
		
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{read_id} = $fragment_id."/".$read_end;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{read_end} = $read_end;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{strand} = $strand;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{start} = $start;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{end} = $end;
	}
	close CLU;
}

sub read_span_pval
{
	my $span_pval_filename = shift;
	my $span_pval_hash_ref = shift;
	
	open SPP, $span_pval_filename or die "Error: Unable to find $span_pval_filename: $!\n";
	while (<SPP>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		
		$span_pval_hash_ref->{$cluster_id}{pvalue} = $fields[1];
	}
	close SPP;
}

sub read_split_pval
{
	my $split_pval_filename = shift;
	my $split_pval_hash_ref = shift;
	
	open SPP, $split_pval_filename or die "Error: Unable to find $split_pval_filename: $!\n";
	while (<SPP>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		
		$split_pval_hash_ref->{$cluster_id}{pos_pvalue} = $fields[1];
		$split_pval_hash_ref->{$cluster_id}{min_pvalue} = $fields[2];
	}
	close SPP;
}

sub read_seq
{
	my $seqs_filename = shift;
	my $seqs_hash_ref = shift;
	
	open SEQ, $seqs_filename or die "Error: Unable to find $seqs_filename: $!\n";
	while (<SEQ>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		
		$seqs_hash_ref->{$cluster_id}{sequence} = $fields[1];
		$seqs_hash_ref->{$cluster_id}{inter_length} = $fields[2];
		$seqs_hash_ref->{$cluster_id}{split_count} = $fields[3];
		$seqs_hash_ref->{$cluster_id}{split_pos_average} = $fields[4];
		$seqs_hash_ref->{$cluster_id}{split_min_average} = $fields[5];
	}
	close SEQ;
}

sub read_breaks
{
	my $breaks_filename = shift;
	my $breaks_hash_ref = shift;
	
	open BR, $breaks_filename or die "Error: Unable to find $breaks_filename: $!\n";
	while (<BR>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $reference = $fields[1];
		my $strand = $fields[2];
		my $breakpos = $fields[3];
		
		push @{$breaks_hash_ref->{$cluster_id}{breakpos}}, [$reference,$strand,$breakpos];
	}
	close BR;
}

sub read_annotations
{
	my $annotations_filename = shift;
	my $annotations_hash_ref = shift;
	
	open ANO, $annotations_filename or die "Error: Unable to find $annotations_filename: $!\n";
	while (<ANO>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $anno_type = $fields[1];
		my $anno_value = $fields[2];
		
		$annotations_hash_ref->{$cluster_id}{$anno_type} = $anno_value;
	}
	close ANO;
}

