#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];

use lib dirname($0);
use configdata;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Coallate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -l, --list      List of Cluster IDs\n";

my $help;
my $config_filename;
my $output_directory;
my $cluster_list;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
	'list=s'      => \$cluster_list,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $output_directory or die @usage;
defined $cluster_list or die @usage;

my $config = configdata->new();
$config->read($config_filename);

my $scripts_directory		= $config->get_value("scripts_directory");
my $denovo_assembly			= $config->get_value("denovo_assembly");

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
my $clusters_sam = $output_directory."/clusters.sc";
read_sam_clusters($clusters_sam, \%clusters);

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
my $mapping_stats_filename = $output_directory."/mapping.stats";
read_annotations($annotations_filename, \%annotations);
read_annotations($mapping_stats_filename, \%annotations);

# Read in split read and denovo break predictions
my %splitr_break;
my %denovo_break;
my $splitr_break_filename = $output_directory."/splitr.break";
my $denovo_break_filename = $output_directory."/denovo.break";
read_breaks($splitr_break_filename, \%splitr_break);
read_breaks($denovo_break_filename, \%denovo_break);

# Read in split read and denovo seq predictions
my %splitr_seq;
my %denovo_seq;
my $splitr_seq_filename = $output_directory."/splitr.seq";
my $denovo_seq_filename = $output_directory."/denovo.seq";
read_splitr_seq($splitr_seq_filename, \%splitr_seq);
read_denovo_seq($denovo_seq_filename, \%denovo_seq);

# Read in spanning pvalue file for split and denovo predictions
my %splitr_span_pval;
my %denovo_span_pval;
my $splitr_span_pval_filename = $output_directory."/splitr.span.pval";
my $denovo_span_pval_filename = $output_directory."/denovo.span.pval";
read_span_pval($splitr_span_pval_filename, \%splitr_span_pval);
read_span_pval($denovo_span_pval_filename, \%denovo_span_pval);

# Read in split pvalue file for split predictions
my %splitr_split_pval;
my $splitr_split_pval_filename = $output_directory."/splitr.split.pval";
read_split_pval($splitr_split_pval_filename, \%splitr_split_pval);

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

if (lc($denovo_assembly) eq "yes")
{
	print "denovo_sequence\t";
	print "denovo_min_count\t";
	print "denovo_span_pvalue\t";
}

foreach my $anno_type (@anno_types)
{
	print $anno_type."\t";
}
print "\n";

foreach my $cluster_id (keys %cluster_ids)
{
	print $cluster_id."\t";

	print $splitr_seq{$cluster_id}{sequence}."\t";
	print $splitr_seq{$cluster_id}{split_count}."\t";
	print $splitr_span_pval{$cluster_id}{pvalue}."\t";
	print $splitr_split_pval{$cluster_id}{pos_pvalue}."\t";
	print $splitr_split_pval{$cluster_id}{min_pvalue}."\t";

	if (lc($denovo_assembly) eq "yes")
	{
		print $denovo_seq{$cluster_id}{sequence}."\t";
		print $denovo_seq{$cluster_id}{min_count}."\t";
		print $denovo_span_pval{$cluster_id}{pvalue}."\t";
	}

	foreach my $anno_type (@anno_types)
	{
		my $anno_value = $annotations{$cluster_id}{$anno_type};
		$anno_value = "" if not defined $anno_value;
		print $anno_value."\t";
	}
	
	print "\n";
}

sub read_sam_clusters
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

sub read_splitr_seq
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

sub read_denovo_seq
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
		$seqs_hash_ref->{$cluster_id}{min_count} = $fields[3];
	}
	close SEQ;
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

