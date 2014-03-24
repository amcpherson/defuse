#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;

use lib dirname($0);
use configdata;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -o, --output    Output Directory\n";

my $help;
my $config_filename;
my $output_directory;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $output_directory or die @usage;

my $config = configdata->new();
$config->read($config_filename);

my $genome_fasta                      = $config->get_value("genome_fasta");
my $est_fasta                         = $config->get_value("est_fasta");
my $cdna_fasta                        = $config->get_value("cdna_fasta");
my $rrna_fasta                        = $config->get_value("rrna_fasta");
my $est_alignments                    = $config->get_value("est_alignments");
my $span_coverage_threshold           = $config->get_value("span_coverage_threshold");
my $span_count_threshold              = $config->get_value("span_count_threshold");
my $max_concordant_ratio              = $config->get_value("max_concordant_ratio");
my $split_count_threshold             = $config->get_value("split_count_threshold");
my $span_length_pvalue_threshold      = $config->get_value("span_length_pvalue_threshold");
my $split_pos_pvalue_threshold        = $config->get_value("split_pos_pvalue_threshold");
my $split_min_pvalue_threshold        = $config->get_value("split_min_pvalue_threshold");
my $percent_identity_threshold        = $config->get_value("percent_identity_threshold");
my $genome_breakseq_threshold         = $config->get_value("genome_breakseq_threshold");
my $cdna_breakseq_threshold           = $config->get_value("cdna_breakseq_threshold");
my $est_breakseq_threshold            = $config->get_value("est_breakseq_threshold");
my $estislands_breakseq_threshold     = $config->get_value("estislands_breakseq_threshold");
my $scripts_directory                 = $config->get_value("scripts_directory");

my %span_coverage_failed;
my %span_count_failed;
my %span_align_failed;
my %split_count_failed;
my %split_pos_pvalue_failed;
my %split_min_pvalue_failed;
my %span_length_pvalue_failed;
my %alignment_failed;
my %estislands_failed;

my $log_filename = $output_directory."/filtering.log";
open LOG, ">".$log_filename or die "Error: Unable to open $log_filename\n";

# Read in cluster sam file
my %clusters;
my $clusters_sam = $output_directory."/clusters.sc.sam";
read_sam_clusters($clusters_sam, \%clusters);

# Find cluster spanning read count
my %fusion_span_count;
foreach my $cluster_id (keys %clusters)
{
	die "Error: fusion $cluster_id does not refer to 2 reference sequences\n" if scalar keys %{$clusters{$cluster_id}} != 2;
	foreach my $ref_name (keys %{$clusters{$cluster_id}})
	{
		$fusion_span_count{$cluster_id} = scalar keys %{$clusters{$cluster_id}{$ref_name}};
	}
}

# Read in annotations file
my %annotations;
my $annotations_filename = $output_directory."/annotations.txt";
read_annotations($annotations_filename, \%annotations);

# Read in split read seq predictions
my %splitr_seq;
my $splitr_seq = $output_directory."/splitr.seq";
read_split_seq($splitr_seq, \%splitr_seq);

# Read in spanning pvalue file for split and denovo predictions
my %splitr_span_pval;
my %denovo_span_pval;
my $splitr_span_pval_filename = $output_directory."/splitr.span.pval";
my $denovo_span_pval_filename = $output_directory."/denovo.span.pval";
read_span_pval($splitr_span_pval_filename, \%splitr_span_pval);
read_span_pval($denovo_span_pval_filename, \%denovo_span_pval);

# Read in split pvalues file for split predictions
my %splitr_split_pval;
my $splitr_split_pval_filename = $output_directory."/splitr.split.pval";
read_split_pval($splitr_split_pval_filename, \%splitr_split_pval);

my %false_pos;
foreach my $cluster_id (keys %annotations)
{
	if ($annotations{$cluster_id}{break_predict} eq "splitr")
	{
		print LOG "$cluster_id failed split_pos_pvalue_threshold\n" and $false_pos{$cluster_id} = 1 and $split_pos_pvalue_failed{$cluster_id} = 1 if $splitr_split_pval{$cluster_id}{pos_pvalue} < $split_pos_pvalue_threshold;
		print LOG "$cluster_id failed split_min_pvalue_threshold\n" and $false_pos{$cluster_id} = 1 and $split_min_pvalue_failed{$cluster_id} = 1 if $splitr_split_pval{$cluster_id}{min_pvalue} < $split_min_pvalue_threshold;
		print LOG "$cluster_id failed span_length_pvalue_threshold\n" and $false_pos{$cluster_id} = 1 and $span_length_pvalue_failed{$cluster_id} = 1 if $splitr_span_pval{$cluster_id}{pvalue} < $span_length_pvalue_threshold;
	}
	elsif ($annotations{$cluster_id}{break_predict} eq "denovo")
	{
		print LOG "$cluster_id failed span_length_pvalue_threshold\n" and $false_pos{$cluster_id} = 1 and $span_length_pvalue_failed{$cluster_id} = 1 if $denovo_span_pval{$cluster_id}{pvalue} < $span_length_pvalue_threshold;
	}
	else
	{
		die "Error: Unrecognized break_predict value ".$annotations{$cluster_id}{break_predict}."\n";
	}
	
	my $gene1 = $annotations{$cluster_id}{gene1};
	my $gene2 = $annotations{$cluster_id}{gene2};	

	print LOG "$cluster_id failed span_coverage_threshold $gene1\n" and $false_pos{$cluster_id} = 1 and $span_coverage_failed{$cluster_id} = 1 if $annotations{$cluster_id}{span_coverage1} < $span_coverage_threshold;
	print LOG "$cluster_id failed span_coverage_threshold $gene2\n" and $false_pos{$cluster_id} = 1 and $span_coverage_failed{$cluster_id} = 1 if $annotations{$cluster_id}{span_coverage2} < $span_coverage_threshold;
	
	print LOG "$cluster_id failed span_count_threshold\n" and $false_pos{$cluster_id} = 1 and $span_count_failed{$cluster_id} = 1 if $fusion_span_count{$cluster_id} < $span_count_threshold;
	print LOG "$cluster_id failed split_count_threshold\n" and $false_pos{$cluster_id} = 1 and $split_count_failed{$cluster_id} = 1 if $splitr_seq{$cluster_id}{split_count} < $split_count_threshold;

	print LOG "$cluster_id failed concordant alignment filter\n" and $false_pos{$cluster_id} = 1 and $span_align_failed{$cluster_id} = 1 if $annotations{$cluster_id}{concordant_ratio} > $max_concordant_ratio;

	print LOG "$cluster_id failed genome alignment\n" and $false_pos{$cluster_id} = 1 and $alignment_failed{$cluster_id} = 1 if $annotations{$cluster_id}{genome_breakseqs_percident} > $genome_breakseq_threshold;
	print LOG "$cluster_id failed cdna alignment\n" and $false_pos{$cluster_id} = 1 and $alignment_failed{$cluster_id} = 1 if $annotations{$cluster_id}{cdna_breakseqs_percident} > $cdna_breakseq_threshold;
	print LOG "$cluster_id failed est alignment\n" and $false_pos{$cluster_id} = 1 and $alignment_failed{$cluster_id} = 1 if $annotations{$cluster_id}{est_breakseqs_percident} > $est_breakseq_threshold;
	print LOG "$cluster_id failed estislands\n" and $false_pos{$cluster_id} = 1 and $estislands_failed{$cluster_id} = 1 if $annotations{$cluster_id}{breakseqs_estislands_percident} > $estislands_breakseq_threshold;
}

foreach my $cluster_id (keys %annotations)
{
	next if defined $false_pos{$cluster_id};
	print $cluster_id."\n";
}

print LOG "Failure Statistics\n";
print LOG "Total\t".scalar(keys(%annotations))."\n";
print LOG "Span Coverage\t".scalar(keys(%span_coverage_failed))."\n";
print LOG "Span Count\t".scalar(keys(%span_count_failed))."\n";
print LOG "Span Alignment\t".scalar(keys(%span_align_failed))."\n";
print LOG "Split Count\t".scalar(keys(%split_count_failed))."\n";
print LOG "Split Position P-Value\t".scalar(keys(%split_pos_pvalue_failed))."\n";
print LOG "Split Minimum P-Value\t".scalar(keys(%split_min_pvalue_failed))."\n";
print LOG "Span Length P-Value\t".scalar(keys(%span_length_pvalue_failed))."\n";
print LOG "Boundary Alignment\t".scalar(keys(%alignment_failed))."\n";
print LOG "EST Islands\t".scalar(keys(%estislands_failed))."\n";

print LOG "Attrition Statistics\n";
my %cumulative_failed;
my $num_failed;
$num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Total\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%span_coverage_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Span Coverage\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%span_align_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Span Align\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%span_count_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Span Count\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%split_count_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Split Count\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%split_pos_pvalue_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Split Position P-Value\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%split_min_pvalue_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Split Minimum P-Value\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%span_length_pvalue_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Span Length P-Value\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%alignment_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "Boundary Alignment\t".$num_failed."\n";
%cumulative_failed = (%cumulative_failed,%estislands_failed); $num_failed = scalar(keys(%annotations)) - scalar(keys(%cumulative_failed)); print LOG "EST Islands\t".$num_failed."\n";

close LOG;

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
		my $read_id = $fields[1];
		my $flag = $fields[2];
		my $ref_name = $fields[3];
		my $start = $fields[4];
		my $end = $start + length($fields[10]) - 1;
		my $sequence = $fields[10];
	
		$read_id =~ /(.*)\/([12])/;
		my $fragment_id = $1;
		my $read_end = $2;
	
		my $strand;
		if ($flag & hex('0x0010'))
		{
			$strand = "-";
		}
		else
		{
			$strand = "+";
		}
		
		if ($strand eq "-")
		{
			$sequence = reverse($sequence);
			$sequence =~ tr/ACGTacgt/TGCAtgca/;
		}
	
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{read_id} = $read_id;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{read_end} = $read_end;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{strand} = $strand;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{sequence} = $sequence;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{start} = $start;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{end} = $end;
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
		
		$span_pval_hash_ref->{$cluster_id}{pvalue} = $fields[2];
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

sub read_split_seq
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

