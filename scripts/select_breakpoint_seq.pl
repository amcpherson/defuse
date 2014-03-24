#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -o, --output    Output Directory\n";

my $help;
my $output_directory;

GetOptions
(
	'help'        => \$help,
	'output=s'    => \$output_directory,
);

not defined $help or usage() and exit;

defined $output_directory or die @usage;

my $splitr_break_filename = $output_directory."/splitr.break";
my $denovo_break_filename = $output_directory."/denovo.break";
my $splitr_seq_filename = $output_directory."/splitr.seq";
my $denovo_seq_filename = $output_directory."/denovo.seq";
my $splitr_span_pval_filename = $output_directory."/splitr.span.pval";
my $denovo_span_pval_filename = $output_directory."/denovo.span.pval";

my $break_predict_filename = $output_directory."/break.predict";
my $break_filename = $output_directory."/break";
my $seq_filename = $output_directory."/seq";
my $span_pval_filename = $output_directory."/span.pval";

my %splitr_span_pval;
my %denovo_span_pval;
read_span_pval($splitr_span_pval_filename, \%splitr_span_pval);
read_span_pval($denovo_span_pval_filename, \%denovo_span_pval);

open BP, ">".$break_predict_filename or die "Error: Unable to find $break_predict_filename: $!\n";

my %splitr_cluster_ids;
my %denovo_cluster_ids;
my %all_cluster_ids = (%splitr_span_pval, %denovo_span_pval);
foreach my $cluster_id (keys %all_cluster_ids)
{
	# Select the breakpoint prediction with the highest spanning pvalue
	my $break_predict;
	if (not defined $denovo_span_pval{$cluster_id})
	{
		$break_predict = "splitr";
		$splitr_cluster_ids{$cluster_id} = 1;
	}
	elsif (not defined $splitr_span_pval{$cluster_id})
	{
		$break_predict = "denovo";
		$denovo_cluster_ids{$cluster_id} = 1;
	}
	elsif ($denovo_span_pval{$cluster_id}{pvalue} > $splitr_span_pval{$cluster_id}{pvalue})
	{
		$break_predict = "denovo";
		$denovo_cluster_ids{$cluster_id} = 1;
	}
	else
	{
		$break_predict = "splitr";
		$splitr_cluster_ids{$cluster_id} = 1;
	}
	
	print BP $cluster_id."\t".$break_predict."\n";
}

close BP;

my $break_file;
open $break_file, ">".$break_filename or die "Error: Unable to find $break_filename: $!\n";
filter_by_cluster_id($splitr_break_filename, $break_file, \%splitr_cluster_ids);
filter_by_cluster_id($denovo_break_filename, $break_file, \%denovo_cluster_ids);
close $break_file;

my $seq_file;
open $seq_file, ">".$seq_filename or die "Error: Unable to find $seq_filename: $!\n";
filter_by_cluster_id($splitr_seq_filename, $seq_file, \%splitr_cluster_ids);
filter_by_cluster_id($denovo_seq_filename, $seq_file, \%denovo_cluster_ids);
close $seq_file;

my $span_pval_file;
open $span_pval_file, ">".$span_pval_filename or die "Error: Unable to find $span_pval_filename: $!\n";
filter_by_cluster_id($splitr_span_pval_filename, $span_pval_file, \%splitr_cluster_ids);
filter_by_cluster_id($denovo_span_pval_filename, $span_pval_file, \%denovo_cluster_ids);
close $span_pval_file;

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

sub filter_by_cluster_id
{
	my $in_filename = shift;
	my $out_file = shift;
	my $cluster_ids = shift;
	
	open IN, $in_filename or die "Error: Unable to find $in_filename: $!\n";
	while (<IN>)
	{
		my $line = $_;
		
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		
		if (defined $cluster_ids->{$cluster_id})
		{
			print $out_file $line;
		}
	}
	close IN;
}

