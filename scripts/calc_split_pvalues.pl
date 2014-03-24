#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Temp qw[tempfile tempdir];

use lib dirname($0);
use configdata;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Calculate spanning fragment length pvalue\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -s, --seqs      Split Sequence Predictions Filename\n";
push @usage, "  -v, --vstat     Covariance Statistics Filename\n";

my $help;
my $config_filename;
my $seqs_filename;
my $vstat_filename;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'seqs=s'      => \$seqs_filename,
	'vstat=s'     => \$vstat_filename,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $seqs_filename or die @usage;
defined $vstat_filename or die @usage;

my $config = configdata->new();
$config->read($config_filename);

my $r_bin = $config->get_value("r_bin");

# Read in the read stats
my %covariance_values;
get_stats($vstat_filename, \%covariance_values);

my $split_pos_covariance = $covariance_values{"split_pos_covariance"};
my $split_min_covariance = $covariance_values{"split_min_covariance"};

open SEQ, $seqs_filename or die "Error: Unable to find $seqs_filename: $!\n";
while (<SEQ>)
{
	chomp;
	my @fields = split /\t/;
	
	my $cluster_id = $fields[0];
	my $split_count = $fields[3];
	my $split_pos_average = $fields[4];
	my $split_min_average = $fields[5];
	
	my $pval_result_split_pos = `echo '2 * pnorm(-1.0 * abs($split_pos_average - 0.5)/(sqrt($split_pos_covariance+(1/(12*$split_count)))))' | $r_bin --vanilla -q`;
	$pval_result_split_pos =~ /\[1\]\s+(.+)/ or die "Error: Unable to interpret pvalue from R\n";
	my $split_pos_pvalue = $1;

	my $pval_result_split_min = `echo 'pnorm(($split_min_average - 0.5)/(sqrt($split_min_covariance+(1/(12*$split_count)))))' | $r_bin --vanilla -q`;
	$pval_result_split_min =~ /\[1\]\s+(.+)/ or die "Error: Unable to interpret pvalue from R\n";
	my $split_min_pvalue = $1;

	print $cluster_id."\t".$split_pos_pvalue."\t".$split_min_pvalue."\n";
}
close SEQ;

sub get_stats
{
	my $stats_filename = shift;
	my $stats_outref = shift;
	
	open STATS, $stats_filename or die "Error: Unable to open $stats_filename\n";
	my @stats = <STATS>;
	chomp(@stats);
	close STATS;

	scalar @stats == 2 or die "Error: Stats file $stats_filename does not have 2 lines\n";

	my @keys = split /\t/, $stats[0];
	my @values = split /\t/, $stats[1];

	scalar @keys == scalar @values or die "Error: Stats file $stats_filename with column mismatch\n";

	foreach my $stat_index (0..$#keys)
	{
		my $key = $keys[$stat_index];
		my $value = $values[$stat_index];

		$stats_outref->{$key} = $value;
	}
}


