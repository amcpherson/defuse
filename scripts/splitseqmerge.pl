#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];

use lib dirname($0);
use configdata;
use cmdrunner;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Predict breakpoint sequences in parallel.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -s, --submit    Submitter Type\n";
push @usage, "  -p, --parallel  Maximum Number of Parallel Jobs (default: 1)\n";

my $help;
my $config_filename;
my $output_directory;
my $submitter_type;
my $max_parallel;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
	'submit=s'    => \$submitter_type,
	'parallel=s'  => \$max_parallel,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
defined $output_directory or die @usage;
defined $submitter_type or die @usage;
defined $max_parallel or die @usage;

my $config = configdata->new();
$config->read($config_filename);

# Config values
my $cdna_gene_fasta       = $config->get_value("cdna_gene_fasta");
my $cdna_gene_regions     = $config->get_value("cdna_gene_regions");
my $scripts_directory     = $config->get_value("scripts_directory");
my $tools_directory       = $config->get_value("tools_directory");
my $regions_per_job       = $config->get_value("regions_per_job");

sub verify_file_exists
{
	my $filename = shift;
	
	if (not -e $filename)
	{
		die "Error: Required file $filename does not exist\n";
	}
}

sub verify_directory_exists
{
	my $filename = shift;
	
	if (not -e $filename)
	{
		die "Error: Required directory $filename does not exist\n";
	}
}

# Ensure required files exist
verify_file_exists($cdna_gene_fasta);
verify_file_exists($cdna_gene_regions);

# Ensure required directories exist
verify_directory_exists($scripts_directory);
verify_directory_exists($output_directory);

# Script and binary paths
my $split_regions_script = "$scripts_directory/split_regions.pl";
my $split_seq_bin = "$tools_directory/splitseq";
my $denovo_seq_bin = "$tools_directory/denovoseq";

$output_directory = abs_path($output_directory);
$config_filename = abs_path($config_filename);

my $job_directory = abs_path($output_directory."/jobs");

mkdir $job_directory if not -e $job_directory;

my $reads_prefix = $output_directory."/reads";
my $reads_end_1_fastq = $reads_prefix.".1.fastq";
my $reads_end_2_fastq = $reads_prefix.".2.fastq";

-e $reads_end_1_fastq and -e $reads_end_2_fastq or die "Error: Unable to find fastq files\n";

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/splitseqmerge";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("splitseqmerge");
$runner->prefix($log_prefix);
$runner->submitter($submitter_type);
$runner->maxparallel($max_parallel);
$runner->jobmem(3000000000);

my $read_stats = $output_directory."/concordant.read.stats";
my $expression = $output_directory."/expression.txt";
my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
my $discordant_aligned_bam = $output_directory."/discordant.aligned.bam";
my $discordant_unaligned_bam = $output_directory."/discordant.unaligned.bam";

# Read in the read stats
my %read_stat_values;
get_stats($read_stats, \%read_stat_values);

my $read_length_min = $read_stat_values{"readlength_min"};
my $read_length_max = $read_stat_values{"readlength_max"};
my $fragment_mean = $read_stat_values{"fraglength_mean"};
my $fragment_stddev = $read_stat_values{"fraglength_stddev"};
my $fragment_max = int($fragment_mean + 3 * $fragment_stddev);

print "Splitting align region pairs\n";
my $clusters_sc_regions = $output_directory."/clusters.sc.regions";
my $split_regions_list = $job_directory."/clusters.sc.regions.splitlist";
my $split_regions_template = $job_directory."/clusters.sc.regions.split_%";
$runner->run("$split_regions_script #<1 $regions_per_job $split_regions_template > #>1", [$clusters_sc_regions], [$split_regions_list]);

open SRL, $split_regions_list or die "Error: Unable to open $split_regions_list\n";
my @split_regions = <SRL>;
close SRL;
chomp(@split_regions);

# Names of job products
my @all_jobs_denovo_break;
my @all_jobs_denovo_seq;
my @all_jobs_denovo_log;

my $common_args = "";
$common_args .= " -m $read_length_min";
$common_args .= " -x $read_length_max";
$common_args .= " -u $fragment_mean";
$common_args .= " -s $fragment_stddev";
$common_args .= " -e $cdna_gene_regions";
$common_args .= " -f $cdna_gene_fasta";
$common_args .= " -d $discordant_aligned_bam";
$common_args .= " -a $discordant_unaligned_bam";
$common_args .= " -r $reads_prefix";

foreach my $split_regions_index (0..$#split_regions)
{
	my $split_regions_filename = $split_regions[$split_regions_index];

	my $split_regions_denovo_break = $split_regions_filename.".denovo.break";
	my $split_regions_denovo_seq = $split_regions_filename.".denovo.seq";
	my $split_regions_denovo_log = $split_regions_filename.".denovo.log";

	push @all_jobs_denovo_break, $split_regions_denovo_break;
	push @all_jobs_denovo_seq, $split_regions_denovo_seq;
	push @all_jobs_denovo_log, $split_regions_denovo_log;

	my $denovo_command = $denovo_seq_bin.$common_args." -i #<1 -b #>1 -q #>2 -l #>3";
	$runner->padd($denovo_command, [$split_regions_filename], [$split_regions_denovo_break, $split_regions_denovo_seq, $split_regions_denovo_log]);
}
$runner->prun();

# Merge products of the breakpoint sequencing process
my $denovo_break = $output_directory."/denovo.break";
my $denovo_seq = $output_directory."/denovo.seq";
my $denovo_log = $output_directory."/denovo.log";

$runner->run("cat #<A > #>1", [@all_jobs_denovo_break], [$denovo_break]);
$runner->run("cat #<A > #>1", [@all_jobs_denovo_seq], [$denovo_seq]);
$runner->run("cat #<A > #>1", [@all_jobs_denovo_log], [$denovo_log]);

unlink $split_regions_list;
unlink @split_regions;
unlink @all_jobs_denovo_break;
unlink @all_jobs_denovo_seq;
unlink @all_jobs_denovo_log;

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

