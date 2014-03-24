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

$| = 1;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Find spanning alignments in parallel.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -l, --local     Job Local Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --name      Library Name\n";
push @usage, "  -s, --submit    Submitter Type\n";
push @usage, "  -p, --parallel  Maximum Number of Parallel Jobs (default: 1)\n";

my $help;
my $config_filename;
my $joblocal_directory;
my $output_directory;
my $library_name;
my $submitter_type;
my $max_parallel;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'local=s'     => \$joblocal_directory,
	'output=s'    => \$output_directory,
	'name=s'      => \$library_name,
	'submit=s'    => \$submitter_type,
	'parallel=s'  => \$max_parallel,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
defined $joblocal_directory or die @usage;
defined $output_directory or die @usage;
defined $library_name or die @usage;
defined $submitter_type or die @usage;
defined $max_parallel or die @usage;

my $config = configdata->new();
$config->read($config_filename);

# Config values
my $scripts_directory     = $config->get_value("scripts_directory");
my $reads_per_job         = $config->get_value("reads_per_job");

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

# Ensure required directories exist
verify_directory_exists($scripts_directory);
verify_directory_exists($output_directory);

# Script and binary paths
my $alignjob_script = "$scripts_directory/alignjob.pl";
my $split_fastq_script = "$scripts_directory/split_fastq.pl";
my $merge_read_stats_script = "$scripts_directory/merge_read_stats.pl";
my $merge_expression_script = "$scripts_directory/merge_expression.pl";

$output_directory = abs_path($output_directory);

my $job_directory = abs_path($output_directory."/jobs");

mkdir $job_directory if not -e $job_directory;

my $reads_end_1_fastq = $output_directory."/reads.1.fastq";
my $reads_end_2_fastq = $output_directory."/reads.2.fastq";

-e $reads_end_1_fastq and -e $reads_end_2_fastq or die "Error: Unable to find fastq files\n";

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/splitalignmerge";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("splitalignmerge");
$runner->prefix($log_prefix);
$runner->submitter($submitter_type);
$runner->maxparallel($max_parallel);
$runner->jobmem(3000000000);

my $reads_end_1_split_filenames = $job_directory."/reads.1.split";
my $reads_end_2_split_filenames = $job_directory."/reads.2.split";

my $reads_end_1_split_fastq_template = $job_directory."/reads_%.1.fastq";
my $reads_end_2_split_fastq_template = $job_directory."/reads_%.2.fastq";

$runner->run("$split_fastq_script #<1 $reads_per_job $reads_end_1_split_fastq_template > #>1", [$reads_end_1_fastq], [$reads_end_1_split_filenames]);
$runner->run("$split_fastq_script #<1 $reads_per_job $reads_end_2_split_fastq_template > #>1", [$reads_end_2_fastq], [$reads_end_2_split_filenames]);

open SF1, $reads_end_1_split_filenames or die "Error: Unable to open $reads_end_1_split_filenames\n";
my @split_filenames_1 = <SF1>;
close SF1;
chomp(@split_filenames_1);

open SF2, $reads_end_2_split_filenames or die "Error: Unable to open $reads_end_2_split_filenames\n";
my @split_filenames_2 = <SF2>;
close SF2;
chomp(@split_filenames_2);

die "Error: Split inconsistent\n" unless scalar @split_filenames_1 == scalar @split_filenames_2;

# Names of job products
my @all_jobs_read_stats;
my @all_jobs_cdna_pair_bam;
my @all_jobs_discordant_aligned_bam;
my @all_jobs_discordant_unaligned_bam;
my @all_jobs_expression;

foreach my $split_filename_index (0..$#split_filenames_1)
{
	my $reads_end_1_split_fastq = $split_filenames_1[$split_filename_index];
	my $reads_end_2_split_fastq = $split_filenames_2[$split_filename_index];
	
	$reads_end_1_split_fastq =~ /\/(reads_\d+)\.1\.fastq$/;
	my $split_prefix_1 = $1;
	
	$reads_end_2_split_fastq =~ /\/(reads_\d+)\.2\.fastq$/;
	my $split_prefix_2 = $1;
	
	die "Error: Split inconsistent\n" unless defined $split_prefix_1 and $split_prefix_1 eq $split_prefix_2;
	my $split_prefix = $split_prefix_1;

	# Prefix for jobs
	my $library_jobs_split_prefix = $job_directory."/".$split_prefix;
	
	# Names of job products
	my $job_read_stats = $library_jobs_split_prefix.".concordant.read.stats";
	my $job_cdna_pair_bam = $library_jobs_split_prefix.".cdna.pair.bam";
	my $job_discordant_aligned_bam = $library_jobs_split_prefix.".discordant.aligned.bam";
	my $job_discordant_unaligned_bam = $library_jobs_split_prefix.".discordant.unaligned.bam";
	my $job_expression = $library_jobs_split_prefix.".expression.txt";
	
	push @all_jobs_read_stats, $job_read_stats;
	push @all_jobs_cdna_pair_bam, $job_cdna_pair_bam;
	push @all_jobs_discordant_aligned_bam, $job_discordant_aligned_bam;
	push @all_jobs_discordant_unaligned_bam, $job_discordant_unaligned_bam;
	push @all_jobs_expression, $job_expression;

	my $job_cmd = $alignjob_script." ";
	$job_cmd .= "-c ".$config_filename." ";
	$job_cmd .= "-j ".$job_directory." ";
	$job_cmd .= "-l ".$joblocal_directory." ";
	$job_cmd .= "-o ".$output_directory." ";
	$job_cmd .= "-n ".$library_name." ";
	$job_cmd .= "-p ".$split_prefix." ";
	
	$runner->padd("$job_cmd", [$reads_end_1_split_fastq, $reads_end_2_split_fastq], [$job_read_stats, $job_cdna_pair_bam, $job_discordant_aligned_bam, $job_discordant_unaligned_bam, $job_expression]); 
}
$runner->prun();

# Merge products of the alignment process
my $read_stats = $output_directory."/concordant.read.stats";
my $expression = $output_directory."/expression.txt";
my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
my $discordant_aligned_bam = $output_directory."/discordant.aligned.bam";
my $discordant_unaligned_bam = $output_directory."/discordant.unaligned.bam";

# Maximum number of files to merge at once
my $merge_max = 100;

# Merge merge_max at a time
my @merge_jobs_read_stats;
my @merge_jobs_expression;
my @merge_jobs_cdna_pair_bam;
my @merge_jobs_discordant_aligned_bam;
my @merge_jobs_discordant_unaligned_bam;
my @all_intermediate_read_stats;
my @all_intermediate_expression;
my @all_intermediate_cdna_pair_bam;
my @all_intermediate_discordant_aligned_bam;
my @all_intermediate_discordant_unaligned_bam;
foreach my $job_index (0..$#all_jobs_read_stats)
{
	my $job_read_stats = $all_jobs_read_stats[$job_index];
	my $job_expression = $all_jobs_expression[$job_index];
	my $job_cdna_pair_bam = $all_jobs_cdna_pair_bam[$job_index];
	my $job_discordant_aligned_bam = $all_jobs_discordant_aligned_bam[$job_index];
	my $job_discordant_unaligned_bam = $all_jobs_discordant_unaligned_bam[$job_index];
	
	push @merge_jobs_read_stats, $job_read_stats;
	push @merge_jobs_expression, $job_expression;
	push @merge_jobs_cdna_pair_bam, $job_cdna_pair_bam;
	push @merge_jobs_discordant_aligned_bam, $job_discordant_aligned_bam;
	push @merge_jobs_discordant_unaligned_bam, $job_discordant_unaligned_bam;
	
	if (scalar @merge_jobs_read_stats == $merge_max or $job_index == $#all_jobs_read_stats)
	{
		my $intermediate_index = scalar(@all_intermediate_read_stats) + 1;
		my $intermediate_prefix = "intermediate.".$intermediate_index;
		
		my $intermediate_read_stats = $job_directory."/$intermediate_prefix.concordant.read.stats";
		my $intermediate_expression = $job_directory."/$intermediate_prefix.expression.txt";
		my $intermediate_cdna_pair_bam = $job_directory."/$intermediate_prefix.cdna.pair.bam";
		my $intermediate_discordant_aligned_bam = $job_directory."/$intermediate_prefix.discordant.aligned.bam";
		my $intermediate_discordant_unaligned_bam = $job_directory."/$intermediate_prefix.discordant.unaligned.bam";
		
		# Merge read statistics
		$runner->padd("$merge_read_stats_script #<A > #>1", [@merge_jobs_read_stats], [$intermediate_read_stats]);
		
		# Merge expression statistics
		$runner->padd("$merge_expression_script #<A > #>1", [@merge_jobs_expression], [$intermediate_expression]);
		
		# Merge bam files or copy if theres only one
		if (scalar @merge_jobs_cdna_pair_bam == 1)
		{
			rename $merge_jobs_cdna_pair_bam[0], $intermediate_cdna_pair_bam;
			rename $merge_jobs_discordant_aligned_bam[0], $intermediate_discordant_aligned_bam;
			rename $merge_jobs_discordant_unaligned_bam[0], $intermediate_discordant_unaligned_bam;
		}
		else
		{
			$runner->padd("samtools merge #>1 #<A", [@merge_jobs_cdna_pair_bam], [$intermediate_cdna_pair_bam]);
			$runner->padd("samtools merge #>1 #<A", [@merge_jobs_discordant_aligned_bam], [$intermediate_discordant_aligned_bam]);
			$runner->padd("samtools merge #>1 #<A", [@merge_jobs_discordant_unaligned_bam], [$intermediate_discordant_unaligned_bam]);
		}
		
		push @all_intermediate_read_stats, $intermediate_read_stats;
		push @all_intermediate_expression, $intermediate_expression;
		push @all_intermediate_cdna_pair_bam, $intermediate_cdna_pair_bam;
		push @all_intermediate_discordant_aligned_bam, $intermediate_discordant_aligned_bam;
		push @all_intermediate_discordant_unaligned_bam, $intermediate_discordant_unaligned_bam;
		
		@merge_jobs_read_stats = ();
		@merge_jobs_expression = ();
		@merge_jobs_cdna_pair_bam = ();
		@merge_jobs_discordant_aligned_bam = ();
		@merge_jobs_discordant_unaligned_bam = ();	
	}
}
$runner->prun();

# Merge intermediate read statistics
$runner->run("$merge_read_stats_script #<A > #>1", [@all_intermediate_read_stats], [$read_stats]);

# Merge intermediate expression statistics
$runner->run("$merge_expression_script #<A > #>1", [@all_intermediate_expression], [$expression]);

# Merge bam files or copy if theres only one
if (scalar @all_intermediate_discordant_aligned_bam == 1)
{
	rename $all_intermediate_cdna_pair_bam[0], $cdna_pair_bam;
	rename $all_intermediate_discordant_aligned_bam[0], $discordant_aligned_bam;
	rename $all_intermediate_discordant_unaligned_bam[0], $discordant_unaligned_bam;
}
else
{
	$runner->run("samtools merge #>1 #<A", [@all_intermediate_cdna_pair_bam], [$cdna_pair_bam]);
	$runner->run("samtools merge #>1 #<A", [@all_intermediate_discordant_aligned_bam], [$discordant_aligned_bam]);
	$runner->run("samtools merge #>1 #<A", [@all_intermediate_discordant_unaligned_bam], [$discordant_unaligned_bam]);
}

unlink @all_jobs_read_stats;
unlink @all_jobs_expression;
unlink @all_jobs_cdna_pair_bam;
unlink @all_jobs_discordant_aligned_bam;
unlink @all_jobs_discordant_unaligned_bam;
unlink @all_intermediate_read_stats;
unlink @all_intermediate_expression;
unlink @all_intermediate_cdna_pair_bam;
unlink @all_intermediate_discordant_aligned_bam;
unlink @all_intermediate_discordant_unaligned_bam;
unlink @split_filenames_1;
unlink @split_filenames_2;
unlink $reads_end_1_split_filenames;
unlink $reads_end_2_split_filenames;

print "Finished Split-Align-Merge\n";

