#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use File::Spec;

use lib dirname($0);
use configdata;
use cmdrunner;
use parsers;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Run the defuse pipeline for fusion discovery.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --data      Source Data Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --name      Library Name (default: Output Directory Suffix)\n";
push @usage, "  -l, --local     Job Local Directory (default: Output Directory)\n";
push @usage, "  -s, --submit    Submitter Type (default: direct)\n";
push @usage, "  -p, --parallel  Maximum Number of Parallel Jobs (default: 1)\n";

my $help;
my $config_filename;
my $source_directory;
my $output_directory;
my $library_name;
my $joblocal_directory;
my $submitter_type;
my $max_parallel;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'data=s'      => \$source_directory,
	'output=s'    => \$output_directory,
	'name=s'      => \$library_name,
	'local=s'     => \$joblocal_directory,
	'submit=s'    => \$submitter_type,
	'parallel=s'  => \$max_parallel,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
defined $source_directory or die @usage;
defined $output_directory or die @usage;

mkdir $output_directory if not -d $output_directory;

$output_directory = abs_path($output_directory);
$config_filename = abs_path($config_filename);

# Guess library name
if (not defined $library_name)
{
	$library_name = "";
	my @output_splitdir = File::Spec->splitdir($output_directory);
	while ($library_name eq "" and scalar @output_splitdir > 0)
	{
		$library_name = pop(@output_splitdir);
	}
	defined $library_name or die "Error: Unable to infer library name from output director $output_directory\n";
}

# Job local directory defaults to temp subdirectory of output directory
if (not defined $joblocal_directory)
{
	$joblocal_directory = $output_directory."/tmp";
	mkdir $joblocal_directory if not -d $joblocal_directory;
}

# Submitter type defaults to direct
if (not defined $submitter_type)
{
	$submitter_type = "direct";
}

# Max parallel defaults to 1 for direct submitter
if (not defined $max_parallel)
{
	if ($submitter_type eq "direct")
	{
		$max_parallel = 1;
	}
	else
	{
		$max_parallel = 200;
	}
}

my $config = configdata->new();
$config->read($config_filename);

# Config values
my $genome_fasta          = $config->get_value("genome_fasta");
my $cdna_gene_fasta       = $config->get_value("cdna_gene_fasta");
my $cdna_gene_regions     = $config->get_value("cdna_gene_regions");
my $cdna_fasta            = $config->get_value("cdna_fasta");
my $cdna_regions          = $config->get_value("cdna_regions");
my $gene_tran_list        = $config->get_value("gene_tran_list");
my $ig_gene_list          = $config->get_value("ig_gene_list");
my $max_insert_size       = $config->get_value("max_insert_size");
my $span_count_threshold  = $config->get_value("span_count_threshold");
my $clustering_precision  = $config->get_value("clustering_precision");
my $scripts_directory     = $config->get_value("scripts_directory");
my $tools_directory       = $config->get_value("tools_directory");
my $reads_per_job         = $config->get_value("reads_per_job");
my $regions_per_job       = $config->get_value("regions_per_job");
my $gene_info_list        = $config->get_value("gene_info_list");
my $samtools_bin          = $config->get_value("samtools_bin");
my $split_min_anchor      = $config->get_value("split_min_anchor");
my $cov_samp_density      = $config->get_value("covariance_sampling_density");

my $cdna_fasta_fai        = $cdna_fasta.".fai";

my $mailto;
if ($config->has_value("mailto"))
{
	$mailto = $config->get_value("mailto");
}

my $status = "failure";
sub mailme
{
	return if not defined $mailto;

	my $text = "Fusion analysis of library $library_name finished with status $status";

	print "Attempting to mail $mailto the result\n";

	system "echo '$text' | mail -s '[AUTO] $text' $mailto";
}

my $main_pid = $$;
END
{
	my $retcode = $?;
	if (defined $main_pid and $$ == $main_pid)
	{
		mailme();
	}
	$? = $retcode;
}

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
verify_file_exists($genome_fasta);
verify_file_exists($cdna_gene_fasta);
verify_file_exists($cdna_fasta);
verify_file_exists($cdna_fasta_fai);
verify_file_exists($gene_tran_list);
verify_file_exists($ig_gene_list);

# Ensure required directories exist
verify_directory_exists($scripts_directory);

# Script and binary paths
my $retreive_fastq_script = "$scripts_directory/retreive_fastq.pl";
my $remove_end_sed_command = "sed 's#/[12]\$##g'";
my $filter_bowtie_script = "$scripts_directory/filter_bowtie.pl";
my $filter_reference_script = "$scripts_directory/filter_bowtie_reference.pl";
my $read_stats_script = "$scripts_directory/read_stats.pl";
my $expression_script = "$scripts_directory/calculate_expression_simple.pl";
my $split_fastq_script = "$scripts_directory/split_fastq.pl";
my $split_regions_script = "$scripts_directory/split_regions.pl";
my $get_align_regions_script = "$scripts_directory/get_align_regions.pl";
my $merge_read_stats_script = "$scripts_directory/merge_read_stats.pl";
my $merge_expression_script = "$scripts_directory/merge_expression.pl";
my $split_align_merge_script = "$scripts_directory/splitalignmerge.pl";
my $split_seq_merge_script = "$scripts_directory/splitseqmerge.pl";
my $get_cluster_alignments_script = "$scripts_directory/get_cluster_alignments.pl";
my $bowtie2sam_script = "$scripts_directory/bowtie2sam.pl";
my $calccov_bin = "$tools_directory/calccov";
my $maxvalclust_bin = "$tools_directory/maximalvalidclusters";
my $setcover_bin = "$tools_directory/setcover";
my $split_seq_bin = "$tools_directory/splitseq";
my $denovo_seq_bin = "$tools_directory/denovoseq";
my $calc_span_pvalues_script = "$scripts_directory/calc_span_pvalues.pl";
my $calc_split_pvalues_script = "$scripts_directory/calc_split_pvalues.pl";
my $annotate_fusions_script = "$scripts_directory/annotate_fusions.pl";
my $filter_fusions_script = "$scripts_directory/filter_fusions.pl";
my $coallate_fusions_script = "$scripts_directory/coallate_fusions.pl";

mkdir $output_directory if not -e $output_directory;

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/defuse";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("defuse");
$runner->prefix($log_prefix);
$runner->maxparallel($max_parallel);

# Fastq files and prefix for index and name map
my $reads_prefix = $output_directory."/reads";
my $reads_end_1_fastq = $reads_prefix.".1.fastq";
my $reads_end_2_fastq = $reads_prefix.".2.fastq";

# Retrieve Fastq Files
$runner->run("$retreive_fastq_script -c $config_filename -d $source_directory -o $output_directory -s $submitter_type", [], [$reads_end_1_fastq,$reads_end_2_fastq]);

# Products of the alignment process
my $read_stats = $output_directory."/concordant.read.stats";
my $expression = $output_directory."/expression.txt";
my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
my $discordant_aligned_bam = $output_directory."/discordant.aligned.bam";
my $discordant_unaligned_bam = $output_directory."/discordant.unaligned.bam";

# Run the alignment process
$runner->run("$split_align_merge_script -c $config_filename -l $joblocal_directory -o $output_directory -n $library_name -s $submitter_type -p $max_parallel", [$reads_end_1_fastq,$reads_end_2_fastq], [$read_stats,$cdna_pair_bam,$discordant_aligned_bam,$discordant_unaligned_bam,$expression]);

# Run remaining commands on cluster nodes
$runner->submitter($submitter_type);
$runner->jobmem(30000000000);

# Index each bam file
my $cdna_pair_bam_bai = $cdna_pair_bam.".bai";
my $discordant_aligned_bam_bai = $discordant_aligned_bam.".bai";
my $discordant_unaligned_bam_bai = $discordant_unaligned_bam.".bai";
$runner->padd("$samtools_bin index #<1", [$cdna_pair_bam], [$cdna_pair_bam_bai]);
$runner->padd("$samtools_bin index #<1", [$discordant_aligned_bam], [$discordant_aligned_bam_bai]);
$runner->padd("$samtools_bin index #<1", [$discordant_unaligned_bam], [$discordant_unaligned_bam_bai]);
$runner->prun();

# Read in the read stats
my %read_stat_values;
parsers::get_stats($read_stats, \%read_stat_values);

my $read_length_min = $read_stat_values{"readlength_min"};
my $read_length_max = $read_stat_values{"readlength_max"};
my $fragment_mean = $read_stat_values{"fraglength_mean"};
my $fragment_stddev = $read_stat_values{"fraglength_stddev"};
my $fragment_max = int($fragment_mean + 3 * $fragment_stddev);

print "Read Stats\n";
print "\tFragment mean $fragment_mean stddev $fragment_stddev\n";
print "\tRead length min $read_length_min max $read_length_max\n";

print "Finding paired end alignment covariance\n";
my $cov_stats = $output_directory."/covariance.stats";
$runner->run("$calccov_bin -m $fragment_max -a $split_min_anchor -d $cov_samp_density -g $gene_tran_list -c #<1 -v #>1", [$cdna_pair_bam], [$cov_stats]);

print "Generating maximal valid clusters\n";
my $clusters = $output_directory."/clusters.txt";
$runner->run("$maxvalclust_bin -m $span_count_threshold -a $clustering_precision -u $fragment_mean -s $fragment_stddev -d #<1 -g $gene_tran_list -c #>1", [$discordant_aligned_bam], [$clusters]);

print "Generating maximum parsimony solution\n";
my $clusters_sc = $output_directory."/clusters.sc.txt";
$runner->run("$setcover_bin -m $span_count_threshold -c #<1 -o #>1", [$clusters], [$clusters_sc]);

print "Generating sam format clusters\n";
my $clusters_sc_sam = $output_directory."/clusters.sc.sam";
$runner->run("$samtools_bin view #<1 | $get_cluster_alignments_script #<2 > #>1", [$discordant_aligned_bam,$clusters_sc], [$clusters_sc_sam]);

print "Generating spanning alignment regions file\n";
my $clusters_sc_regions = $output_directory."/clusters.sc.regions";
$runner->run("$get_align_regions_script < #<1 > #>1", [$clusters_sc_sam], [$clusters_sc_regions]);

print "Breakpoint sequence prediction\n";
my $splitr_break = $output_directory."/splitr.break";
my $splitr_seq = $output_directory."/splitr.seq";
my $splitr_log = $output_directory."/splitr.log";
my $denovo_break = $output_directory."/denovo.break";
my $denovo_seq = $output_directory."/denovo.seq";
my $denovo_log = $output_directory."/denovo.log";
$runner->submitter("direct");
$runner->run("$split_seq_merge_script -c $config_filename -o $output_directory -s $submitter_type -p $max_parallel", [$clusters_sc_regions,$discordant_aligned_bam,$discordant_unaligned_bam], [$splitr_break,$splitr_seq,$denovo_break,$denovo_seq]);
$runner->submitter($submitter_type);

print "Calculating spanning pvalues\n";
my $splitr_span_pval = $output_directory."/splitr.span.pval";
my $denovo_span_pval = $output_directory."/denovo.span.pval";
$runner->padd("$calc_span_pvalues_script -c $config_filename -p #<1 -b #<2 -s #<3 -r #<4 -v #<5 > #>1", [$clusters_sc_sam,$splitr_break,$splitr_seq,$read_stats,$cov_stats], [$splitr_span_pval]);
$runner->padd("$calc_span_pvalues_script -c $config_filename -p #<1 -b #<2 -s #<3 -r #<4 -v #<5 > #>1", [$clusters_sc_sam,$denovo_break,$denovo_seq,$read_stats,$cov_stats], [$denovo_span_pval]);
$runner->prun();

# Calculate split read pvalues
my $splitr_split_pval = $output_directory."/splitr.split.pval";
$runner->run("$calc_split_pvalues_script -c $config_filename -s #<1 -v #<2 > #>1", [$splitr_seq,$cov_stats], [$splitr_split_pval]);

# Annotate fusions
my $annotations_filename = $output_directory."/annotations.txt";
$runner->run("$annotate_fusions_script -c $config_filename -o $output_directory -n $library_name > #>1", [$clusters_sc_sam], [$annotations_filename]);

# Run the filtering script
my $filtered_filename = $output_directory."/filtered.txt";
$runner->run("$filter_fusions_script -c $config_filename -o $output_directory > #>1", [$annotations_filename], [$filtered_filename]);

# Coallate the filtered and unfiltered results
my $results_filename = $output_directory."/results.txt";
my $filtered_results_filename = $output_directory."/results.filtered.txt";
$runner->run("$coallate_fusions_script -c $config_filename -o $output_directory -l #<1 > #>1", [$annotations_filename], [$results_filename]);
$runner->run("$coallate_fusions_script -c $config_filename -o $output_directory -l #<1 > #>1", [$filtered_filename], [$filtered_results_filename]);

$status = "success";

