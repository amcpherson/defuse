#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
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
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --name      Library Name (default: Output Directory Suffix)\n";
push @usage, "  -s, --submit    Submitter Type (default: direct)\n";
push @usage, "  -p, --parallel  Maximum Number of Parallel Jobs (default: 1)\n";

my $help;
my $config_filename;
my $output_directory;
my $library_name;
my $submitter_type;
my $max_parallel;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
	'name=s'      => \$library_name,
	'submit=s'    => \$submitter_type,
	'parallel=s'  => \$max_parallel,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
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
my $rscript_bin           = $config->get_value("rscript_bin");
my $split_min_anchor      = $config->get_value("split_min_anchor");
my $cov_samp_density      = $config->get_value("covariance_sampling_density");
my $denovo_assembly       = $config->get_value("denovo_assembly");
my $remove_job_files      = $config->get_value("remove_job_files");
my $positive_controls     = $config->get_value("positive_controls");

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
my $alignjob_script = "$scripts_directory/alignjob.pl";
my $filter_bowtie_script = "$scripts_directory/filter_bowtie.pl";
my $filter_reference_script = "$scripts_directory/filter_bowtie_reference.pl";
my $read_stats_script = "$scripts_directory/read_stats.pl";
my $expression_script = "$scripts_directory/calculate_expression_simple.pl";
my $split_fastq_script = "$scripts_directory/split_fastq.pl";
my $split_regions_script = "$scripts_directory/split_regions.pl";
my $split_candidate_reads_script = "$scripts_directory/split_candidate_reads.pl";
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
my $findcandidatereads_bin = "$tools_directory/findcandidatereads";
my $dosplitalign_bin = "$tools_directory/dosplitalign";
my $evalsplitalign_bin = "$tools_directory/evalsplitalign";
my $calc_span_pvalues_script = "$scripts_directory/calc_span_pvalues.pl";
my $calc_split_pvalues_script = "$scripts_directory/calc_split_pvalues.pl";
my $calc_map_stats_script = "$scripts_directory/calculate_mapping_stats.pl";
my $annotate_fusions_script = "$scripts_directory/annotate_fusions.pl";
my $filter_fusions_script = "$scripts_directory/filter_fusions.pl";
my $coallate_fusions_script = "$scripts_directory/coallate_fusions.pl";
my $adaboost_rscript = "$rscript_bin ".$scripts_directory."/run_adaboost.R";

mkdir $output_directory if not -e $output_directory;

my $job_directory = $output_directory."/jobs";
my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/defuse";

mkdir $log_directory if not -e $log_directory;
mkdir $job_directory if not -e $job_directory;

my $runner = cmdrunner->new();
$runner->name("reannotate");
$runner->prefix($log_prefix);
$runner->maxparallel($max_parallel);

print "Annotating fusions\n";
my $splitr_span_pval = $output_directory."/splitr.span.pval";
my $denovo_span_pval = $output_directory."/denovo.span.pval";
my $splitr_split_pval = $output_directory."/splitr.split.pval";
my $clusters_sc_sam = $output_directory."/clusters.sc.sam";
my $annotations_filename = $output_directory."/annotations.txt";
my $mapping_stats_filename = $output_directory."/mapping.stats";
$runner->run("$annotate_fusions_script -c $config_filename -o $output_directory -n $library_name > #>1", [$splitr_span_pval,$denovo_span_pval,$splitr_split_pval], [$annotations_filename]);
$runner->run("$calc_map_stats_script -c $config_filename -o $output_directory > #>1", [$clusters_sc_sam], [$mapping_stats_filename]);

print "Filtering fusions\n";
my $filtered_filename = $output_directory."/filtered.txt";
$runner->run("$filter_fusions_script -c $config_filename -o $output_directory > #>1", [$annotations_filename], [$filtered_filename]);

print "Coallating fusions\n";
my $results_filename = $output_directory."/results.txt";
my $filtered_results_filename = $output_directory."/results.filtered.txt";
$runner->run("$coallate_fusions_script -c $config_filename -o $output_directory -l #<1 > #>1", [$annotations_filename], [$results_filename]);
$runner->run("$coallate_fusions_script -c $config_filename -o $output_directory -l #<1 > #>1", [$filtered_filename], [$filtered_results_filename]);

print "Running adaboost classifier\n";
my $results_classify = $output_directory."/results.classify.txt";
$runner->run("$adaboost_rscript $positive_controls #<1 #>1", [$results_filename], [$results_classify]);

print "Success\n";
$status = "success";

