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
my $scripts_directory     = $config->get_value("scripts_directory");
my $tools_directory       = $config->get_value("tools_directory");
my $samtools_bin          = $config->get_value("samtools_bin");
my $rscript_bin           = $config->get_value("rscript_bin");
my $remove_job_files      = $config->get_value("remove_job_files");
my $positive_controls     = $config->get_value("positive_controls");
my $probability_threshold = $config->get_value("probability_threshold");

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

# Script and binary paths
my $calc_map_stats_script = "$scripts_directory/calculate_mapping_stats.pl";
my $annotate_fusions_script = "$scripts_directory/annotate_fusions.pl";
my $filter_fusions_script = "$scripts_directory/filter_fusions.pl";
my $coallate_fusions_script = "$scripts_directory/coallate_fusions.pl";
my $adaboost_rscript = "$rscript_bin ".$scripts_directory."/run_adaboost.R";
my $filter_script = "$scripts_directory/filter.pl";

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
$runner->submitter($submitter_type);

print "Annotating fusions\n";
my $splitr_span_pval = $output_directory."/splitr.span.pval";
my $denovo_span_pval = $output_directory."/denovo.span.pval";
my $splitr_split_pval = $output_directory."/splitr.split.pval";
my $clusters_sc_sam = $output_directory."/clusters.sc.sam";
my $annotations_filename = $output_directory."/annotations.txt";
my $mapping_stats_filename = $output_directory."/mapping.stats";
$runner->run("$annotate_fusions_script -c $config_filename -o $output_directory -n $library_name > #>1", [$splitr_span_pval,$denovo_span_pval,$splitr_split_pval], [$annotations_filename]);
$runner->run("$calc_map_stats_script -c $config_filename -o $output_directory > #>1", [$clusters_sc_sam], [$mapping_stats_filename]);

print "Coallating fusions\n";
my $results_filename = $output_directory."/results.txt";
$runner->run("$coallate_fusions_script -c $config_filename -o $output_directory -l #<1 > #>1", [$annotations_filename], [$results_filename]);

print "Running adaboost classifier\n";
my $results_classify = $output_directory."/results.classify.txt";
$runner->run("$adaboost_rscript $positive_controls #<1 #>1", [$results_filename], [$results_classify]);

print "Filtering fusions\n";
my $filtered_results_filename = $output_directory."/results.filtered.txt";
$runner->run("$filter_script probability '> $probability_threshold' < #<1 > #>1", [$results_classify], [$filtered_results_filename]);

print "Success\n";
$status = "success";

