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
use cmdrunner;

my $main_pid = $$;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Run the bowtie alignment pipeline for fusions.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -l, --local     Job Local Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --name      Library Name\n";
push @usage, "  -j, --job       Job Name\n";
push @usage, "  -p, --prefix    Prefix of Job Input/Output\n";

my $help;
my $config_filename;
my $dataset_directory;
my $local_directory;
my $output_directory;
my $library_name;
my $job_name;
my $job_prefix;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'local=s'     => \$local_directory,
	'output=s'    => \$output_directory,
	'name=s'      => \$library_name,
	'job=s'       => \$job_name,
	'prefix=s'    => \$job_prefix,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
defined $dataset_directory or die @usage;
defined $local_directory or die @usage;
defined $output_directory or die @usage;
defined $library_name or die @usage;
defined $job_name or die @usage;
defined $job_prefix or die @usage;

print "Starting alignjob $job_name\n";

my $source_directory = abs_path("$FindBin::RealBin/../");

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

# Config values
my $gene_models              = $config->get_value("gene_models");
my $cdna_fasta               = $config->get_value("cdna_fasta");
my $genome_fasta             = $config->get_value("genome_fasta");
my $reference_fasta          = $config->get_value("reference_fasta");
my $rrna_fasta               = $config->get_value("rrna_fasta");
my $cdna_regions             = $config->get_value("cdna_regions");
my $ig_gene_list             = $config->get_value("ig_gene_list");
my $max_insert_size          = $config->get_value("max_insert_size");
my $discord_read_trim        = $config->get_value("discord_read_trim");
my $dna_concordant_length    = $config->get_value("dna_concordant_length");
my $scripts_directory        = $config->get_value("scripts_directory");
my $tools_directory          = $config->get_value("tools_directory");
my $bowtie_bin               = $config->get_value("bowtie_bin");
my $samtools_bin             = $config->get_value("samtools_bin");
my $bowtie_threads           = $config->get_value("bowtie_threads");
my @prefilter_fastas         = $config->get_list("prefilter");
my $remove_job_temp_files    = $config->get_value("remove_job_temp_files");
my $bowtie_quals             = $config->get_value("bowtie_quals");
my $bowtie_params            = $config->get_value("bowtie_params");
my $split_min_anchor         = $config->get_value("split_min_anchor");
my $cov_samp_density         = $config->get_value("covariance_sampling_density");
my $max_paired_alignments    = $config->get_value("max_paired_alignments");

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
verify_file_exists($cdna_fasta);
verify_file_exists($reference_fasta);
verify_file_exists($ig_gene_list);

# Ensure required directories exist
verify_directory_exists($scripts_directory);
verify_directory_exists($local_directory);
verify_directory_exists($output_directory);

my $align_pair_bin = "$bowtie_bin $bowtie_params $bowtie_quals --sam-nohead -S -t -X $max_insert_size";
my $align_single_bin = "$bowtie_bin $bowtie_params $bowtie_quals --sam-nohead -S -t -k 100 -m 100";
my $sort_cmd = "sort";

my $filter_fastq_script = "$scripts_directory/filter_fastq.pl";
my $filter_sam_readids_script = "$scripts_directory/filter_sam_readids.pl";
my $filter_sam_genes_script = "$scripts_directory/filter_sam_genes.pl";
my $read_stats_script = "$scripts_directory/read_stats.pl";
my $expression_script = "$scripts_directory/calculate_expression_simple.pl";
my $find_unmappable_script = "$scripts_directory/find_unmappable.pl";
my $find_concordant_gene_script = "$scripts_directory/find_concordant_gene.pl";
my $find_concordant_region_script = "$scripts_directory/find_concordant_region.pl";
my $sam_readids_script = "$scripts_directory/sam_readids.pl";
my $filter_sam_concordant_script = "$scripts_directory/filter_sam_concordant.pl";
my $filter_sam_mapped_script = "$scripts_directory/filter_sam_mapped.pl";
my $trim_fastq_script = "$scripts_directory/trim_fastq.pl";
my $intersect_script = "$scripts_directory/intersect.pl";
my $divide_sam_chr_pairs_script = "$scripts_directory/divide_sam_chr_pairs.pl";
my $match_paired_alignments_script = "$scripts_directory/match_paired_alignments.pl";
my $filter_unmatched_script = "$scripts_directory/filter_unmatched.pl";
my $calccov_bin = "$tools_directory/calccov";

-d $local_directory or die "Error: Could not find local_directory $local_directory\n";

my $local_prefix = $local_directory."/".$library_name.".".$job_name;

# List of local files to remove on control-c
my @local_filenames;

# Cleanup method to remove files
sub cleanup
{
	return if scalar @local_filenames == 0;
	print "Cleaning Up\n";
	if (lc($remove_job_temp_files) eq "yes")
	{
		unlink @local_filenames;
	}
}

# Cleanup when dieing
# Ensure we are passing the correct return code
# Ensure this END block gets called from the main thread/process only
END
{
	my $retcode = $?;
	if ($$ == $main_pid)
	{
		cleanup();
	}
	$? = $retcode;
}

# SIGINT handler
sub inthandler
{
	die "Interrupted\n";
}

$SIG{'INT'} = 'inthandler';

sub get_local_filename
{
	my $suffix = shift;
	my $local_filename = $local_prefix.".".$suffix;
	push @local_filenames, $local_filename;
	return $local_filename;
}

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/".$job_name;

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("alignjob.$job_name");
$runner->prefix($log_prefix);
$runner->submitter("direct");

my $reads_end_1_fastq = $job_prefix.".1.fastq";
my $reads_end_2_fastq = $job_prefix.".2.fastq";

-e $reads_end_1_fastq or die "Error: Unable to find fastq file $reads_end_1_fastq\n";
-e $reads_end_2_fastq or die "Error: Unable to find fastq file $reads_end_2_fastq\n";

# Products of the alignment process
my $read_stats = $job_prefix.".concordant.read.stats";
my $spanlength_samples = $job_prefix.".spanlength.samples";
my $splitpos_samples = $job_prefix.".splitpos.samples";
my $splitmin_samples = $job_prefix.".splitmin.samples";
my $expression = $job_prefix.".expression.txt";
my $cdna_pair_sam = $job_prefix.".cdna.pair.sam";
my $spanning_filelist = $job_prefix.".spanning.filelist";
my $improper_sam = $job_prefix.".improper.sam";

print "Finding concordant alignments to cdna\n";
my $cdna_concordant_readids = get_local_filename("cdna.concordant.readids");
$runner->run("$align_pair_bin $cdna_fasta -1 #<1 -2 #<2 #>1", [$reads_end_1_fastq, $reads_end_2_fastq], [$cdna_pair_sam]);
$runner->run("$filter_sam_concordant_script < #<1 | $sam_readids_script > #>1", [$cdna_pair_sam], [$cdna_concordant_readids]);

print "Finding concordant alignments to dna\n";
my $dna_pair_sam = get_local_filename("dna.pair.sam");
my $dna_concordant_readids = get_local_filename("dna.concordant.readids");
$runner->run("$align_pair_bin $genome_fasta -1 #<1 -2 #<2 #>1", [$reads_end_1_fastq, $reads_end_2_fastq], [$dna_pair_sam]);
$runner->run("$filter_sam_concordant_script < #<1 | $sam_readids_script > #>1", [$dna_pair_sam], [$dna_concordant_readids]);

print "Calculating read statistics\n";
$runner->run("$read_stats_script < #<1 > #>1", [$cdna_pair_sam], [$read_stats]);

print "Calculating covariance samples\n";
my $multi_exon_transcripts_arg = (lc($config->get_value("multi_exon_transcripts_stats")) eq "yes") ? "--multiexon" : "";
$runner->run("$calccov_bin -t $discord_read_trim -a $split_min_anchor -d $cov_samp_density -g $cdna_regions $multi_exon_transcripts_arg -c #<1 -l #>1 -p #>2 -m #>3", [$cdna_pair_sam], [$spanlength_samples, $splitpos_samples, $splitmin_samples]);

print "Calculating expression\n";
$runner->run("$expression_script < #<1 > #>1", [$cdna_pair_sam], [$expression]);

my $pair_concordant_readids = get_local_filename("pair.concordant.readids");
$runner->run("cat #<1 #<2 > #>1", [$cdna_concordant_readids,$dna_concordant_readids], [$pair_concordant_readids]);

print "Trimming and filtering reads\n";
my $trim_reads_end_1_fastq = get_local_filename("trim.1.fastq");
my $trim_reads_end_2_fastq = get_local_filename("trim.2.fastq");
$runner->run("$trim_fastq_script $discord_read_trim < #<1 | $filter_fastq_script -i #<2 > #>1", [$reads_end_1_fastq,$pair_concordant_readids], [$trim_reads_end_1_fastq]);
$runner->run("$trim_fastq_script $discord_read_trim < #<1 | $filter_fastq_script -i #<2 > #>1", [$reads_end_2_fastq,$pair_concordant_readids], [$trim_reads_end_2_fastq]);

print "Finding all cdna alignments\n";
my $cdna_end_1_sam = get_local_filename("cdna.end.1.sam");
my $cdna_end_2_sam = get_local_filename("cdna.end.2.sam");
$runner->run("$align_single_bin $cdna_fasta #<1 #>1", [$trim_reads_end_1_fastq], [$cdna_end_1_sam]);
$runner->run("$align_single_bin $cdna_fasta #<1 #>1", [$trim_reads_end_2_fastq], [$cdna_end_2_sam]);

print "Finding all dna alignments\n";
my $dna_end_1_sam = get_local_filename("dna.end.1.sam");
my $dna_end_2_sam = get_local_filename("dna.end.2.sam");
$runner->run("$align_single_bin $genome_fasta #<1 #>1", [$trim_reads_end_1_fastq], [$dna_end_1_sam]);
$runner->run("$align_single_bin $genome_fasta #<1 #>1", [$trim_reads_end_2_fastq], [$dna_end_2_sam]);

my @concordant_readids;

print "Finding unmappable read ids\n";
my $unmappable_readids = get_local_filename("unmappable.readids");
$runner->run("cat #<A | $find_unmappable_script $gene_models $max_paired_alignments > #>1", [$dna_end_1_sam, $dna_end_2_sam, $cdna_end_1_sam, $cdna_end_2_sam], [$unmappable_readids]);
push @concordant_readids, $unmappable_readids;

print "Finding same gene concordant read ids for cdna\n";
my $gene_concordant_readids = get_local_filename("gene.concordant.readids");
$runner->run("cat #<A | $find_concordant_gene_script $gene_models > #>1", [$dna_end_1_sam, $dna_end_2_sam, $cdna_end_1_sam, $cdna_end_2_sam], [$gene_concordant_readids]);
push @concordant_readids, $gene_concordant_readids;

print "Finding same region concordant read ids for dna\n";
my $region_concordant_readids = get_local_filename("region.concordant.readids");
$runner->run("cat #<A | $find_concordant_region_script $gene_models $dna_concordant_length > #>1", [$dna_end_1_sam, $dna_end_2_sam, $cdna_end_1_sam, $cdna_end_2_sam], [$region_concordant_readids]);
push @concordant_readids, $region_concordant_readids;

print "Finding alignments to rRNA\n";
my $rrna_end_1_sam = get_local_filename("rrna.end.1.sam");
my $rrna_end_2_sam = get_local_filename("rrna.end.2.sam");
$runner->run("$align_single_bin $rrna_fasta #<1 #>1", [$reads_end_1_fastq], [$rrna_end_1_sam]);
$runner->run("$align_single_bin $rrna_fasta #<1 #>1", [$reads_end_2_fastq], [$rrna_end_2_sam]);

print "Finding rrna anchored concordant read ids\n";
my $rrna_end_1_readids = get_local_filename("rrna.end.1.readids");
my $rrna_end_2_readids = get_local_filename("rrna.end.2.readids");
$runner->run("$filter_sam_mapped_script < #<1 | $sam_readids_script > #>1", [$rrna_end_1_sam], [$rrna_end_1_readids]);
$runner->run("$filter_sam_mapped_script < #<1 | $sam_readids_script > #>1", [$rrna_end_2_sam], [$rrna_end_2_readids]);
push @concordant_readids, $rrna_end_1_readids;
push @concordant_readids, $rrna_end_2_readids;
	
my $prefilter_out_num = 1;
foreach my $prefilter_fasta (@prefilter_fastas)
{
	print "Finding concordant alignments to prefilter fasta $prefilter_fasta\n";
	my $prefilter_pair_sam = get_local_filename("prefilter.".$prefilter_out_num.".pair.sam");
	my $prefilter_pair_readids = get_local_filename("prefilter.".$prefilter_out_num.".pair.readids");
	$runner->run("$align_pair_bin $prefilter_fasta -1 #<1 -2 #<2 #>1", [$reads_end_1_fastq, $reads_end_2_fastq], [$prefilter_pair_sam]);
	$runner->run("$filter_sam_concordant_script < #<1 | $sam_readids_script > #>1", [$prefilter_pair_sam], [$prefilter_pair_readids]);
	
	push @concordant_readids, $prefilter_pair_readids;
	$prefilter_out_num++;
	
	if (lc($remove_job_temp_files) eq "yes")
	{
		unlink $prefilter_pair_sam;
	}
}

print "Excluding IG rearrangements\n";
my $cdna_end_1_ig_readids = get_local_filename("cdna.ig.end.1.readids");
my $cdna_end_2_ig_readids = get_local_filename("cdna.ig.end.2.readids");
my $ig_readids = get_local_filename("ig.readids");
$runner->run("$filter_sam_mapped_script < #<1 | $filter_sam_genes_script $ig_gene_list | $sam_readids_script > #>1", [$cdna_end_1_sam], [$cdna_end_1_ig_readids]);
$runner->run("$filter_sam_mapped_script < #<1 | $filter_sam_genes_script $ig_gene_list | $sam_readids_script > #>1", [$cdna_end_2_sam], [$cdna_end_2_ig_readids]);
$runner->run("$intersect_script #<1 #<2 > #>1", [$cdna_end_1_ig_readids, $cdna_end_2_ig_readids], [$ig_readids]);
push @concordant_readids, $ig_readids;

print "Coallating concordant read ids\n";
my $all_concordant_readids = get_local_filename("concordant.readids");
$runner->run("cat #<A > #>1", [@concordant_readids], [$all_concordant_readids]);

print "Creating list of improper reads\n";
my $dna_cdna_end_1_sam = get_local_filename("dna.cdna.end.1.sam");
my $dna_cdna_end_2_sam = get_local_filename("dna.cdna.end.2.sam");
$runner->run("$match_paired_alignments_script #<1 #<2 > #>1", [$dna_end_1_sam, $cdna_end_1_sam], [$dna_cdna_end_1_sam]);
$runner->run("$match_paired_alignments_script #<1 #<2 > #>1", [$dna_end_2_sam, $cdna_end_2_sam], [$dna_cdna_end_2_sam]);
$runner->run("$match_paired_alignments_script #<1 #<2 | $filter_sam_mapped_script | $filter_sam_readids_script -i #<3 > #>1", [$dna_cdna_end_1_sam, $dna_cdna_end_2_sam, $all_concordant_readids], [$improper_sam]);

print "Dividing sam output files\n";
my $spanning_align_prefix = $job_prefix.".spanning/";
mkdir $spanning_align_prefix if not -d $spanning_align_prefix;
$runner->run("$filter_unmatched_script < #<1 | $divide_sam_chr_pairs_script -t $cdna_regions -p $spanning_align_prefix > #>1", [$improper_sam], [$spanning_filelist]);

print "Finished Alignments\n";

