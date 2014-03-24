#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];

use lib dirname($0);
use configdata;
use cmdrunner;

my $main_pid = $$;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Run the bowtie alignment pipeline for fusions.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -l, --local     Job Local Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --name      Library Name\n";
push @usage, "  -j, --job       Job Directory\n";
push @usage, "  -p, --prefix    Prefix of Input/Output\n";

my $help;
my $config_filename;
my $local_directory;
my $output_directory;
my $library_name;
my $job_directory;
my $prefix;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'local=s'     => \$local_directory,
	'output=s'    => \$output_directory,
	'name=s'      => \$library_name,
	'job=s'       => \$job_directory,
	'prefix=s'    => \$prefix,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
defined $local_directory or die @usage;
defined $output_directory or die @usage;
defined $library_name or die @usage;
defined $job_directory or die @usage;
defined $prefix or die @usage;

print "Starting alignjob $prefix\n";

my $config = configdata->new();
$config->read($config_filename);

# Config values
my $cdna_gene_fasta      = $config->get_value("cdna_gene_fasta_cluster");
my $cdna_fasta           = $config->get_value("cdna_fasta_cluster");
my $rrna_fasta           = $config->get_value("rrna_fasta");
my $ig_gene_list         = $config->get_value("ig_gene_list");
my $max_insert_size      = $config->get_value("max_insert_size");
my $scripts_directory    = $config->get_value("scripts_directory");
my $bowtie_bin           = $config->get_value("bowtie_bin");
my $samtools_bin         = $config->get_value("samtools_bin");
my $bowtie_threads       = $config->get_value("bowtie_threads");
my @prefilter_fastas     = $config->get_list("prefilter");

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

# Samtools fasta indices required
my $cdna_fasta_index = $cdna_fasta.".fai";
my $cdna_gene_fasta_index = $cdna_gene_fasta.".fai";

# Ensure required files exist
verify_file_exists($cdna_fasta);
verify_file_exists($cdna_fasta_index);
verify_file_exists($cdna_gene_fasta);
verify_file_exists($cdna_gene_fasta_index);
verify_file_exists($ig_gene_list);

# Ensure required directories exist
verify_directory_exists($scripts_directory);
verify_directory_exists($local_directory);
verify_directory_exists($output_directory);

my $align_pair_bin = "$bowtie_bin --sam-nosq -S --mm -t -X $max_insert_size";
my $align_single_bin = "$bowtie_bin --sam-nosq -S --mm -t -k 100 -m 100";
my $sort_cmd = "sort";

my $filter_sam_readids_script = "$scripts_directory/filter_sam_readids.pl";
my $filter_sam_genes_script = "$scripts_directory/filter_sam_reference.pl";
my $read_stats_script = "$scripts_directory/read_stats.pl";
my $expression_script = "$scripts_directory/calculate_expression_simple.pl";
my $find_concordant_ensembl_script = "$scripts_directory/find_concordant_ensembl.pl";
my $sam_readids_script = "$scripts_directory/sam_readids.pl";
my $filter_sam_concordant_script = "$scripts_directory/filter_sam_concordant.pl";
my $filter_sam_mapped_script = "$scripts_directory/filter_sam_mapped.pl";

my $job_prefix = $job_directory."/".$prefix;
my $local_prefix = abs_path($local_directory)."/".$library_name.".".$prefix;

# List of local files to remove on control-c
my @local_filenames;

# Cleanup method to remove files
sub cleanup
{
	return if scalar @local_filenames == 0;
	print "Cleaning Up\n";
	foreach my $filename (@local_filenames)
	{
		unlink $filename if -e $filename;
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
my $log_prefix = $log_directory."/".$prefix;

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("alignjob.$prefix");
$runner->prefix($log_prefix);
$runner->submitter("direct");

my $reads_end_1_fastq = $job_prefix.".1.fastq";
my $reads_end_2_fastq = $job_prefix.".2.fastq";

-e $reads_end_1_fastq or die "Error: Unable to find fastq file $reads_end_1_fastq\n";
-e $reads_end_2_fastq or die "Error: Unable to find fastq file $reads_end_2_fastq\n";

# Products of the alignment process
my $read_stats = $job_prefix.".concordant.read.stats";
my $expression = $job_prefix.".expression.txt";
my $cdna_pair_bam = $job_prefix.".cdna.pair.bam";
my $discordant_aligned_bam = $job_prefix.".discordant.aligned.bam";
my $discordant_unaligned_bam = $job_prefix.".discordant.unaligned.bam";

print "Finding concordant alignments to cdna\n";
my $cdna_pair_sam = get_local_filename("cdna.pair.sam");
$runner->run("$align_pair_bin $cdna_fasta -1 #<1 -2 #<2 #>1", [$reads_end_1_fastq, $reads_end_2_fastq], [$cdna_pair_sam]);

print "Calculating read statistics\n";
$runner->run("$read_stats_script < #<1 > #>1", [$cdna_pair_sam], [$read_stats]);

print "Calculating expression\n";
$runner->run("$expression_script < #<1 > #>1", [$cdna_pair_sam], [$expression]);

if (not $runner->uptodate([$reads_end_1_fastq, $reads_end_2_fastq], [$discordant_aligned_bam, $discordant_unaligned_bam]))
{
	print "Finding alignable reads\n";
	my $cdna_end_1_sam = get_local_filename("cdna.end.1.sam");
	my $cdna_end_2_sam = get_local_filename("cdna.end.2.sam");
	$runner->run("$align_single_bin $cdna_gene_fasta #<1 #>1", [$reads_end_1_fastq], [$cdna_end_1_sam]);
	$runner->run("$align_single_bin $cdna_gene_fasta #<1 #>1", [$reads_end_2_fastq], [$cdna_end_2_sam]);
	
	my $concordant_readids = get_local_filename("concordant.readids");

	if (not -e $concordant_readids)
	{
		print "Finding same gene concordant read ids for cdna\n";
		my $cdna_concordant_readids = get_local_filename("cdna.concordant.readids");
		$runner->run("$find_concordant_ensembl_script #<1 #<2 | $sort_cmd > #>1", [$cdna_end_1_sam, $cdna_end_2_sam], [$cdna_concordant_readids]);

		print "Finding concordant read ids for cdna\n";	
		my $cdna_pair_readids = get_local_filename("cdna.pair.readids");
		$runner->run("$filter_sam_concordant_script < #<1 | $sam_readids_script | $sort_cmd -u > #>1", [$cdna_pair_sam], [$cdna_pair_readids]);

		print "Finding alignments to rRNA\n";
		my $rrna_end_1_sam = get_local_filename("rrna.end.1.sam");
		my $rrna_end_2_sam = get_local_filename("rrna.end.2.sam");
		$runner->run("$align_single_bin $rrna_fasta #<1 #>1", [$reads_end_1_fastq], [$rrna_end_1_sam]);
		$runner->run("$align_single_bin $rrna_fasta #<1 #>1", [$reads_end_2_fastq], [$rrna_end_2_sam]);

		print "Finding rrna anchored concordant read ids\n";
		my $rrna_end_1_readids = get_local_filename("rrna.end.1.readids");
		my $rrna_end_2_readids = get_local_filename("rrna.end.2.readids");
		$runner->run("$filter_sam_mapped_script < #<1 | $sam_readids_script | $sort_cmd -u > #>1", [$rrna_end_1_sam], [$rrna_end_1_readids]);
		$runner->run("$filter_sam_mapped_script < #<1 | $sam_readids_script | $sort_cmd -u > #>1", [$rrna_end_2_sam], [$rrna_end_2_readids]);
		
		my @all_prefilter_pair_readids;
		my $prefilter_out_num = 1;
		foreach my $prefilter_fasta (@prefilter_fastas)
		{
			print "Finding concordant alignments to prefilter fasta $prefilter_fasta\n";
			my $prefilter_pair_sam = get_local_filename("prefilter.".$prefilter_out_num.".pair.sam");
			my $prefilter_pair_readids = get_local_filename("prefilter.".$prefilter_out_num.".pair.readids");
			$runner->run("$align_pair_bin $prefilter_fasta -1 #<1 -2 #<2 #>1", [$reads_end_1_fastq, $reads_end_2_fastq], [$prefilter_pair_sam]);
			$runner->run("$filter_sam_concordant_script < #<1 | $sam_readids_script | $sort_cmd -u > #>1", [$prefilter_pair_sam], [$prefilter_pair_readids]);
			
			push @all_prefilter_pair_readids, $prefilter_pair_readids;
			unlink $prefilter_pair_sam;
			$prefilter_out_num++;
		}

		print "Excluding IG rearrangements\n";
		my $cdna_end_1_ig_readids = get_local_filename("cdna.ig.end.1.readids");
		my $cdna_end_2_ig_readids = get_local_filename("cdna.ig.end.2.readids");
		my $ig_readids = get_local_filename("ig.readids");
		$runner->run("$filter_sam_mapped_script < #<1 | $filter_sam_genes_script $ig_gene_list | $sam_readids_script | uniq | $sort_cmd -u > #>1", [$cdna_end_1_sam], [$cdna_end_1_ig_readids]);
		$runner->run("$filter_sam_mapped_script < #<1 | $filter_sam_genes_script $ig_gene_list | $sam_readids_script | uniq | $sort_cmd -u > #>1", [$cdna_end_2_sam], [$cdna_end_2_ig_readids]);
		$runner->run("join #<1 #<2 > #>1", [$cdna_end_1_ig_readids, $cdna_end_2_ig_readids], [$ig_readids]);
	
		print "Coallating concordant read ids\n";
		$runner->run("$sort_cmd -u -m #<A > #>1", [$cdna_concordant_readids, $rrna_end_1_readids, $rrna_end_2_readids, $cdna_pair_readids, $ig_readids, @all_prefilter_pair_readids], [$concordant_readids]);
		
		unlink $cdna_pair_readids;
		unlink @all_prefilter_pair_readids;
		unlink $cdna_end_1_ig_readids;
		unlink $cdna_end_2_ig_readids;
		unlink $ig_readids;
	}
	
	print "Extracting aligned read ids\n";
	my $cdna_end_1_readids = get_local_filename("cdna.end.1.readids");
	my $cdna_end_2_readids = get_local_filename("cdna.end.2.readids");
	my $cdna_all_readids = get_local_filename("cdna.all.readids");
	my $cdna_aligned_readids = get_local_filename("cdna.aligned.readids");
	my $cdna_unaligned_readids = get_local_filename("cdna.unaligned.readids");
	my $discordant_aligned_readids = get_local_filename("discordant.aligned.readids");
	my $discordant_unaligned_readids = get_local_filename("discordant.unaligned.readids");	
	$runner->run("$filter_sam_mapped_script < #<1 | $sam_readids_script | uniq | $sort_cmd -u > #>1", [$cdna_end_1_sam], [$cdna_end_1_readids]);
	$runner->run("$filter_sam_mapped_script < #<1 | $sam_readids_script | uniq | $sort_cmd -u > #>1", [$cdna_end_2_sam], [$cdna_end_2_readids]);
	$runner->run("$sort_cmd -m #<1 #<2 | uniq > #>1", [$cdna_end_1_readids, $cdna_end_2_readids], [$cdna_all_readids]);
	$runner->run("join #<1 #<2 > #>1", [$cdna_end_1_readids, $cdna_end_2_readids], [$cdna_aligned_readids]);
	$runner->run("join -v 1 -v 2 #<1 #<2 > #>1", [$cdna_end_1_readids, $cdna_end_2_readids], [$cdna_unaligned_readids]);
	$runner->run("join -v 1 #<1 #<2 > #>1", [$cdna_aligned_readids, $concordant_readids], [$discordant_aligned_readids]);
	$runner->run("join -v 1 #<1 #<2 > #>1", [$cdna_unaligned_readids, $concordant_readids], [$discordant_unaligned_readids]);
	
	print "Dividing sam output files\n";
	my $discordant_aligned_sam = get_local_filename("discordant.aligned.sam");
	my $discordant_unaligned_sam = get_local_filename("discordant.unaligned.sam");
	$runner->run("cat #<1 #<2 | $filter_sam_readids_script #<3 | $sort_cmd > #>1", [$cdna_end_1_sam, $cdna_end_2_sam, $discordant_aligned_readids], [$discordant_aligned_sam]);
	$runner->run("cat #<1 #<2 | $filter_sam_readids_script #<3 | $sort_cmd > #>1", [$cdna_end_1_sam, $cdna_end_2_sam, $discordant_unaligned_readids], [$discordant_unaligned_sam]);
	
	print "Converting sam to bam\n";
	my $cdna_pair_bam_prefix = $cdna_pair_bam;
	my $discordant_aligned_bam_prefix = $discordant_aligned_bam;
	my $discordant_unaligned_bam_prefix = $discordant_unaligned_bam;
	$cdna_pair_bam_prefix =~ s/.bam$//;
	$discordant_aligned_bam_prefix =~ s/.bam$//;
	$discordant_unaligned_bam_prefix =~ s/.bam$//;
	$runner->run("$samtools_bin view -bt $cdna_fasta_index #<1 | $samtools_bin sort - $cdna_pair_bam_prefix", [$cdna_pair_sam], [$cdna_pair_bam]);
	$runner->run("$samtools_bin view -bt $cdna_gene_fasta_index #<1 | $samtools_bin sort - $discordant_aligned_bam_prefix", [$discordant_aligned_sam], [$discordant_aligned_bam]);
	$runner->run("$samtools_bin view -bt $cdna_gene_fasta_index #<1 | $samtools_bin sort - $discordant_unaligned_bam_prefix", [$discordant_unaligned_sam], [$discordant_unaligned_bam]);

	unlink $cdna_end_1_sam;
	unlink $cdna_end_2_sam;
	unlink $cdna_end_1_readids;
	unlink $cdna_end_2_readids;
	unlink $cdna_all_readids;
	unlink $cdna_aligned_readids;
	unlink $cdna_unaligned_readids;
	unlink $discordant_aligned_readids;
	unlink $discordant_unaligned_readids;
	unlink $concordant_readids;
}

print "Finished Alignments\n";

