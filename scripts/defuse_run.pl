#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use File::Spec;
use POSIX qw(strftime);

use FindBin;
use lib "$FindBin::RealBin";
use configdata;
use cmdrunner;
use parsers;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Run the deFuse pipeline for fusion discovery.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -r, --res       Main results filename (default: results.tsv in Output Directory)\n";
push @usage, "  -a, --rescla    Results with a probability column filename (default: results.classify.tsv in Output Directory)\n";
push @usage, "  -b, --resfil    Filtered by the probability threshold results filename (default: results.filtered.tsv in Output Directory)\n";
push @usage, "  -1, --1fastq    Fastq filename 1\n";
push @usage, "  -2, --2fastq    Fastq filename 2\n";
push @usage, "  -n, --name      Library Name (default: Output Directory Suffix)\n";
push @usage, "  -l, --local     Job Local Directory (default: Output Directory)\n";
push @usage, "  -s, --submit    Submitter Type (default: direct)\n";
push @usage, "  -p, --parallel  Maximum Number of Parallel Jobs (default: 1)\n";

my $help;
my $config_filename;
my $dataset_directory;
my $fastq1_filename;
my $fastq2_filename;
my $output_directory;
my $results_filename;
my $results_classify_filename;
my $results_filtered_filename;
my $library_name;
my $joblocal_directory;
my $submitter_type;
my $max_parallel;
my $timestamp;
my $datestring;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'1fastq=s'    => \$fastq1_filename,
	'2fastq=s'    => \$fastq2_filename,
	'output=s'    => \$output_directory,
	'res=s'       => \$results_filename,
	'rescla=s'    => \$results_classify_filename,
	'resfil=s'    => \$results_filtered_filename,
	'name=s'      => \$library_name,
	'local=s'     => \$joblocal_directory,
	'submit=s'    => \$submitter_type,
	'parallel=s'  => \$max_parallel,
);

not defined $help or die @usage;

defined $dataset_directory or die @usage;
defined $output_directory or die @usage;

$timestamp = time();
my $timestamp0 = $timestamp;
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Starting deFuse v0.7.0 for files ".basename($fastq1_filename)." and ".basename($fastq2_filename)."\n";

mkdir $output_directory if not -d $output_directory;

my $source_directory = abs_path("$FindBin::RealBin/../");

if (not defined $config_filename)
{
	$config_filename = $source_directory."/scripts/config.txt";
}

-e $config_filename or die "Error: Unable to find config file $config_filename\n";

$output_directory = abs_path($output_directory);
$config_filename = abs_path($config_filename);

# Main results filename
if (not defined $results_filename)
{
	$results_filename = $output_directory."/results.tsv";
}

# Results with a probability column filename
if (not defined $results_classify_filename)
{
	$results_classify_filename = $output_directory."/results.classify.tsv";
}

# Filtered according to the probability threshold results filename
if (not defined $results_filtered_filename)
{
	$results_filtered_filename = $output_directory."/results.filtered.tsv";
}

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
$config->read($config_filename, $dataset_directory, $source_directory);

# Config values
my $gene_models           = $config->get_value("gene_models");
my $genome_fasta          = $config->get_value("genome_fasta");
my $chromosome_prefix     = $config->get_value("chromosome_prefix");
my $reference_fasta       = $config->get_value("reference_fasta");
my $cdna_fasta            = $config->get_value("cdna_fasta");
my $cdna_regions          = $config->get_value("cdna_regions");
my $exons_fasta           = $config->get_value("exons_fasta");
my $cds_fasta             = $config->get_value("cds_fasta");
my $est_fasta             = $config->get_value("est_fasta");
my @est_split_fastas      = $config->get_list("est_split_fasta");
my $ig_gene_list          = $config->get_value("ig_gene_list");
my $max_insert_size       = $config->get_value("max_insert_size");
my $span_count_threshold  = $config->get_value("span_count_threshold");
my $probability_threshold = $config->get_value("probability_threshold");
my $clustering_precision  = $config->get_value("clustering_precision");
my $scripts_directory     = $config->get_value("scripts_directory");
my $tools_directory       = $config->get_value("tools_directory");
my $reads_per_job         = $config->get_value("reads_per_job");
my $samtools_bin          = $config->get_value("samtools_bin");
my $rscript_bin           = $config->get_value("rscript_bin");
my $split_min_anchor      = $config->get_value("split_min_anchor");
my $cov_samp_density      = $config->get_value("covariance_sampling_density");
my $remove_job_files      = $config->get_value("remove_job_files");
my $positive_controls     = $config->get_value("positive_controls");
my $discord_read_trim     = $config->get_value("discord_read_trim");
my %chromosomes           = $config->get_hash("chromosomes");
my $mt_chromosome         = $config->get_value("mt_chromosome");
my $num_blat_sequences    = $config->get_value("num_blat_sequences");
my $blat_bin              = $config->get_value("blat_bin");
my $dna_concordant_len    = $config->get_value("dna_concordant_length");
my $gmap_bin              = $config->get_value("gmap_bin");
my $gmap_index_directory  = $config->get_value("gmap_index_directory");
my $qsub_params           = $config->get_value("qsub_params");

my $cdna_fasta_fai        = $cdna_fasta.".fai";

my $mailto;
if ($config->has_value("mailto"))
{
	$mailto = $config->get_value("mailto");
}

my $status = "failure";
sub mailme
{
	return if not defined $mailto or not $mailto;

	my $text = "Fusion analysis of library $library_name finished with status $status";

	my $datestring = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print "[$datestring]  Attempting to mail $mailto the result\n";

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

# Ensure required directories exist
verify_directory_exists($scripts_directory);

# Script and binary paths
my $index_paired_fastq_script = "$scripts_directory/index_paired_fastq.pl";
my $alignjob_script = "$scripts_directory/alignjob.pl";
my $read_stats_script = "$scripts_directory/read_stats.pl";
my $expression_script = "$scripts_directory/calculate_expression_simple.pl";
my $split_fastq_script = "$scripts_directory/split_fastq.pl";
my $get_align_regions_script = "$scripts_directory/get_align_regions.pl";
my $merge_read_stats_script = "$scripts_directory/merge_read_stats.pl";
my $merge_cov_samples_script = "$scripts_directory/merge_cov_samples.pl";
my $merge_expression_script = "$scripts_directory/merge_expression.pl";
my $merge_clusters_script = "$scripts_directory/merge_clusters.pl";
my $clustermatepairs_bin = "$tools_directory/clustermatepairs";
my $segregate_mitochondrial_script = "$scripts_directory/segregate_mitochondrial.pl";
my $select_fusion_clusters_script = "$scripts_directory/select_fusion_clusters.pl";
my $prep_local_alignment_seqs_script = "$scripts_directory/prep_local_alignment_seqs.pl";
my $filter_column_script = "$scripts_directory/filter_column.pl";
my $remove_duplicates_script = "$scripts_directory/remove_duplicates.pl";
my $setcover_bin = "$tools_directory/setcover";
my $localalign_bin = "$tools_directory/localalign";
my $dosplitalign_bin = "$tools_directory/dosplitalign";
my $evalsplitalign_bin = "$tools_directory/evalsplitalign";
my $calc_span_stats_script = "$scripts_directory/calc_span_stats.pl";
my $create_read_fasta_script = "$scripts_directory/create_read_fasta.pl";
my $make_fasta_script = "$scripts_directory/make_fasta.pl";
my $split_fasta_script = "$scripts_directory/split_fasta.pl";
my $annotate_fusions_script = "$scripts_directory/annotate_fusions.pl";
my $filter_script = "$scripts_directory/filter.pl";
my $coallate_fusions_script = "$scripts_directory/coallate_fusions.pl";
my $evaluate_fraglength_rscript = "$rscript_bin ".$scripts_directory."/evaluate_fraglength_mean.R";
my $evaluate_split_rscript = "$rscript_bin ".$scripts_directory."/evaluate_split.R";
my $adaboost_rscript = "$rscript_bin ".$scripts_directory."/run_adaboost.R";

mkdir $output_directory if not -e $output_directory;

my $job_directory = $output_directory."/jobs";
my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/defuse";

mkdir $log_directory if not -e $log_directory;
mkdir $job_directory if not -e $job_directory;

my $runner = cmdrunner->new();
$runner->name("defuse");
$runner->prefix($log_prefix);
$runner->maxparallel($max_parallel);
$runner->submitter($submitter_type);
$runner->qsub_params($qsub_params);
$runner->jobmem(6000000000);

# Fastq files and prefix for index and name map
my $reads_prefix = $output_directory."/reads";
my $reads_end_1_fastq = $reads_prefix.".1.fastq";
my $reads_end_2_fastq = $reads_prefix.".2.fastq";
my $reads_index_filename = $reads_prefix.".fqi";
my $reads_names_filename = $reads_prefix.".names";
my $reads_sources_filename = $reads_prefix.".sources";

# Optionally retrieve fastq files
if (defined $fastq1_filename and defined $fastq2_filename)
{
	$timestamp = time();
	$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
	print "[$datestring]  Importing fastq files";
	-e $fastq1_filename or die "Error: Unable to find fastq $fastq1_filename\n";
	-e $fastq2_filename or die "Error: Unable to find fastq $fastq2_filename\n";
	$fastq1_filename = abs_path($fastq1_filename);
	$fastq2_filename = abs_path($fastq2_filename);
	$runner->run("$index_paired_fastq_script #<1 #<2 #>1 #>2 #>3 #>4", [$fastq1_filename,$fastq2_filename], [$reads_end_1_fastq,$reads_end_2_fastq,$reads_index_filename,$reads_names_filename]);
}

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Splitting fastq files";
my $reads_split_prefix = $job_directory."/reads";
my $reads_split_catalog = $output_directory."/reads.split.catalog";
$runner->run("$split_fastq_script $reads_prefix $reads_per_job $reads_split_prefix > #>1", [$reads_end_1_fastq,$reads_end_2_fastq], [$reads_split_catalog]);

# Read split catalog
my @split_fastq_prefixes = readsplitcatalog($reads_split_catalog);

# Create job info
my @job_infos;
foreach my $split_fastq_prefix (@split_fastq_prefixes)
{
	my %job_info;
	$job_info{prefix} = $split_fastq_prefix;
	$job_info{name} = basename($split_fastq_prefix);
	$job_info{fastq1} = $split_fastq_prefix.".1.fastq";
	$job_info{fastq2} = $split_fastq_prefix.".2.fastq";
	
	push @job_infos, \%job_info;
}

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Discordant alignments";
my @job_read_stats;
my @job_spanlength_samples;
my @job_splitpos_samples;
my @job_splitmin_samples;
my @job_expression;
my @job_cdna_pair_sam;
foreach my $job_info (@job_infos)
{
	# Names of job products
	$job_info->{read_stats} = $job_info->{prefix}.".concordant.read.stats";
	$job_info->{spanlength_samples} = $job_info->{prefix}.".spanlength.samples";
	$job_info->{splitpos_samples} = $job_info->{prefix}.".splitpos.samples";
	$job_info->{splitmin_samples} = $job_info->{prefix}.".splitmin.samples";
	$job_info->{expression} = $job_info->{prefix}.".expression.txt";
	$job_info->{cdna_pair_sam} = $job_info->{prefix}.".cdna.pair.sam";
	$job_info->{improper_sam} = $job_info->{prefix}.".improper.sam";
	$job_info->{spanning_filelist} = $job_info->{prefix}.".spanning.filelist";
	
	my $job_cmd = $alignjob_script." ";
	$job_cmd .= "-c ".$config_filename." ";
	$job_cmd .= "-d ".$dataset_directory." ";
	$job_cmd .= "-j ".$job_info->{name}." ";
	$job_cmd .= "-l ".$joblocal_directory." ";
	$job_cmd .= "-o ".$output_directory." ";
	$job_cmd .= "-n ".$library_name." ";
	$job_cmd .= "-p ".$job_info->{prefix}." ";
	
	my @job_products;
	push @job_products, $job_info->{read_stats};
	push @job_products, $job_info->{spanlength_samples};
	push @job_products, $job_info->{splitpos_samples};
	push @job_products, $job_info->{splitmin_samples};
	push @job_products, $job_info->{expression};
	push @job_products, $job_info->{cdna_pair_sam};
	push @job_products, $job_info->{improper_sam};
	push @job_products, $job_info->{spanning_filelist};
	
	$runner->padd("$job_cmd", [$job_info->{fastq1}, $job_info->{fastq2}], [@job_products]); 
	
	push @job_read_stats, $job_info->{read_stats};
	push @job_spanlength_samples, $job_info->{spanlength_samples};
	push @job_splitpos_samples, $job_info->{splitpos_samples};
	push @job_splitmin_samples, $job_info->{splitmin_samples};
	push @job_expression, $job_info->{expression};
	push @job_cdna_pair_sam, $job_info->{cdna_pair_sam};
}
$runner->prun();

# Merge read stats and expression
my $read_stats = $output_directory."/concordant.read.stats";
my $spanlength_cov = $output_directory."/spanlength.cov";
my $splitpos_cov = $output_directory."/splitpos.cov";
my $splitmin_cov = $output_directory."/splitmin.cov";
my $expression = $output_directory."/expression.txt";
my $cdna_pair_sam = $output_directory."/cdna.pair.sam";
$runner->padd("$merge_read_stats_script #<A > #>1", [@job_read_stats], [$read_stats]);
$runner->padd("$merge_cov_samples_script #<A > #>1", [@job_spanlength_samples], [$spanlength_cov]);
$runner->padd("$merge_cov_samples_script #<A > #>1", [@job_splitpos_samples], [$splitpos_cov]);
$runner->padd("$merge_cov_samples_script #<A > #>1", [@job_splitmin_samples], [$splitmin_cov]);
$runner->padd("$merge_expression_script #<A > #>1", [@job_expression], [$expression]);
$runner->padd("cat #<A > #>1", [@job_cdna_pair_sam], [$cdna_pair_sam]);
$runner->prun();

# Read in the read stats
my %read_stat_values;
parsers::get_stats($read_stats, \%read_stat_values);

my $read_length_min = $read_stat_values{"readlength_min"};
my $read_length_max = $read_stat_values{"readlength_max"};
my $fragment_mean = $read_stat_values{"fraglength_mean"};
my $fragment_stddev = $read_stat_values{"fraglength_stddev"};
my $fragment_max = int($fragment_mean + 3 * $fragment_stddev);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Read Stats\n";
print "                         Fragment mean $fragment_mean stddev $fragment_stddev\n";
print "                         Read length min $read_length_min max $read_length_max\n";

if ($fragment_mean / $discord_read_trim < 3)
{
	my $suggested_discord_read_trim = int($fragment_mean / 3);
	if ($suggested_discord_read_trim < 25)
	{
		print "WARNING: cdna fragments too short, deFuse may produce poor results\n";
	}
	else
	{
		print "WARNING: discord_read_trim likely set too high, should be at most $suggested_discord_read_trim\n";
	}
}

# Read divided spanning alignment lists
my %spanning_filenames;
foreach my $job_info (@job_infos)
{
	open SFL, $job_info->{spanning_filelist} or die;
	while (<SFL>)
	{
		chomp;
		my ($chr1,$chr2,$filename) = split /\t/;
		push @{$spanning_filenames{$chr1}{$chr2}}, $filename;
	}
	close SFL;
}

$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Generating discordant alignment clusters";
my @chr_cluster_filenames;
foreach my $chr1 (keys %spanning_filenames)
{
	foreach my $chr2 (keys %{$spanning_filenames{$chr1}})
	{
		my $chr_cluster_filename = $output_directory."/clusters.".$chr1."-".$chr2;
		$runner->padd("cat #<A | $clustermatepairs_bin -m $span_count_threshold -p $clustering_precision -u $fragment_mean -s $fragment_stddev -a - -c #>1", [@{$spanning_filenames{$chr1}{$chr2}}], [$chr_cluster_filename]);
		push @chr_cluster_filenames, $chr_cluster_filename;
	}
}
$runner->prun();
my $clusters_all = $output_directory."/clusters.all";
$runner->run("$merge_clusters_script #<A > #>1", [@chr_cluster_filenames], [$clusters_all]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Remove mitochondrial-genomic clusters";
my $clusters = $output_directory."/clusters";
$runner->run("$segregate_mitochondrial_script $gene_models $mt_chromosome < #<1 > #>1", [$clusters_all], [$clusters]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Generating maximum parsimony solution";
my $clusters_sc_all = $output_directory."/clusters.sc.all";
$runner->jobmem(32000000000);
$runner->run("$setcover_bin -m $span_count_threshold -c #<1 -o #>1", [$clusters], [$clusters_sc_all]);
$runner->jobmem(8000000000);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Selecting fusion clusters";
my $clusters_sc_unfilt = $output_directory."/clusters.sc.unfilt";
$runner->run("$select_fusion_clusters_script $gene_models < #<1 > #>1", [$clusters_sc_all], [$clusters_sc_unfilt]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Preparing sequences for local realignment";
my $clusters_sc_local_seq = $output_directory."/clusters.sc.local.seq";
$runner->run("$prep_local_alignment_seqs_script -r $reference_fasta -g $gene_models -c #<1 -s $dna_concordant_len > #>1", [$clusters_sc_unfilt], [$clusters_sc_local_seq]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Performing local realignment";
my $clusters_sc_local_align = $output_directory."/clusters.sc.local.align";
$runner->run("$localalign_bin -m 10 -x -5 -g -5 -t 0.8 < #<1 > #>1", [$clusters_sc_local_seq], [$clusters_sc_local_align]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Filtering concordant clusters";
my $clusters_sc = $output_directory."/clusters.sc";
$runner->run("cat #<1 | $filter_column_script #<2 0 1 | $remove_duplicates_script $span_count_threshold > #>1", [$clusters_sc_unfilt,$clusters_sc_local_align], [$clusters_sc]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Generating spanning alignment regions file";
my $clusters_sc_regions = $output_directory."/clusters.sc.regions";
$runner->run("$get_align_regions_script < #<1 > #>1", [$clusters_sc], [$clusters_sc_regions]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Calculating split read alignments";
foreach my $job_info (@job_infos)
{
	$job_info->{splitalign} = $job_info->{prefix}.".splitreads.alignments";
	$runner->padd("$dosplitalign_bin -r #<1 -n $read_length_min -x $read_length_max -u $fragment_mean -s $fragment_stddev -e $cdna_regions -f $reference_fasta -1 #<2 -2 #<3 -i #<4 -a #>1", [$clusters_sc_regions,$job_info->{fastq1},$job_info->{fastq2},$job_info->{improper_sam}], [$job_info->{splitalign}]);
}
$runner->prun();
my @job_splitreads_alignments;
foreach my $job_info (@job_infos)
{
	$job_info->{splitalignsorted} = $job_info->{prefix}.".splitreads.alignments.sorted";
	$runner->padd("sort -T $joblocal_directory -n -k 1 #<1 > #>1", [$job_info->{splitalign}], [$job_info->{splitalignsorted}]);
	push @job_splitreads_alignments, $job_info->{splitalignsorted};
}
my $splitreads_alignments = $output_directory."/splitreads.alignments";
$runner->prun();
$runner->run("sort -m -n -k 1 #<A > #>1", [@job_splitreads_alignments], [$splitreads_alignments]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Evaluating split reads";
my $splitreads_break = $output_directory."/splitreads.break";
my $splitreads_seq = $output_directory."/splitreads.seq";
my $splitreads_predalign = $output_directory."/splitreads.predalign";
$runner->run("$evalsplitalign_bin -r #<1 -n $read_length_min -x $read_length_max -u $fragment_mean -s $fragment_stddev -e $cdna_regions -f $reference_fasta -a #<2 -b #>1 -q #>2 -p #>3", [$clusters_sc_regions,$splitreads_alignments], [$splitreads_break,$splitreads_seq,$splitreads_predalign]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Calculating spanning stats";
my $splitreads_span_stats = $output_directory."/splitreads.span.stats";
$runner->padd("$calc_span_stats_script -c #<1 -b #<2 -s #<3 > #>1", [$clusters_sc,$splitreads_break,$splitreads_seq], [$splitreads_span_stats]);
$runner->prun();

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Calculating spanning p-values";
my $splitreads_span_pval = $output_directory."/splitreads.span.pval";
$runner->padd("$evaluate_fraglength_rscript #<1 #<2 $discord_read_trim #<3 #>1", [$read_stats,$spanlength_cov,$splitreads_span_stats], [$splitreads_span_pval]);
$runner->prun();

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Calculating split read pvalues";
my $splitreads_split_pval = $output_directory."/splitreads.split.pval";
$runner->run("$evaluate_split_rscript #<1 #<2 #<3 #>1", [$splitpos_cov,$splitmin_cov,$splitreads_seq], [$splitreads_split_pval]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Creating fastas";
my $breakpoints_fasta = $output_directory."/breakpoints.fa";
$runner->run("cut -f1,2 #<1 | $make_fasta_script > #>1", [$splitreads_seq], [$breakpoints_fasta]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Splitting fastas";
my $breakpoints_split_prefix = $job_directory."/breakpoints";
my $breakpoints_split_catalog = $output_directory."/breakpoints.split.catalog";
$runner->run("$split_fasta_script #<1 $num_blat_sequences $breakpoints_split_prefix > #>1", [$breakpoints_fasta], [$breakpoints_split_catalog]);
my @breakpoints_split_fastas = readsplitcatalog($breakpoints_split_catalog);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Breakpoint alignments";
my $breakpoints_genome_psl = $output_directory."/breakpoints.genome.psl";
my $breakpoints_genome_nointron_psl = $output_directory."/breakpoints.genome.nointron.psl";
my $breakpoints_cdna_psl = $output_directory."/breakpoints.cdna.psl";
my $breakpoints_est_psl = $output_directory."/breakpoints.est.psl";
my $breakpoints_exons_psl = $output_directory."/breakpoints.exons.psl";
my $breakpoints_cds_psl = $output_directory."/breakpoints.cds.psl";

my @breakpoint_jobs;
push @breakpoint_jobs, [$breakpoints_fasta, $breakpoints_genome_psl, \@breakpoints_split_fastas, "genome", "", "genome"];
push @breakpoint_jobs, [$breakpoints_fasta, $breakpoints_genome_nointron_psl, \@breakpoints_split_fastas, "genome", "--nosplicing", "genome.nointron"];
push @breakpoint_jobs, [$breakpoints_fasta, $breakpoints_cdna_psl, \@breakpoints_split_fastas, "cdna", "", "cdna"];
push @breakpoint_jobs, [$breakpoints_fasta, $breakpoints_est_psl, \@breakpoints_split_fastas, "est", "", "est"];
push @breakpoint_jobs, [$breakpoints_fasta, $breakpoints_exons_psl, \@breakpoints_split_fastas, "exons", "", "exons"];
push @breakpoint_jobs, [$breakpoints_fasta, $breakpoints_cds_psl, \@breakpoints_split_fastas, "cds", "", "cds"];

do_breakpoint_jobs(@breakpoint_jobs);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Annotating fusions";
$runner->jobmem(16000000000);
my $annotations_filename = $output_directory."/annotations";
$runner->run("$annotate_fusions_script -c $config_filename -d $dataset_directory -o $output_directory -n $library_name > #>1", [$splitreads_span_pval,$splitreads_break,$splitreads_seq,$clusters_sc], [$annotations_filename]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Coallating fusions";
$runner->run("$coallate_fusions_script -c $config_filename -d $dataset_directory -o $output_directory -l #<1 > #>1", [$annotations_filename], [$results_filename]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Running adaboost classifier";
$runner->run("$adaboost_rscript $positive_controls #<1 #>1", [$results_filename], [$results_classify_filename]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Filtering fusions";
$runner->run("$filter_script probability '> $probability_threshold' < #<1 > #>1", [$results_classify_filename], [$results_filtered_filename]);

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  Success";

# Remove job files
if (lc($remove_job_files) eq "yes")
{
	foreach my $split_fastq_prefix (@split_fastq_prefixes)
	{
		unlink $split_fastq_prefix.".1.fastq";
		unlink $split_fastq_prefix.".2.fastq";
		unlink $split_fastq_prefix.".candidate.reads.alignments";
		unlink $split_fastq_prefix.".candidate.reads";
	}
}

$status = "success";

print " [".(time() - $timestamp)." sec]\n";
$timestamp = time();
$datestring = strftime("%Y-%m-%d %H:%M:%S", localtime($timestamp));
print "[$datestring]  deFuse finished. Total time: ".(time() - $timestamp0)." sec.\n";

sub readsplitcatalog
{
	my $split_catalog = shift;
	my @split_prefixes;

	open SC, $split_catalog or die "Error: Unable to open $split_catalog: $!\n";
	while (<SC>)
	{
		chomp;
		my @fields = split /\t/;

		push @split_prefixes, $fields[0];
	}
	close SC;

	return @split_prefixes;
}

sub do_breakpoint_jobs
{
	my @merge_jobs;
	foreach my $breakpoint_job (@_)
	{
		my @breakpoint_job_info = @{$breakpoint_job};
		
		my $fasta = shift @breakpoint_job_info;
		my $psl = shift @breakpoint_job_info;
		if (not cmdrunner::uptodate([$fasta], [$psl]))
		{
			my @split_psls = padd_do_breakpoint_alignments(@breakpoint_job_info);
			push @merge_jobs, [\@split_psls, $psl];
		}
	}
	$runner->prun();
	
	foreach my $merge_job (@merge_jobs)
	{
		if (@{$merge_job->[0]} == 0)
		{
			system "touch $merge_job->[1]";
		}
		else
		{
			$runner->padd("cat #<A > #>1", $merge_job->[0], [$merge_job->[1]]);
		}
	}
	$runner->prun();
}

sub padd_do_breakpoint_alignments
{
	my $split_fastas_ref = shift;
	my $ref_name = shift;
	my $parameters = shift;
	my $align_name = shift;
	
	if ($ref_name eq "genome")
	{
		return padd_do_gmap_alignments_genome($split_fastas_ref, $parameters, $align_name);
	}
	elsif ($ref_name eq "est")
	{
		return padd_do_gmap_alignments_est($split_fastas_ref, $parameters, $align_name);
	}
	elsif ($ref_name eq "cds")
	{
		return padd_do_blat_alignments($split_fastas_ref, $cds_fasta, $parameters, $align_name);
	}
	elsif ($ref_name eq "exons")
	{
		return padd_do_blat_alignments($split_fastas_ref, $exons_fasta, $parameters, $align_name);
	}
	else
	{
		return padd_do_gmap_alignments($split_fastas_ref, $ref_name, $parameters, $align_name);
	}
}
	
sub padd_do_gmap_alignments
{
	my $split_fastas_ref = shift;
	my $ref_name = shift;
	my $parameters = shift;
	my $align_name = shift;
	
	my @output_psls;
	foreach my $split_fasta (@{$split_fastas_ref})
	{
		my $output_psl = $split_fasta.".".$align_name.".psl";
		$runner->padd("$gmap_bin -D $gmap_index_directory -d $ref_name -f psl $parameters #<1 > #>1", [$split_fasta], [$output_psl]);
		push @output_psls, $output_psl;
	}
	
	return @output_psls;
}

sub padd_do_blat_alignments
{
	my $split_fastas_ref = shift;
	my $reference = shift;
	my $parameters = shift;
	my $align_name = shift;
	
	my @output_psls;
	foreach my $split_fasta (@{$split_fastas_ref})
	{
		my $output_psl = $split_fasta.".".$align_name.".psl";
		$runner->padd("$blat_bin -noHead $parameters $reference.2bit #<1 #>1", [$split_fasta], [$output_psl]);
		push @output_psls, $output_psl;
	}
	
	return @output_psls;
}

sub padd_do_gmap_alignments_genome
{
	my $split_fastas_ref = shift;
	my $gmap_parameters = shift;
	my $align_name = shift;
	
	my @output_psls;
	foreach my $chromosome (keys %chromosomes)
	{
		push @output_psls, padd_do_gmap_alignments($split_fastas_ref, "chr".$chromosome, $gmap_parameters, $align_name.".".$chromosome);
	}
	
	return @output_psls;
}

sub padd_do_gmap_alignments_est
{
	my $split_fastas_ref = shift;
	my $gmap_parameters = shift;
	my $align_name = shift;
	
	my @output_psls;
	foreach my $est_split_index (0..$#est_split_fastas)
	{
		push @output_psls, padd_do_gmap_alignments($split_fastas_ref, "est".$est_split_index, $gmap_parameters, $align_name.".".$est_split_index);
	}
	
	return @output_psls;
}



