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
push @usage, "  -d, --data      Source Data Directory (default: Skip data retrieval)\n";
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
my $gene_models           = $config->get_value("gene_models");
my $genome_fasta          = $config->get_value("genome_fasta");
my $chromosome_prefix     = $config->get_value("chromosome_prefix");
my $reference_fasta       = $config->get_value("reference_fasta");
my $cdna_fasta            = $config->get_value("cdna_fasta");
my $cdna_regions          = $config->get_value("cdna_regions");
my $exons_fasta           = $config->get_value("exons_fasta");
my $cds_fasta             = $config->get_value("cds_fasta");
my $rrna_fasta            = $config->get_value("rrna_fasta");
my $est_fasta             = $config->get_value("est_fasta");
my $est_split_catalog     = $config->get_value("est_split_catalog");
my $ig_gene_list          = $config->get_value("ig_gene_list");
my $max_insert_size       = $config->get_value("max_insert_size");
my $span_count_threshold  = $config->get_value("span_count_threshold");
my $probability_threshold = $config->get_value("probability_threshold");
my $clustering_precision  = $config->get_value("clustering_precision");
my $scripts_directory     = $config->get_value("scripts_directory");
my $tools_directory       = $config->get_value("tools_directory");
my $reads_per_job         = $config->get_value("reads_per_job");
my $regions_per_job       = $config->get_value("regions_per_job");
my $samtools_bin          = $config->get_value("samtools_bin");
my $rscript_bin           = $config->get_value("rscript_bin");
my $split_min_anchor      = $config->get_value("split_min_anchor");
my $cov_samp_density      = $config->get_value("covariance_sampling_density");
my $denovo_assembly       = $config->get_value("denovo_assembly");
my $remove_job_files      = $config->get_value("remove_job_files");
my $positive_controls     = $config->get_value("positive_controls");
my $discord_read_trim     = $config->get_value("discord_read_trim");
my %chromosomes           = $config->get_hash("chromosomes");
my $mt_chromosome         = $config->get_value("mt_chromosome");
my $num_blat_sequences    = $config->get_value("num_blat_sequences");
my $blat_bin              = $config->get_value("blat_bin");
my $percident_threshold   = $config->get_value("percent_identity_threshold");
my $dna_concordant_len    = $config->get_value("dna_concordant_length");

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
verify_file_exists($reference_fasta);
verify_file_exists($cdna_fasta);
verify_file_exists($cdna_fasta_fai);
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
my $clustermatepairs_bin = "$tools_directory/clustermatepairs";
my $segregate_mitochondrial_script = "$scripts_directory/segregate_mitochondrial.pl";
my $select_fusion_clusters_script = "$scripts_directory/select_fusion_clusters.pl";
my $prep_local_alignment_seqs_script = "$scripts_directory/prep_local_alignment_seqs.pl";
my $get_fusion_sequences_script = "$scripts_directory/get_fusion_sequences.pl";
my $filter_psl_concordant_script = "$scripts_directory/filter_psl_concordant.pl";
my $filter_column_script = "$scripts_directory/filter_column.pl";
my $select_breakpoint_seq_script = "$scripts_directory/select_breakpoint_seq.pl";
my $remove_duplicates_script = "$scripts_directory/remove_duplicates.pl";
my $setcover_bin = "$tools_directory/setcover";
my $localalign_bin = "$tools_directory/localalign";
my $split_seq_bin = "$tools_directory/splitseq";
my $denovo_seq_bin = "$tools_directory/denovoseq";
my $findcandidatereads_bin = "$tools_directory/findcandidatereads";
my $dosplitalign_bin = "$tools_directory/dosplitalign";
my $evalsplitalign_bin = "$tools_directory/evalsplitalign";
my $calc_span_stats_script = "$scripts_directory/calc_span_stats.pl";
my $create_read_fasta_script = "$scripts_directory/create_read_fasta.pl";
my $make_fasta_script = "$scripts_directory/make_fasta.pl";
my $split_fasta_script = "$scripts_directory/split_fasta.pl";
my $annotate_fusions_script = "$scripts_directory/annotate_fusions.pl";
my $calc_map_stats_script = "$scripts_directory/calculate_mapping_stats.pl";
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
$runner->jobmem(3000000000);

# Fastq files and prefix for index and name map
my $reads_prefix = $output_directory."/reads";
my $reads_end_1_fastq = $reads_prefix.".1.fastq";
my $reads_end_2_fastq = $reads_prefix.".2.fastq";
my $reads_index_filename = $reads_prefix.".fqi";
my $reads_names_filename = $reads_prefix.".names";
my $reads_sources_filename = $reads_prefix.".sources";

# Optionally retrieve fastq files
if (defined $source_directory)
{
	print "Retreiving fastq files\n";
	-d $source_directory or die "Error: Unable to find source directory $source_directory\n";
	$source_directory = abs_path($source_directory);
	$runner->run("$retreive_fastq_script -c $config_filename -d $source_directory -1 $reads_end_1_fastq -2 $reads_end_2_fastq -i $reads_index_filename -n $reads_names_filename -s $reads_sources_filename -f", [], []);
}

print "Splitting fastq files\n";
my $reads_split_prefix = $job_directory."/reads";
my $reads_split_catalog = $output_directory."/reads.split.catalog";
$runner->run("$split_fastq_script $reads_prefix $reads_per_job $reads_split_prefix > #>1", [$reads_end_1_fastq,$reads_end_2_fastq], [$reads_split_catalog]);

# Read split catalog
my @split_fastq_prefixes = readsplitcatalog($reads_split_catalog);

# Products of the alignment process
my $read_stats = $output_directory."/concordant.read.stats";
my $expression = $output_directory."/expression.txt";
my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
my $spanning_bam = $output_directory."/spanning.bam";
my $spanning_byread_bam = $output_directory."/spanning.byread.bam";
my $anchored_bam = $output_directory."/anchored.bam";

print "Split align merge for discordant alignments\n";
if (not cmdrunner::uptodate([$reads_end_1_fastq,$reads_end_2_fastq], [$read_stats,$cdna_pair_bam,$spanning_bam,$spanning_byread_bam,$anchored_bam,$expression]))
{
	splitalignmerge(@split_fastq_prefixes);
}

print "Indexing discordant bam files\n";
my $cdna_pair_bam_bai = $cdna_pair_bam.".bai";
my $spanning_bam_bai = $spanning_bam.".bai";
my $anchored_bam_bai = $anchored_bam.".bai";
$runner->padd("$samtools_bin index #<1 #>1", [$cdna_pair_bam], [$cdna_pair_bam_bai]);
$runner->padd("$samtools_bin index #<1 #>1", [$spanning_bam], [$spanning_bam_bai]);
$runner->padd("$samtools_bin index #<1 #>1", [$anchored_bam], [$anchored_bam_bai]);
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

# Run remaining commands on high memory cluster nodes
$runner->jobmem(10000000000);

print "Finding paired end alignment covariance\n";
my $cov_stats = $output_directory."/covariance.stats";
$runner->run("$calccov_bin -m $fragment_max -a $split_min_anchor -d $cov_samp_density -g $cdna_regions -c #<1 -v #>1", [$cdna_pair_bam], [$cov_stats]);

print "Generating maximal valid clusters\n";
my $clusters_all = $output_directory."/clusters.all";
$runner->run("$clustermatepairs_bin -m $span_count_threshold -p $clustering_precision -u $fragment_mean -s $fragment_stddev -r #<1 -c #>1", [$spanning_byread_bam], [$clusters_all]);

print "Remove mitochondrial-genomic clusters\n";
my $clusters = $output_directory."/clusters";
$runner->run("$segregate_mitochondrial_script $gene_models $mt_chromosome < #<1 > #>1", [$clusters_all], [$clusters]);

print "Generating maximum parsimony solution\n";
my $clusters_sc_all = $output_directory."/clusters.sc.all";
$runner->run("$setcover_bin -m $span_count_threshold -c #<1 -o #>1", [$clusters], [$clusters_sc_all]);

print "Selecting fusion clusters\n";
my $clusters_sc_unfilt = $output_directory."/clusters.sc.unfilt";
$runner->run("$select_fusion_clusters_script $gene_models < #<1 > #>1", [$clusters_sc_all], [$clusters_sc_unfilt]);

print "Preparing sequences for local realignment\n";
my $clusters_sc_local_seq = $output_directory."/clusters.sc.local.seq";
$runner->run("$prep_local_alignment_seqs_script -r $reference_fasta -g $gene_models -c #<1 -s $dna_concordant_len > #>1", [$clusters_sc_unfilt], [$clusters_sc_local_seq]);

print "Performing local realignment\n";
my $clusters_sc_local_align = $output_directory."/clusters.sc.local.align";
$runner->run("$localalign_bin -m 10 -x -5 -g -5 -t 0.8 < #<1 > #>1", [$clusters_sc_local_seq], [$clusters_sc_local_align]);

print "Filtering concordant clusters\n";
my $clusters_sc = $output_directory."/clusters.sc";
$runner->run("cat #<1 | $filter_column_script #<2 0 1 | $remove_duplicates_script $span_count_threshold > #>1", [$clusters_sc_unfilt,$clusters_sc_local_align], [$clusters_sc]);

print "Generating spanning alignment regions file\n";
my $clusters_sc_regions = $output_directory."/clusters.sc.regions";
$runner->run("$get_align_regions_script < #<1 > #>1", [$clusters_sc], [$clusters_sc_regions]);

print "Finding candidate split reads\n";
my $candidate_sequences = $output_directory."/candidate.seqs";
my $candidate_reads = $output_directory."/candidate.reads";
$runner->run("$findcandidatereads_bin -i #<1 -c #>1 -r #>2 -m $read_length_min -x $read_length_max -u $fragment_mean -s $fragment_stddev -e $cdna_regions -f $reference_fasta -d $spanning_bam -a $anchored_bam", [$clusters_sc_regions], [$candidate_sequences,$candidate_reads]);

print "Calculating split reads\n";
my $candidate_alignments = $output_directory."/candidate.alignments";
if (not cmdrunner::uptodate([$candidate_sequences,$candidate_reads], [$candidate_alignments]))
{
	dosplitalignments(\@split_fastq_prefixes, $candidate_alignments);
}

print "Evaluating split reads\n";
my $splitr_break = $output_directory."/splitr.break";
my $splitr_seq = $output_directory."/splitr.seq";
my $splitr_log = $output_directory."/splitr.log";
$runner->run("$evalsplitalign_bin -c #<1 -a #<2 -b #>1 -q #>2", [$candidate_sequences,$candidate_alignments], [$splitr_break,$splitr_seq]);

my $denovo_break = $output_directory."/denovo.break";
my $denovo_seq = $output_directory."/denovo.seq";
my $denovo_log = $output_directory."/denovo.log";
if (lc($denovo_assembly) eq "yes")
{
	print "Denovo Breakpoint sequence prediction\n";
	$runner->submitter("direct");
	$runner->run("$split_seq_merge_script -c $config_filename -o $output_directory -s $submitter_type -p $max_parallel", [$clusters_sc_regions,$spanning_bam,$anchored_bam], [$splitr_break,$splitr_seq,$denovo_break,$denovo_seq]);
	$runner->submitter($submitter_type);
}
else
{
	$runner->run("touch #>1 #>2", [$clusters_sc_regions,$spanning_bam,$anchored_bam], [$denovo_break,$denovo_seq]);
}

print "Calculating spanning stats\n";
my $splitr_span_stats = $output_directory."/splitr.span.stats";
my $denovo_span_stats = $output_directory."/denovo.span.stats";
$runner->padd("$calc_span_stats_script -c #<1 -b #<2 -s #<3 > #>1", [$clusters_sc,$splitr_break,$splitr_seq], [$splitr_span_stats]);
$runner->padd("$calc_span_stats_script -c #<1 -b #<2 -s #<3 > #>1", [$clusters_sc,$denovo_break,$denovo_seq], [$denovo_span_stats]);
$runner->prun();

print "Calculating spanning p-values\n";
my $splitr_span_pval = $output_directory."/splitr.span.pval";
my $denovo_span_pval = $output_directory."/denovo.span.pval";
$runner->padd("$evaluate_fraglength_rscript #<1 #<2 $discord_read_trim #<3 #>1", [$read_stats,$cov_stats,$splitr_span_stats], [$splitr_span_pval]);
$runner->padd("$evaluate_fraglength_rscript #<1 #<2 $discord_read_trim #<3 #>1", [$read_stats,$cov_stats,$denovo_span_stats], [$denovo_span_pval]);
$runner->prun();

print "Calculating split read pvalues\n";
my $splitr_split_pval = $output_directory."/splitr.split.pval";
$runner->run("$evaluate_split_rscript #<1 #<2 #>1", [$cov_stats,$splitr_seq], [$splitr_split_pval]);

print "Selecting breakpoint sequence\n";
my $break = $output_directory."/break";
my $seq = $output_directory."/seq";
my $span_pval = $output_directory."/span.pval";
my $break_predict = $output_directory."/break.predict";
$runner->run("$select_breakpoint_seq_script -o $output_directory", [$splitr_span_pval,$splitr_break,$splitr_seq,$denovo_span_pval,$denovo_break,$denovo_seq], [$span_pval,$break,$seq,$break_predict]);

print "Creating fastas\n";
my $fusionreads_fasta = $output_directory."/fusionreads.fa";
my $breakpoints_fasta = $output_directory."/breakpoints.fa";
$runner->run("$create_read_fasta_script #<1 $reads_prefix > #>1", [$clusters_sc], [$fusionreads_fasta]);
$runner->run("cut -f1,2 #<1 | $make_fasta_script > #>1", [$seq], [$breakpoints_fasta]);

print "Splitting fastas\n";
my $fusionreads_split_prefix = $job_directory."/fusionreads";
my $fusionreads_split_catalog = $output_directory."/fusionreads.split.catalog";
my $breakpoints_split_prefix = $job_directory."/breakpoints";
my $breakpoints_split_catalog = $output_directory."/breakpoints.split.catalog";
$runner->run("$split_fasta_script #<1 $num_blat_sequences $fusionreads_split_prefix > #>1", [$fusionreads_fasta], [$fusionreads_split_catalog]);
$runner->run("$split_fasta_script #<1 $num_blat_sequences $breakpoints_split_prefix > #>1", [$breakpoints_fasta], [$breakpoints_split_catalog]);
my @fusionreads_split_fastas = readsplitcatalog($fusionreads_split_catalog);
my @breakpoints_split_fastas = readsplitcatalog($breakpoints_split_catalog);

print "Blat alignments\n";
my $fusionreads_genome_psl = $output_directory."/fusionreads.genome.psl";
my $fusionreads_cdna_psl = $output_directory."/fusionreads.cdna.psl";
my $fusionreads_rrna_psl = $output_directory."/fusionreads.rrna.psl";
my $breakpoints_genome_psl = $output_directory."/breakpoints.genome.psl";
my $breakpoints_genome_nointron_psl = $output_directory."/breakpoints.genome.nointron.psl";
my $breakpoints_cdna_psl = $output_directory."/breakpoints.cdna.psl";
my $breakpoints_est_psl = $output_directory."/breakpoints.est.psl";
my $breakpoints_exons_psl = $output_directory."/breakpoints.exons.psl";
my $breakpoints_cds_psl = $output_directory."/breakpoints.cds.psl";

my @blat_jobs;
push @blat_jobs, [$fusionreads_fasta, $fusionreads_genome_psl, \@fusionreads_split_fastas, $genome_fasta, "", "genome"];
push @blat_jobs, [$fusionreads_fasta, $fusionreads_cdna_psl, \@fusionreads_split_fastas, $cdna_fasta, "", "cdna"];
push @blat_jobs, [$fusionreads_fasta, $fusionreads_rrna_psl, \@fusionreads_split_fastas, $rrna_fasta, "", "rrna"];
push @blat_jobs, [$breakpoints_fasta, $breakpoints_genome_psl, \@breakpoints_split_fastas, $genome_fasta, "", "genome"];
push @blat_jobs, [$breakpoints_fasta, $breakpoints_genome_nointron_psl, \@breakpoints_split_fastas, $genome_fasta, "-maxIntron=0", "genome.nointron"];
push @blat_jobs, [$breakpoints_fasta, $breakpoints_cdna_psl, \@breakpoints_split_fastas, $cdna_fasta, "", "cdna"];
push @blat_jobs, [$breakpoints_fasta, $breakpoints_est_psl, \@breakpoints_split_fastas, $est_fasta, "", "est"];
push @blat_jobs, [$breakpoints_fasta, $breakpoints_exons_psl, \@breakpoints_split_fastas, $exons_fasta, "", "exons"];
push @blat_jobs, [$breakpoints_fasta, $breakpoints_cds_psl, \@breakpoints_split_fastas, $cds_fasta, "", "cds"];

doblatjobs(@blat_jobs);

print "Annotating fusions\n";
my $annotations_filename = $output_directory."/annotations";
my $mapping_stats_filename = $output_directory."/mapping.stats";
$runner->run("$annotate_fusions_script -c $config_filename -o $output_directory -n $library_name > #>1", [$span_pval,$break,$seq,$break_predict,$clusters_sc], [$annotations_filename]);
$runner->run("$calc_map_stats_script -c $config_filename -o $output_directory > #>1", [$clusters_sc], [$mapping_stats_filename]);

print "Coallating fusions\n";
my $results_filename = $output_directory."/results.tsv";
$runner->run("$coallate_fusions_script -c $config_filename -o $output_directory -l #<1 > #>1", [$annotations_filename,$mapping_stats_filename], [$results_filename]);

print "Running adaboost classifier\n";
my $results_classify = $output_directory."/results.classify.tsv";
$runner->run("$adaboost_rscript $positive_controls #<1 #>1", [$results_filename], [$results_classify]);

print "Filtering fusions\n";
my $filtered_results_filename = $output_directory."/results.filtered.tsv";
$runner->run("$filter_script probability '> $probability_threshold' < #<1 > #>1", [$results_classify], [$filtered_results_filename]);

print "Success\n";

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

sub splitalignmerge
{
	# Names of job products
	my @all_jobs_read_stats;
	my @all_jobs_cdna_pair_bam;
	my @all_jobs_spanning_bam;
	my @all_jobs_spanning_byread_bam;
	my @all_jobs_anchored_bam;
	my @all_jobs_expression;
	
	foreach my $split_fastq_prefix (@_)
	{
		# Job name
		my $job_name = basename($split_fastq_prefix);

		# Fastq files from split
		my $reads_end_1_split_fastq = $split_fastq_prefix.".1.fastq";
		my $reads_end_2_split_fastq = $split_fastq_prefix.".2.fastq";
		
		# Names of job products
		my $job_read_stats = $split_fastq_prefix.".concordant.read.stats";
		my $job_cdna_pair_bam = $split_fastq_prefix.".cdna.pair.bam";
		my $job_spanning_bam = $split_fastq_prefix.".spanning.bam";
		my $job_spanning_byread_bam = $split_fastq_prefix.".spanning.byread.bam";
		my $job_anchored_bam = $split_fastq_prefix.".anchored.bam";
		my $job_expression = $split_fastq_prefix.".expression.txt";
		
		push @all_jobs_read_stats, $job_read_stats;
		push @all_jobs_cdna_pair_bam, $job_cdna_pair_bam;
		push @all_jobs_spanning_bam, $job_spanning_bam;
		push @all_jobs_spanning_byread_bam, $job_spanning_byread_bam;
		push @all_jobs_anchored_bam, $job_anchored_bam;
		push @all_jobs_expression, $job_expression;
	
		my $job_cmd = $alignjob_script." ";
		$job_cmd .= "-c ".$config_filename." ";
		$job_cmd .= "-j ".$job_name." ";
		$job_cmd .= "-l ".$joblocal_directory." ";
		$job_cmd .= "-o ".$output_directory." ";
		$job_cmd .= "-n ".$library_name." ";
		$job_cmd .= "-p ".$split_fastq_prefix." ";
		
		$runner->padd("$job_cmd", [$reads_end_1_split_fastq, $reads_end_2_split_fastq], [$job_read_stats, $job_cdna_pair_bam, $job_spanning_bam, $job_spanning_byread_bam, $job_anchored_bam, $job_expression]); 
	}
	$runner->prun();
	
	# Merge products of the alignment process
	my $read_stats = $output_directory."/concordant.read.stats";
	my $expression = $output_directory."/expression.txt";
	my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
	my $spanning_bam = $output_directory."/spanning.bam";
	my $spanning_byread_bam = $output_directory."/spanning.byread.bam";
	my $anchored_bam = $output_directory."/anchored.bam";
	
	# Maximum number of files to merge at once
	my $merge_max = 100;
	
	# Merge merge_max at a time
	my @merge_jobs_read_stats;
	my @merge_jobs_expression;
	my @merge_jobs_cdna_pair_bam;
	my @merge_jobs_spanning_bam;
	my @merge_jobs_spanning_byread_bam;
	my @merge_jobs_anchored_bam;
	my @all_intermediate_read_stats;
	my @all_intermediate_expression;
	my @all_intermediate_cdna_pair_bam;
	my @all_intermediate_spanning_bam;
	my @all_intermediate_spanning_byread_bam;
	my @all_intermediate_anchored_bam;
	foreach my $job_index (0..$#all_jobs_read_stats)
	{
		my $job_read_stats = $all_jobs_read_stats[$job_index];
		my $job_expression = $all_jobs_expression[$job_index];
		my $job_cdna_pair_bam = $all_jobs_cdna_pair_bam[$job_index];
		my $job_spanning_bam = $all_jobs_spanning_bam[$job_index];
		my $job_spanning_byread_bam = $all_jobs_spanning_byread_bam[$job_index];
		my $job_anchored_bam = $all_jobs_anchored_bam[$job_index];
		
		push @merge_jobs_read_stats, $job_read_stats;
		push @merge_jobs_expression, $job_expression;
		push @merge_jobs_cdna_pair_bam, $job_cdna_pair_bam;
		push @merge_jobs_spanning_bam, $job_spanning_bam;
		push @merge_jobs_spanning_byread_bam, $job_spanning_byread_bam;
		push @merge_jobs_anchored_bam, $job_anchored_bam;
		
		if (scalar @merge_jobs_read_stats == $merge_max or $job_index == $#all_jobs_read_stats)
		{
			my $intermediate_index = scalar(@all_intermediate_read_stats) + 1;
			my $intermediate_prefix = "intermediate.".$intermediate_index;
			
			my $intermediate_read_stats = $job_directory."/$intermediate_prefix.concordant.read.stats";
			my $intermediate_expression = $job_directory."/$intermediate_prefix.expression.txt";
			my $intermediate_cdna_pair_bam = $job_directory."/$intermediate_prefix.cdna.pair.bam";
			my $intermediate_spanning_bam = $job_directory."/$intermediate_prefix.spanning.bam";
			my $intermediate_spanning_byread_bam = $job_directory."/$intermediate_prefix.spanning.byread.bam";
			my $intermediate_anchored_bam = $job_directory."/$intermediate_prefix.anchored.bam";
			
			# Merge read statistics
			$runner->padd("$merge_read_stats_script #<A > #>1", [@merge_jobs_read_stats], [$intermediate_read_stats]);
			
			# Merge expression statistics
			$runner->padd("$merge_expression_script #<A > #>1", [@merge_jobs_expression], [$intermediate_expression]);
			
			# Merge bam files or copy if theres only one
			if (scalar @merge_jobs_cdna_pair_bam == 1)
			{
				rename $merge_jobs_cdna_pair_bam[0], $intermediate_cdna_pair_bam;
				rename $merge_jobs_spanning_bam[0], $intermediate_spanning_bam;
				rename $merge_jobs_spanning_byread_bam[0], $intermediate_spanning_byread_bam;
				rename $merge_jobs_anchored_bam[0], $intermediate_anchored_bam;
			}
			else
			{
				$runner->padd("$samtools_bin merge #>1 #<A", [@merge_jobs_cdna_pair_bam], [$intermediate_cdna_pair_bam]);
				$runner->padd("$samtools_bin merge #>1 #<A", [@merge_jobs_spanning_bam], [$intermediate_spanning_bam]);
				$runner->padd("$samtools_bin merge -n #>1 #<A", [@merge_jobs_spanning_byread_bam], [$intermediate_spanning_byread_bam]);
				$runner->padd("$samtools_bin merge #>1 #<A", [@merge_jobs_anchored_bam], [$intermediate_anchored_bam]);
			}
			
			push @all_intermediate_read_stats, $intermediate_read_stats;
			push @all_intermediate_expression, $intermediate_expression;
			push @all_intermediate_cdna_pair_bam, $intermediate_cdna_pair_bam;
			push @all_intermediate_spanning_bam, $intermediate_spanning_bam;
			push @all_intermediate_spanning_byread_bam, $intermediate_spanning_byread_bam;
			push @all_intermediate_anchored_bam, $intermediate_anchored_bam;
			
			@merge_jobs_read_stats = ();
			@merge_jobs_expression = ();
			@merge_jobs_cdna_pair_bam = ();
			@merge_jobs_spanning_bam = ();
			@merge_jobs_spanning_byread_bam = ();
			@merge_jobs_anchored_bam = ();	
		}
	}
	$runner->prun();
	
	# Merge intermediate read statistics
	$runner->run("$merge_read_stats_script #<A > #>1", [@all_intermediate_read_stats], [$read_stats]);
	
	# Merge intermediate expression statistics
	$runner->run("$merge_expression_script #<A > #>1", [@all_intermediate_expression], [$expression]);
	
	# Merge bam files or copy if theres only one
	if (scalar @all_intermediate_spanning_bam == 1)
	{
		rename $all_intermediate_cdna_pair_bam[0], $cdna_pair_bam;
		rename $all_intermediate_spanning_bam[0], $spanning_bam;
		rename $all_intermediate_spanning_byread_bam[0], $spanning_byread_bam;
		rename $all_intermediate_anchored_bam[0], $anchored_bam;
	}
	else
	{
		$runner->run("$samtools_bin merge #>1 #<A", [@all_intermediate_cdna_pair_bam], [$cdna_pair_bam]);
		$runner->run("$samtools_bin merge #>1 #<A", [@all_intermediate_spanning_bam], [$spanning_bam]);
		$runner->run("$samtools_bin merge -n #>1 #<A", [@all_intermediate_spanning_byread_bam], [$spanning_byread_bam]);
		$runner->run("$samtools_bin merge #>1 #<A", [@all_intermediate_anchored_bam], [$anchored_bam]);
	}

	if (lc($remove_job_files) eq "yes")	
	{
		unlink @all_jobs_read_stats;
		unlink @all_jobs_expression;
		unlink @all_jobs_cdna_pair_bam;
		unlink @all_jobs_spanning_bam;
		unlink @all_jobs_spanning_byread_bam;
		unlink @all_jobs_anchored_bam;
		unlink @all_intermediate_read_stats;
		unlink @all_intermediate_expression;
		unlink @all_intermediate_cdna_pair_bam;
		unlink @all_intermediate_spanning_bam;
		unlink @all_intermediate_spanning_byread_bam;
		unlink @all_intermediate_anchored_bam;
	}
}

sub dosplitalignments
{
	my @split_fastq_prefixes = @{$_[0]};
	my $candidate_alignments = $_[1];
	
	my $candidate_reads_split_catalog = $output_directory."/candidates.split.catalog";
	$runner->run("$split_candidate_reads_script #<1 #<2 > #>1", [$reads_split_catalog,$candidate_reads], [$candidate_reads_split_catalog]);
	my @candidate_reads_split_filenames = readsplitcatalog($candidate_reads_split_catalog);

	scalar @split_fastq_prefixes == scalar @candidate_reads_split_filenames or die "Error: problem with candidate reads: $!\n";
	
	# Split alignments for each job
	my @all_job_candidate_alignments;
	
	foreach my $split_fastq_index (0..$#split_fastq_prefixes)
	{
		# Fastq files from split
		my $split_fastq_prefix = $split_fastq_prefixes[$split_fastq_index];
		my $reads_end_1_split_fastq = $split_fastq_prefix.".1.fastq";
		my $reads_end_2_split_fastq = $split_fastq_prefix.".2.fastq";
		
		# Candidate reads for each split
		my $candidate_reads_split_filename = $candidate_reads_split_filenames[$split_fastq_index];
		
		# Candidate read alignments
		my $job_candidate_alignments = $candidate_reads_split_filename.".alignments";
		push @all_job_candidate_alignments, $job_candidate_alignments;
		
		# Calculate split read alignment
		$runner->padd("$dosplitalign_bin -1 #<1 -2 #<2 -c #<3 -r #<4 -a #>1", [$reads_end_1_split_fastq,$reads_end_2_split_fastq,$candidate_sequences,$candidate_reads_split_filename], [$job_candidate_alignments]);
	}
	$runner->prun();

	# Merge the alignments
	if (scalar @all_job_candidate_alignments > 0)
	{
		$runner->run("cat #<A > #>1", [@all_job_candidate_alignments], [$candidate_alignments]);
	}
	else
	{
		$runner->run("touch #>1", [], [$candidate_alignments]);
	}
}

sub doblatjobs
{
	my @merge_jobs;
	foreach my $blat_job (@_)
	{
		my @blat_job_info = @{$blat_job};
		
		my $fasta = shift @blat_job_info;
		my $psl = shift @blat_job_info;
		if (not cmdrunner::uptodate([$fasta], [$psl]))
		{
			my @split_psls = padd_doblatalignments(@blat_job_info);
			push @merge_jobs, [\@split_psls, $psl];
		}
	}
	$runner->prun();
	
	foreach my $merge_job (@merge_jobs)
	{
		$runner->padd("cat #<A > #>1", $merge_job->[0], [$merge_job->[1]]);
	}
	$runner->prun();
}

sub padd_doblatalignments
{
	my $split_fastas_ref = shift;
	my $reference = shift;
	my $blat_parameters = shift;
	my $align_name = shift;
	
	if ($reference eq $genome_fasta)
	{
		return padd_doblatalignments_genome($split_fastas_ref, $blat_parameters, $align_name);
	}
	elsif ($reference eq $est_fasta)
	{
		return padd_doblatalignments_est($split_fastas_ref, $blat_parameters, $align_name);
	}
	
	my @output_psls;
	foreach my $split_fasta (@{$split_fastas_ref})
	{
		my $output_psl = $split_fasta.".".$align_name.".psl";
		$runner->padd("$blat_bin -noHead -stepSize=12 -tileSize=12 -minMatch=4 $blat_parameters $reference.2bit #<1 #>1", [$split_fasta], [$output_psl]);
		push @output_psls, $output_psl;
	}
	
	return @output_psls;
}

sub padd_doblatalignments_genome
{
	my $split_fastas_ref = shift;
	my $blat_parameters = shift;
	my $align_name = shift;
	
	my @output_psls;
	foreach my $chromosome (keys %chromosomes)
	{
		my $chromosome_fasta = $chromosome_prefix.".".$chromosome.".fa";
		push @output_psls, padd_doblatalignments($split_fastas_ref, $chromosome_fasta, $blat_parameters, $align_name.".".$chromosome);
	}
	
	return @output_psls;
}

sub padd_doblatalignments_est
{
	my $split_fastas_ref = shift;
	my $blat_parameters = shift;
	my $align_name = shift;
	
	my @est_split_fastas = readsplitcatalog($est_split_catalog);
	
	my @output_psls;
	foreach my $est_split_fasta_index (0..$#est_split_fastas)
	{
		my $est_split_fasta = $est_split_fastas[$est_split_fasta_index];
		push @output_psls, padd_doblatalignments($split_fastas_ref, $est_split_fasta, $blat_parameters, $align_name.".".$est_split_fasta_index);
	}
	
	return @output_psls;
}



