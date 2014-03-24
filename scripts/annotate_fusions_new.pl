#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

$| = 1;

use lib dirname($0);
use configdata;
use cmdrunner;
use gene_models;

use lib dirname($0)."/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --name      Library Name\n";

my $help;
my $config_filename;
my $output_directory;
my $library_name;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
	'name=s'      => \$library_name,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $output_directory or die @usage;
defined $library_name or die @usage;

my $config = configdata->new();
$config->read($config_filename);

# Config values
my $reference_fasta			= $config->get_value("reference_fasta");
my $gene_models_filename	= $config->get_value("gene_models");
my $genome_fasta			= $config->get_value("genome_fasta");
my $cdna_fasta				= $config->get_value("cdna_fasta");
my $exons_fasta				= $config->get_value("exons_fasta");
my $cds_fasta				= $config->get_value("cds_fasta");
my $rrna_fasta				= $config->get_value("rrna_fasta");
my $est_fasta				= $config->get_value("est_fasta");
my $est_alignments			= $config->get_value("est_alignments");
my $repeats_regions			= $config->get_value("repeats_regions");
my $splice_bias             = $config->get_value("splice_bias");
my $denovo_assembly			= $config->get_value("denovo_assembly");
my $tools_directory			= $config->get_value("tools_directory");
my $scripts_directory		= $config->get_value("scripts_directory");
my $samtools_bin			= $config->get_value("samtools_bin");
my $percident_threshold		= $config->get_value("percent_identity_threshold");

my $genome_max_ins = 2000;
my $est_max_ins = 10000;
my $cdna_max_ins = 10000000;

my $entropy_breakpoint_adjacent_size = 40;

# Create an indexable genome
my $genome_db = Bio::DB::Fasta->new($genome_fasta);

sub revcomp
{
	my $sequence = shift;
	my $revcomp = reverse($sequence);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

my $estislandsbin = $tools_directory."/estislands";
my $calc_map_stats_script = $scripts_directory."/calculate_mapping_stats.pl";
my $break_concordant_script = $scripts_directory."/calc_break_concordant.pl";
my $interrupted_script = $scripts_directory."/calc_interrupted.pl";

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/annotate";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("annotate");
$runner->prefix($log_prefix);

# Read in the read stats
my $read_stats = $output_directory."/concordant.read.stats";
my %read_stat_values;
get_stats($read_stats, \%read_stat_values);

# Read in the gene models
my $gene_models = gene_models->new($gene_models_filename);

# Main fusion info hash
my %fusion_info;

# Read in split read breakpoints 
my $break_filename = $output_directory."/splitreads.break";
read_breaks($break_filename, \%fusion_info);

# Read in split read sequences 
my $seq_filename = $output_directory."/splitreads.seq";
read_seq($seq_filename, \%fusion_info);

# Read in cluster info
my $clusters_sc = $output_directory."/clusters.sc";
read_cluster_info($clusters_sc, \%fusion_info);

# Filter fusions with an invalid sequence prediction
foreach my $cluster_id (keys %fusion_info)
{
	my $sequence = $fusion_info{$cluster_id}{sequence};
	my $break_in_seq = index $sequence, "|";
	
	if ($sequence eq "N" or $sequence eq "" or $break_in_seq < 0)
	{
		delete $fusion_info{$cluster_id};
	}
}

# Find fusion breakpoint information
foreach my $cluster_id (keys %fusion_info)
{
	for my $cluster_end ("0","1")
	{
		my $ref_name = $fusion_info{$cluster_id}{$cluster_end}{reference};
		my $breakpos = $fusion_info{$cluster_id}{$cluster_end}{breakpos};
		my $strand = $fusion_info{$cluster_id}{$cluster_end}{strand};
		my $alignregion = $fusion_info{$cluster_id}{$cluster_end}{alignregion};
		
		my $gene = $gene_models->calc_gene($ref_name, $breakpos);
		my $chromosome = $gene_models->calc_genomic_chromosome($ref_name);
		my $genomic_strand = $gene_models->calc_genomic_strand($ref_name, $strand);
		my $genomic_break_pos = $gene_models->calc_genomic_position($ref_name, $breakpos);
		my @genomic_regions = $gene_models->calc_genomic_regions($ref_name, $alignregion);
		my $gene_location = $gene_models->calc_gene_location($gene, $genomic_break_pos);
		my $gene_align_strand = ($genomic_strand eq $gene_models->{genes}{$gene}{strand}) ? "+" : "-";
		
		$fusion_info{$cluster_id}{$cluster_end}{gene} = $gene;
		$fusion_info{$cluster_id}{$cluster_end}{chromosome} = $chromosome;
		$fusion_info{$cluster_id}{$cluster_end}{genomic_strand} = $genomic_strand;
		$fusion_info{$cluster_id}{$cluster_end}{genomic_break_pos} = $genomic_break_pos;
		$fusion_info{$cluster_id}{$cluster_end}{genomic_regions} = [@genomic_regions];
		$fusion_info{$cluster_id}{$cluster_end}{gene_location} = $gene_location;
		$fusion_info{$cluster_id}{$cluster_end}{gene_align_strand} = $gene_align_strand;
	}
}

calculate_exonboundaries_feature_adjust_breakpos(\%fusion_info);

calculate_orf_feature(\%fusion_info);

calculate_expression_feature(\%fusion_info);

calculate_mapping_stats(\%fusion_info);

my $cdna_fasta_index = $cdna_fasta.".fai";
my $cdna_pair_sam = $output_directory."/cdna.pair.sam";
my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
my $cdna_pair_bam_bai = $cdna_pair_bam.".bai";
my $cdna_pair_bam_prefix = $cdna_pair_bam.".sort";
$runner->run("$samtools_bin view -bt $cdna_fasta_index #<1 | $samtools_bin sort -o - $cdna_pair_bam_prefix > #>1", [$cdna_pair_sam], [$cdna_pair_bam]);
$runner->padd("$samtools_bin index #<1 #>1", [$cdna_pair_bam], [$cdna_pair_bam_bai]);

calculate_splicing_index_feature($break_filename, \%fusion_info);

calculate_interrupted_index_feature($break_filename, \%fusion_info);

calculate_altsplice_feature(\%fusion_info);

calculate_break_adj_entropy_feature(\%fusion_info);

calculate_percident_feature(\%fusion_info);

calculate_breakpoint_homology_feature(\%fusion_info);

calculate_span_coverage_feature(\%fusion_info);

calculate_num_splice_variants_feature(\%fusion_info);

calculate_splice_score_feature(\%fusion_info);

calculate_repeat_proportion_feature(\%fusion_info);

foreach my $cluster_id (sort {$a <=> $b} keys %fusion_info)
{
	my $gene1 = $fusion_info{$cluster_id}{0}{gene};
	my $gene2 = $fusion_info{$cluster_id}{1}{gene};
	
	my $genome_strand1 = $fusion_info{$cluster_id}{0}{genomic_strand};
	my $genome_strand2 = $fusion_info{$cluster_id}{1}{genomic_strand};
	
	my $adjacent = "N";
	$adjacent = "Y" if defined $gene_models->{adjacent_gene}{$gene1}{$gene2};

	my $interchromosomal = "N";
	$interchromosomal = "Y" if $gene_models->{genes}{$gene1}{chromosome} ne $gene_models->{genes}{$gene2}{chromosome};

	my $inversion = "N";
	$inversion = "Y" if $interchromosomal eq "N" and $genome_strand1 eq $genome_strand2;

	my $eversion = "N";
	$eversion = "Y" if $interchromosomal eq "N" and $gene_models->{genes}{$gene1}{region}->[0] < $gene_models->{genes}{$gene2}{region}->[0] and $genome_strand1 eq "-" and $genome_strand2 eq "+";
	$eversion = "Y" if $interchromosomal eq "N" and $gene_models->{genes}{$gene1}{region}->[0] > $gene_models->{genes}{$gene2}{region}->[0] and $genome_strand1 eq "+" and $genome_strand2 eq "-";
	
	my $deletion = "N";
	$deletion = "Y" if $interchromosomal eq "N" and $inversion eq "N" and $eversion eq "N";

	my $read_through = "N";
	$read_through = "Y" if $deletion eq "Y" and $adjacent eq "Y";

	print $cluster_id."\tlibrary_name\t".$library_name."\n";

	print $cluster_id."\tgene1\t".$gene1."\n";
	print $cluster_id."\tgene_name1\t".$gene_models->{genes}{$gene1}{name}."\n";
	print $cluster_id."\tgene_chromosome1\t".$gene_models->{genes}{$gene1}{chromosome}."\n";
	print $cluster_id."\tgene_strand1\t".$gene_models->{genes}{$gene1}{strand}."\n";
	print $cluster_id."\tgene_start1\t".$gene_models->{genes}{$gene1}{region}->[0]."\n";
	print $cluster_id."\tgene_end1\t".$gene_models->{genes}{$gene1}{region}->[1]."\n";

	print $cluster_id."\tgene2\t".$gene2."\n";
	print $cluster_id."\tgene_name2\t".$gene_models->{genes}{$gene2}{name}."\n";
	print $cluster_id."\tgene_chromosome2\t".$gene_models->{genes}{$gene2}{chromosome}."\n";
	print $cluster_id."\tgene_strand2\t".$gene_models->{genes}{$gene2}{strand}."\n";
	print $cluster_id."\tgene_start2\t".$gene_models->{genes}{$gene2}{region}->[0]."\n";
	print $cluster_id."\tgene_end2\t".$gene_models->{genes}{$gene2}{region}->[1]."\n";
	
	print $cluster_id."\tadjacent\t".$adjacent."\n";
	print $cluster_id."\tinterchromosomal\t".$interchromosomal."\n";
	print $cluster_id."\tinversion\t".$inversion."\n";
	print $cluster_id."\teversion\t".$eversion."\n";
	print $cluster_id."\tdeletion\t".$deletion."\n";
	print $cluster_id."\tread_through\t".$read_through."\n";
	
	my @paired_features = (
	"gene_align_strand",
	"genomic_break_pos",
	"genomic_strand",
	"splicing_index",
	"interrupted_index",
	"span_coverage",
	"expression",
	"gene_location",
	"break_adj_entropy",
	"repeat_proportion",
	"repeat_list",
	);
	
	foreach my $feature (@paired_features)
	{
		foreach my $cluster_end ("0","1")
		{
			my $readable_cluster_end = $cluster_end + 1;
			print $cluster_id."\t".$feature.$readable_cluster_end."\t".$fusion_info{$cluster_id}{$cluster_end}{$feature}."\n";
		}
	}
	
	my @unpaired_features = (
	"orf",
	"exonboundaries",
	"altsplice",
	"span_count",
	"breakpoint_homology",
	"splice_score",
	"num_splice_variants",
	"genome_breakseqs_percident",
	"cdna_breakseqs_percident",
	"est_breakseqs_percident",
	"estisland_breakseqs_percident",
	"min_map_count",
	"max_map_count",
	"mean_map_count",
	"num_multi_map",
	);
	
	foreach my $feature (@unpaired_features)
	{
		print $cluster_id."\t".$feature."\t".$fusion_info{$cluster_id}{$feature}."\n";
	}
	
	print $cluster_id."\tbreak_adj_entropy_min\t".min($fusion_info{$cluster_id}{0}{break_adj_entropy},$fusion_info{$cluster_id}{1}{break_adj_entropy})."\n";
	print $cluster_id."\tspan_coverage_min\t".min($fusion_info{$cluster_id}{0}{span_coverage},$fusion_info{$cluster_id}{1}{span_coverage})."\n";
	print $cluster_id."\tspan_coverage_max\t".max($fusion_info{$cluster_id}{0}{span_coverage},$fusion_info{$cluster_id}{1}{span_coverage})."\n";
	print $cluster_id."\tmax_repeat_proportion\t".max($fusion_info{$cluster_id}{0}{repeat_proportion},$fusion_info{$cluster_id}{1}{repeat_proportion})."\n";
}

sub read_seq
{
	my $seqs_filename = shift;
	my $seqs_hash_ref = shift;
	
	open SEQ, $seqs_filename or die "Error: Unable to find $seqs_filename: $!\n";
	while (<SEQ>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		
		$seqs_hash_ref->{$cluster_id}{sequence} = $fields[1];
	}
	close SEQ;
}

sub read_breaks
{
	my $breaks_filename = shift;
	my $breaks_hash_ref = shift;
	
	open BR, $breaks_filename or die "Error: Unable to find $breaks_filename: $!\n";
	while (<BR>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $cluster_end = $fields[1];
		my $reference = $fields[2];
		my $strand = $fields[3];
		my $breakpos = $fields[4];
		
		$breaks_hash_ref->{$cluster_id}{$cluster_end}{reference} = $reference;
		$breaks_hash_ref->{$cluster_id}{$cluster_end}{strand} = $strand;
		$breaks_hash_ref->{$cluster_id}{$cluster_end}{breakpos} = $breakpos;
	}
	close BR;
}

sub read_clusters
{
	my $clusters_filename = shift;
	my $clusters_hash_ref = shift;
	
	open CLU, $clusters_filename or die "Error: Unable to find $clusters_filename: $!\n";
	while (<CLU>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $cluster_end = $fields[1];
		my $fragment_id = $fields[2];
		my $read_end = $fields[3];
		my $ref_name = $fields[4];
		my $strand = $fields[5];
		my $start = $fields[6];
		my $end = $fields[7];
		
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{read_id} = $fragment_id."/".$read_end;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{read_end} = $read_end;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{strand} = $strand;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{start} = $start;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{$fragment_id}{end} = $end;
	}
	close CLU;
}

sub find_breakseqs_overlap
{
	my $breakpoints_psl = shift;
	my $min_percident = shift;
	my $max_left_end_ref = shift;
	my $min_right_start_ref = shift;
	
	open PSL, $breakpoints_psl or die "Error: Unable to open $breakpoints_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
	
		my $num_matches = $psl_fields[0];
		my $num_target_bases_inserted = $psl_fields[7];
		my $strand = $psl_fields[8];
		my $cluster_id = $psl_fields[9];
		my $query_size = $psl_fields[10];
		my $query_start = $psl_fields[11] + 1;
		my $query_end = $psl_fields[12];
		my $target_seq_name = $psl_fields[13];
		my $target_size = $psl_fields[14];
		my $target_start = $psl_fields[15] + 1;
		my $target_end = $psl_fields[16];
		my @block_sizes = split /,/, $psl_fields[18];
		my @block_query_starts = split /,/, $psl_fields[19];
		my @block_target_starts = split /,/, $psl_fields[20];
		
		my $percident = $num_matches / ($query_end - $query_start + 1);
		
		next if $percident < $min_percident;
		
		if ($query_start == 1)
		{
			$max_left_end_ref->{$cluster_id} = $query_end if not defined $max_left_end_ref->{$cluster_id};
			$max_left_end_ref->{$cluster_id} = max($max_left_end_ref->{$cluster_id}, $query_end);
		}
		
		if ($query_end == $query_size)
		{
			$min_right_start_ref->{$cluster_id} = $query_start if not defined $min_right_start_ref->{$cluster_id};
			$min_right_start_ref->{$cluster_id} = min($min_right_start_ref->{$cluster_id}, $query_start);
		}
	}
}

sub find_breakseqs_estislands_percident
{
	my $genome_breakpoints_psl = shift;
	my $max_percident_ref = shift;
	
	my $estislands_psl = $genome_breakpoints_psl.".estisl.psl";	
	$runner->run("$estislandsbin -b #<1 -e $est_alignments -o #>1", [$genome_breakpoints_psl], [$estislands_psl]);
	
	open PSL, $estislands_psl or die "Error: Unable to open $estislands_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
	
		my $num_matches = $psl_fields[0];
		my $cluster_id = $psl_fields[9];
		my $breakpoint_seq_length = $psl_fields[10];
	
		my $percent_identity = $num_matches / $breakpoint_seq_length;
		
		$max_percident_ref->{$cluster_id} = 0 if not defined $max_percident_ref->{$cluster_id};
		$max_percident_ref->{$cluster_id} = max($max_percident_ref->{$cluster_id}, $percent_identity);
	}
	close PSL;
}

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

sub read_expression
{
	my $expression_filename = shift;
	my $expression_ref = shift;

	open EXPR, $expression_filename or die "Error: Unable to open $expression_filename: $!\n";
	while (<EXPR>)
	{
		chomp;
		my @fields = split /\t/;
	
		my $ensgene = $fields[0];
		
		$expression_ref->{$ensgene}{expression} = $fields[1];
	}
	close EXPR;
}

sub calculate_mapping_stats
{
	my $fusion_info_ref = shift;
	
	my $mapping_stats_filename = $output_directory."/mapping.stats";
	$runner->run("$calc_map_stats_script -c $config_filename -o $output_directory > #>1", [$clusters_sc], [$mapping_stats_filename]);
	
	open MPS, $mapping_stats_filename or die "Error: Unable to find $mapping_stats_filename: $!\n";
	while (<MPS>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $anno_type = $fields[1];
		my $anno_value = $fields[2];
		
		next if not defined $fusion_info_ref->{$cluster_id};
		
		$fusion_info_ref->{$cluster_id}{$anno_type} = $anno_value;
	}
	close MPS;
}


# 
# Fusion splice index feature
#

sub calculate_splicing_index_feature
{
	my $breaks_filename = shift;
	my $fusion_info_ref = shift;

	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		foreach my $cluster_end ("0","1")
		{
			$fusion_info_ref->{$cluster_id}{$cluster_end}{splicing_index} = "-";
		}
	}
	
	my $break_concordant_filename = $breaks_filename.".concordant.counts";

	$runner->run("$break_concordant_script -c $config_filename -o $output_directory -b #<1 > #>1", [$breaks_filename], [$break_concordant_filename]);
	
	open BRC, $break_concordant_filename or die "Error: Unable to open $break_concordant_filename: $!\n";
	while (<BRC>)
	{
		chomp;
		my ($cluster_id, $cluster_end, $break_concordant) = split /\t/;
		
		$fusion_info_ref->{$cluster_id}{$cluster_end}{splicing_index} = $break_concordant / $fusion_info_ref->{$cluster_id}{span_count};
	}
	close BRC;
}


# 
# Interrupted index feature
#

sub calculate_interrupted_index_feature
{
	my $breaks_filename = shift;
	my $fusion_info_ref = shift;
	
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		foreach my $cluster_end ("0","1")
		{
			$fusion_info_ref->{$cluster_id}{$cluster_end}{interrupted_index} = "-";
		}
	}
	
	my $interrupted_filename = $breaks_filename.".interrupted.counts";
	
	$runner->run("$interrupted_script -c $config_filename -o $output_directory -b #<1 > #>1", [$breaks_filename], [$interrupted_filename]);
	
	open BIN, $interrupted_filename or die "Error: Unable to open $interrupted_filename: $!\n";
	while (<BIN>)
	{
		chomp;
		my ($cluster_id, $cluster_end, $gene, $size_before, $size_after, $count_before, $count_after) = split /\t/;
		
		my $expression_before = $count_before / ($size_before + 1) + 1;
		my $expression_after = $count_after / ($size_after + 1) + 1;
		
		$fusion_info_ref->{$cluster_id}{$cluster_end}{interrupted_index} = $expression_before / $expression_after;
	}
	close BIN;
}

# Find the combined length of a set of regions
sub regions_length
{
	my @regions = @_;

	my $length = 0;
	foreach my $region (@regions)
	{
		$length += $region->[1] - $region->[0] + 1;
	}

	return $length;
}

# Merge overlapping regions
sub merge_regions
{
	my @regions = @_;
	my @merged;

	@regions = sort { $a->[0] <=> $b->[0] } (@regions);

	my $merged_start;
	my $merged_end;
	foreach my $region (@regions)
	{
		$merged_start = $region->[0] if not defined $merged_start;
		$merged_end = $region->[1] if not defined $merged_end;

		if ($region->[0] > $merged_end + 1)
		{
			push @merged, [$merged_start, $merged_end];

			$merged_start = $region->[0];
			$merged_end = $region->[1];
		}
		else
		{
			$merged_end = max($merged_end, $region->[1]);
		}
	}
	push @merged, [$merged_start, $merged_end];

	return @merged;
}

sub overlapsize
{
	my $region1 = shift;
	my $region2 = shift;
	
	my $size1 = $region1->[1] - $region1->[0] + 1;
	my $size2 = $region2->[1] - $region2->[0] + 1;
	
	my $overlap = min($region2->[1] - $region1->[0] + 1, $region1->[1] - $region2->[0] + 1, $size1, $size2);
	
	return max(0, $overlap);
}

sub get_bins
{
	my $start = shift;
	my $end = shift;
	my $bin_spacing = shift;

	my $start_bin = int($start / $bin_spacing);
	my $end_bin = int($end / $bin_spacing);
	
	return ($start_bin .. $end_bin);
}

sub read_repeats
{
	my $repeats_filename = shift;
	my $repeats_ref = shift;
	
	# Read in repeats
	my @repeat_list;
	my $length_sum = 0;
	open REP, $repeats_filename or die "Error: Unable to open $repeats_filename: $!\n";
	while (<REP>)
	{
		chomp;
		my ($chromosome,$start,$end,$type) = split /\t/;
		
		push @repeat_list, [$chromosome,$start,$end,$type];
		$length_sum += $end - $start;
	}
	close REP;
	
	my $length_mean = $length_sum / scalar @repeat_list;
	
	# Bin spacing so that most repeats are contained in only a few bins
	$repeats_ref->{bin_spacing} = int($length_mean * 5);

	# Binning for gene lookup
	foreach my $repeat (@repeat_list)
	{
		foreach my $bin (get_bins($repeat->[1],$repeat->[2],$repeats_ref->{bin_spacing}))
		{
			push @{$repeats_ref->{binned_repeats}{$repeat->[0]}{$bin}}, [$repeat->[1],$repeat->[2],$repeat->[3]];
		}
	}
}

sub overlap
{
	my $region1 = $_[0];
	my $region2 = $_[1];
	
	if ($region1->[1] < $region2->[0] or $region1->[0] > $region2->[1])
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

sub read_cluster_info
{
	my $clusters_filename = shift;
	my $fusion_info_ref = shift;
	
	my $current_cluster_id;
	my %align_starts;
	my %align_ends;
	my %matched;
	open CLU, $clusters_filename or die "Error: Unable to open $clusters_filename: $!\n";
	while (<CLU>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $cluster_end = $fields[1];
		my $start = $fields[6];
		my $end = $fields[7];
		
		if (defined $current_cluster_id and $current_cluster_id != $cluster_id)
		{
			$fusion_info_ref->{$current_cluster_id}{span_count} = scalar(@{$align_ends{0}});
			
			foreach my $current_cluster_end (keys %align_starts)
			{
				my $align_region_start = min(@{$align_starts{$current_cluster_end}});
				my $align_region_end = max(@{$align_ends{$current_cluster_end}});
				my $num_matched = scalar keys %{$matched{$current_cluster_end}};
				
				$fusion_info_ref->{$current_cluster_id}{$current_cluster_end}{alignregion} = [$align_region_start,$align_region_end];
				$fusion_info_ref->{$current_cluster_id}{$current_cluster_end}{num_matched} = $num_matched;
			}
			
			%align_starts = ();
			%align_ends = ();
			%matched = ();
		}
		
		push @{$align_starts{$cluster_end}}, $start;
		push @{$align_ends{$cluster_end}}, $end;
		
		foreach my $pos ($start..$end)
		{
			$matched{$cluster_end}{$pos} = 1;
		}
		
		$current_cluster_id = $cluster_id;
	}
	close CLU;
	
	if (defined $current_cluster_id)
	{
		$fusion_info_ref->{$current_cluster_id}{span_count} = scalar(@{$align_ends{0}});
		
		foreach my $current_cluster_end (keys %align_starts)
		{
			my $align_region_start = min(@{$align_starts{$current_cluster_end}});
			my $align_region_end = max(@{$align_ends{$current_cluster_end}});
			my $num_matched = scalar keys %{$matched{$current_cluster_end}};
			
			$fusion_info_ref->{$current_cluster_id}{$current_cluster_end}{alignregion} = [$align_region_start,$align_region_end];
			$fusion_info_ref->{$current_cluster_id}{$current_cluster_end}{num_matched} = $num_matched;
		}
	}
}

sub find_alignregion
{
	my $breakpoints_psl = shift;
	my $align_ref = shift;
	
	open PSL, $breakpoints_psl or die "Error: Unable to open $breakpoints_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
	
		my $num_matches = $psl_fields[0];
		my $num_target_bases_inserted = $psl_fields[7];
		my $strand = $psl_fields[8];
		my $cluster_id = $psl_fields[9];
		my $query_size = $psl_fields[10];
		my $target_seq_name = $psl_fields[13];
		my $target_size = $psl_fields[14];
		my @block_sizes = split /,/, $psl_fields[18];
		my @query_starts = split /,/, $psl_fields[19];
		my @target_starts = split /,/, $psl_fields[20];
	
		$target_seq_name =~ /([^|]*)\|([^|]*)\|?([^|]*)?/;
		my $gene = $1;

		warn $target_seq_name." not found\n" if not defined $gene;	
		next unless $gene eq $fusion_info{$cluster_id}{0}{gene} or $gene eq $fusion_info{$cluster_id}{1}{gene};
		
		foreach my $block_index (0..$#block_sizes)
		{
			my $block_size = $block_sizes[$block_index];
			
			my $query_start = $query_starts[$block_index] + 1;
			my $query_end = $query_starts[$block_index] + $block_size;

			if ($strand eq "-")
			{
				$query_start = $query_size - $query_starts[$block_index] - $block_size + 1;
				$query_end = $query_size - $query_starts[$block_index];
			}

			my $target_start = $target_starts[$block_index] + 1;
			my $target_end = $target_starts[$block_index] + $block_size;
			
			push @{$align_ref->{$cluster_id}{$gene}{align_strand}}, $strand;
			push @{$align_ref->{$cluster_id}{$gene}{query_region}}, [$query_start, $query_end, $query_size];
			push @{$align_ref->{$cluster_id}{$gene}{target_region}}, [$target_start, $target_end, $target_size, $target_seq_name];
		}
	}
	close PSL;
}

#
# Expression feature
# 

sub calculate_expression_feature
{
	my $fusion_info_ref = shift;
	
	my %expression;
	my $expression_filename = $output_directory."/expression.txt";
	read_expression($expression_filename, \%expression);

	foreach my $cluster_id (keys %fusion_info)
	{
		foreach my $cluster_end ("0","1")
		{
			my $gene = $fusion_info_ref->{$cluster_id}{$cluster_end}{gene};
			
			if (not defined $expression{$gene})
			{
				$fusion_info_ref->{$cluster_id}{$cluster_end}{expression} = 0;
			}
			else
			{
				$fusion_info_ref->{$cluster_id}{$cluster_end}{expression} = $expression{$gene}{expression};
			}
		}
	}
}


#
# Altsplice feature
#

sub calculate_altsplice_feature
{
	my $fusion_info_ref = shift;
	
	my $breakpoints_genome_psl = $output_directory."/breakpoints.genome.psl";
	
	my %genome_breakseqs_percident;
	find_breakseqs_percident($breakpoints_genome_psl, \%genome_breakseqs_percident);
	
	foreach my $cluster_id (keys %fusion_info)
	{
		if (defined $genome_breakseqs_percident{$cluster_id} and $genome_breakseqs_percident{$cluster_id} > $percident_threshold)
		{
			$fusion_info_ref->{$cluster_id}{altsplice} = "Y";
		}
		else
		{
			$fusion_info_ref->{$cluster_id}{altsplice} = "N";
		}
	}
}

# 
# Exon boundaries and adjust breakpos
#

sub calculate_exonboundaries_feature_adjust_breakpos
{
	my $fusion_info_ref = shift;
	
	my $breakpoints_exons_psl = $output_directory."/breakpoints.exons.psl";
	
	my %exon_align;
	find_alignregion($breakpoints_exons_psl, \%exon_align);
	
	foreach my $cluster_id (keys %fusion_info)
	{
		my $gene1 = $fusion_info_ref->{$cluster_id}{0}{gene};
		my $gene2 = $fusion_info_ref->{$cluster_id}{1}{gene};
		
		$fusion_info_ref->{$cluster_id}{exonboundaries} = "N";
		
		foreach my $start_index1 (0..$#{$exon_align{$cluster_id}{$gene1}{align_strand}})
		{
			my $strand1 = $exon_align{$cluster_id}{$gene1}{align_strand}->[$start_index1];
			my $query_region1 = $exon_align{$cluster_id}{$gene1}{query_region}->[$start_index1];
			my $target_region1 = $exon_align{$cluster_id}{$gene1}{target_region}->[$start_index1];
	
			my $target_size1 = $target_region1->[2];
			my $target_name1 = $target_region1->[3];
	
			foreach my $start_index2 (0..$#{$exon_align{$cluster_id}{$gene2}{align_strand}})
			{	
				my $strand2 = $exon_align{$cluster_id}{$gene2}{align_strand}->[$start_index2];
				my $query_region2 = $exon_align{$cluster_id}{$gene2}{query_region}->[$start_index2];
				my $target_region2 = $exon_align{$cluster_id}{$gene2}{target_region}->[$start_index2];
	
				my $target_size2 = $target_region2->[2];
				my $target_name2 = $target_region2->[3];
				
				my $query_end1_start2 = ($query_region1->[1] + 1 == $query_region2->[0]);
				my $query_end2_start1 = ($query_region2->[1] + 1 == $query_region1->[0]);
				
				if ($query_end1_start2)
				{
					my $query_end1_boundary;
					my $query_end1_targetpos;
					if ($strand1 eq "+")
					{
						$query_end1_boundary = $target_region1->[1] == $target_size1;
						$query_end1_targetpos = $target_size1;
					}
					else
					{
						$query_end1_boundary = $target_region1->[0] == 1;
						$query_end1_targetpos = 1;
					}
					
					my $query_start2_boundary;
					my $query_start2_targetpos;
					if ($strand2 eq "+")
					{
						$query_start2_boundary = $target_region2->[0] == 1;
						$query_start2_targetpos = 1;
					}
					else
					{
						$query_start2_boundary = $target_region2->[1] == $target_size2;
						$query_start2_targetpos = $target_size2;
					}
					
					if ($query_end1_boundary and $query_start2_boundary)
					{
						$fusion_info_ref->{$cluster_id}{exonboundaries} = "Y";
						
						$fusion_info_ref->{$cluster_id}{0}{genomic_break_pos} = $gene_models->exon_to_genome($target_name1, $query_end1_targetpos);
						$fusion_info_ref->{$cluster_id}{1}{genomic_break_pos} = $gene_models->exon_to_genome($target_name2, $query_start2_targetpos);
						
						last;
					}
				}
				elsif ($query_end2_start1)
				{
					my $query_end2_boundary;
					my $query_end2_targetpos;
					if ($strand2 eq "+")
					{
						$query_end2_boundary = $target_region2->[1] == $target_size2;
						$query_end2_targetpos = $target_size2;
					}
					else
					{
						$query_end2_boundary = $target_region2->[0] == 1;
						$query_end2_targetpos = 1;
					}
					
					my $query_start1_boundary;
					my $query_start1_targetpos;
					if ($strand1 eq "+")
					{
						$query_start1_boundary = $target_region1->[0] == 1;
						$query_start1_targetpos = 1;
					}
					else
					{
						$query_start1_boundary = $target_region1->[1] == $target_size1;
						$query_start1_targetpos = $target_size1;
					}
					
					if ($query_end2_boundary and $query_start1_boundary)
					{
						$fusion_info_ref->{$cluster_id}{exonboundaries} = "Y";
						
						$fusion_info_ref->{$cluster_id}{0}{genomic_break_pos} = $gene_models->exon_to_genome($target_name1, $query_start1_targetpos);
						$fusion_info_ref->{$cluster_id}{1}{genomic_break_pos} = $gene_models->exon_to_genome($target_name2, $query_end2_targetpos);
						
						last;
					}
				}
			}
			
			last if $fusion_info_ref->{$cluster_id}{exonboundaries} eq "Y";
		}
	}
}

# 
# ORF feature
#

sub calculate_orf_feature
{
	my $fusion_info_ref = shift;
	
	my $breakpoints_cds_psl = $output_directory."/breakpoints.cds.psl";
	
	my %cds_align;
	find_alignregion($breakpoints_cds_psl, \%cds_align);
	
	foreach my $cluster_id (keys %fusion_info)
	{
		my $gene1 = $fusion_info_ref->{$cluster_id}{0}{gene};
		my $gene2 = $fusion_info_ref->{$cluster_id}{1}{gene};
		
		$fusion_info_ref->{$cluster_id}{orf} = "N";
		
		foreach my $start_index1 (0..$#{$cds_align{$cluster_id}{$gene1}{align_strand}})
		{
			my $strand1 = $cds_align{$cluster_id}{$gene1}{align_strand}->[$start_index1];
			my $query_region1 = $cds_align{$cluster_id}{$gene1}{query_region}->[$start_index1];
			my $target_region1 = $cds_align{$cluster_id}{$gene1}{target_region}->[$start_index1];
		
			foreach my $start_index2 (0..$#{$cds_align{$cluster_id}{$gene2}{align_strand}})
			{
				my $strand2 = $cds_align{$cluster_id}{$gene2}{align_strand}->[$start_index2];
				my $query_region2 = $cds_align{$cluster_id}{$gene2}{query_region}->[$start_index2];
				my $target_region2 = $cds_align{$cluster_id}{$gene2}{target_region}->[$start_index2];
	
				next unless $strand1 eq $strand2;
	
				my $query_phase = ($query_region1->[0] - $query_region2->[0]) % 3;
				
				my $target_phase;
				if ($strand1 eq "+")
				{
					$target_phase = ($target_region1->[0] - $target_region2->[0]) % 3;
				}
				else
				{
					$target_phase = ($target_region1->[1] - $target_region2->[1]) % 3;
				}
	
				$fusion_info_ref->{$cluster_id}{orf} = "Y" if $query_phase == $target_phase;
			}
		}
	}
}


# 
# Entropy feature
#

sub calculate_entropy
{
	my $seq = shift;
	
	my $entropy = 0;
	foreach my $n1 ("A","C","T","G")
	{
		foreach my $n2 ("A","C","T","G")
		{
			my $npair = $n1.$n2;
			my $count = 0; $count++ while $seq =~ /$npair/g;
			
			next if $count == 0;
			
			my $p = $count / (length($seq) - 1);
			my $logp = log($p)/log(2);
			
			$entropy -= $p * $logp;
		}
	}
	
	return $entropy;
}

sub calculate_break_adj_entropy_feature
{
	my $fusion_info_ref = shift;
	
	my $entropy_breakpoint_adjacent_size = 40;
	
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		my $sequence = $fusion_info_ref->{$cluster_id}{sequence};
		my $break_in_seq = index $sequence, "|";
		$sequence =~ s/\|//g;
		
		my $break_adj_seq1 = substr $sequence, max(0, $break_in_seq - $entropy_breakpoint_adjacent_size), min($break_in_seq, $entropy_breakpoint_adjacent_size);
		my $break_adj_seq2 = substr $sequence, $break_in_seq, min(length($sequence) - $break_in_seq, $entropy_breakpoint_adjacent_size);
		
		$fusion_info_ref->{$cluster_id}{0}{break_adj_entropy} = calculate_entropy($break_adj_seq1);
		$fusion_info_ref->{$cluster_id}{1}{break_adj_entropy} = calculate_entropy($break_adj_seq2);
	}
}

#
# Percent Identity Feature
#

sub find_breakseqs_percident
{
	my $breakpoints_psl = shift;
	my $max_percident_ref = shift;
	my $max_ins = shift;
	
	open PSL, $breakpoints_psl or die "Error: Unable to open $breakpoints_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
	
		my $num_matches = $psl_fields[0];
		my $num_target_bases_inserted = $psl_fields[7];
		my $cluster_id = $psl_fields[9];
		my $breakpoint_seq_length = $psl_fields[10];
	
		next if defined $max_ins and $num_target_bases_inserted > $max_ins;
	
		my $percent_identity = $num_matches / $breakpoint_seq_length;
		
		$max_percident_ref->{$cluster_id} = 0 if not defined $max_percident_ref->{$cluster_id};
		$max_percident_ref->{$cluster_id} = max($max_percident_ref->{$cluster_id}, $percent_identity);
	}
	close PSL;
}

sub calculate_percident_feature
{
	my $fusion_info_ref = shift;
	
	my $breakpoints_genome_psl = $output_directory."/breakpoints.genome.psl";
	my $breakpoints_cdna_psl = $output_directory."/breakpoints.cdna.psl";
	my $breakpoints_est_psl = $output_directory."/breakpoints.est.psl";
	
	my %breakseqs_percident;
	find_breakseqs_percident($breakpoints_genome_psl, \%{$breakseqs_percident{genome}}, $genome_max_ins);
	find_breakseqs_percident($breakpoints_cdna_psl, \%{$breakseqs_percident{cdna}}, $cdna_max_ins);
	find_breakseqs_percident($breakpoints_est_psl, \%{$breakseqs_percident{est}}, $est_max_ins);
	find_breakseqs_estislands_percident($breakpoints_genome_psl, \%{$breakseqs_percident{estisland}});
	
	foreach my $ref_type ("genome", "cdna", "est", "estisland")
	{
		foreach my $cluster_id (keys %{$fusion_info_ref})
		{
			$fusion_info_ref->{$cluster_id}{$ref_type."_breakseqs_percident"} = 0;
			
			next if not defined $breakseqs_percident{$ref_type}{$cluster_id};
			
			my $sequence = $fusion_info_ref->{$cluster_id}{sequence};
			my $break_in_seq = index $sequence, "|";
			$sequence =~ s/\|//g;
			
			my $seq_length = length($sequence);
			my $seq1_length = $break_in_seq;
			my $seq2_length = length($sequence) - $break_in_seq;
			
			my $mismatches = (1 - $breakseqs_percident{$ref_type}{$cluster_id}) * $seq_length;
			my $adjusted_percident = 1 - ($mismatches / min($seq1_length, $seq2_length));
			$adjusted_percident = max(0,$adjusted_percident);
		
			$fusion_info_ref->{$cluster_id}{$ref_type."_breakseqs_percident"} = $adjusted_percident;
		}
	}
}


#
# Breakpoint Homology Feature
#

sub calculate_breakpoint_homology_feature
{
	my $fusion_info_ref = shift;
	
	my $breakpoints_genome_nointron_psl = $output_directory."/breakpoints.genome.nointron.psl";
	my $breakpoints_cdna_psl = $output_directory."/breakpoints.cdna.psl";
	
	my %max_left_end;
	my %min_right_start;
	find_breakseqs_overlap($breakpoints_genome_nointron_psl, $percident_threshold, \%max_left_end, \%min_right_start);
	find_breakseqs_overlap($breakpoints_cdna_psl, $percident_threshold, \%max_left_end, \%min_right_start);
	
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		$fusion_info_ref->{$cluster_id}{breakpoint_homology} = 0;
		if (defined $max_left_end{$cluster_id} and defined $min_right_start{$cluster_id})
		{
			$fusion_info_ref->{$cluster_id}{breakpoint_homology} = max(0, $max_left_end{$cluster_id} - $min_right_start{$cluster_id} + 1);
		}
	}
}


#
# Coverage Feature
#

sub calculate_span_coverage_feature
{
	my $fusion_info_ref = shift;
	
	my $minimum_coverage = $read_stat_values{"fraglength_mean"} - $read_stat_values{"readlength_min"};
	
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		foreach my $cluster_end ("0","1")
		{
			my $num_matched = $fusion_info_ref->{$cluster_id}{$cluster_end}{num_matched};
			my $normalized_coverage = $num_matched / $minimum_coverage;
			$fusion_info_ref->{$cluster_id}{$cluster_end}{span_coverage} = $normalized_coverage;
		}
	}
}


#
# Number of Splice Variants Feature
#

sub calculate_num_splice_variants_feature
{
	my $fusion_info_ref = shift;
	
	my %fusion_splice_variants;
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		my $gene1 = $fusion_info_ref->{$cluster_id}{0}{gene};
		my $gene2 = $fusion_info_ref->{$cluster_id}{1}{gene};
		
		my $strand1 = $fusion_info_ref->{$cluster_id}{0}{genomic_strand};
		my $strand2 = $fusion_info_ref->{$cluster_id}{1}{genomic_strand};
		
		my $breakpos1 = $fusion_info_ref->{$cluster_id}{0}{genomic_break_pos};
		my $breakpos2 = $fusion_info_ref->{$cluster_id}{1}{genomic_break_pos};
		
		my $gene_strand_a = ($gene1 lt $gene2) ? $gene1.$strand1 : $gene2.$strand2;
		my $gene_strand_b = ($gene1 lt $gene2) ? $gene2.$strand2 : $gene1.$strand1;
		my $breakpos_a = ($gene1 lt $gene2) ? $breakpos1 : $breakpos2;
		my $breakpos_b = ($gene1 lt $gene2) ? $breakpos2 : $breakpos1;
		
		$fusion_splice_variants{$gene_strand_a}{$gene_strand_b}{$breakpos_a."-".$breakpos_b} = 1;
	}
	
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		my $gene1 = $fusion_info_ref->{$cluster_id}{0}{gene};
		my $gene2 = $fusion_info_ref->{$cluster_id}{1}{gene};
		
		my $strand1 = $fusion_info_ref->{$cluster_id}{0}{genomic_strand};
		my $strand2 = $fusion_info_ref->{$cluster_id}{1}{genomic_strand};
		
		my $gene_strand_a = ($gene1 lt $gene2) ? $gene1.$strand1 : $gene2.$strand2;
		my $gene_strand_b = ($gene1 lt $gene2) ? $gene2.$strand2 : $gene1.$strand1;
		
		$fusion_info_ref->{$cluster_id}{num_splice_variants} = scalar keys %{$fusion_splice_variants{$gene_strand_a}{$gene_strand_b}};
	}
}


#
# Splice Score Feature
#

sub get_splice_seq
{
	my $fusion_info_end_ref = shift;
	
	my $chromosome = $fusion_info_end_ref->{chromosome};
	my $strand = $fusion_info_end_ref->{genomic_strand};
	my $position = $fusion_info_end_ref->{genomic_break_pos};
	
	if ($strand eq "+")
	{
		return "NN" if $position + 1 < 1;
		return "NN" if $position + 2 > $genome_db->length($chromosome);
		
		return $genome_db->seq($chromosome, $position + 1, $position + 2);
	}
	elsif ($strand eq "-")
	{
		return "NN" if $position - 2 < 1;
		return "NN" if $position - 1 > $genome_db->length($chromosome);
		
		return revcomp($genome_db->seq($chromosome, $position - 2, $position - 1));
	}
}

sub calc_edit_dist
{
	my $seq1 = shift;
	my $seq2 = shift;
	
	warn "Error: expected $seq1 and $seq2 to be same size\n" if length($seq1) != length$seq2;
	
	my @nt1 = split //, $seq1;
	my @nt2 = split //, $seq2;
	
	my $dist = 0;
	foreach my $ntindex (0..min($#nt1,$#nt2))
	{
		$dist++ if $nt1[$ntindex] ne $nt2[$ntindex];
	}
	
	return $dist;
}

sub calculate_splice_score_feature
{
	my $fusion_info_ref = shift;
	
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		my $splice_seq1 = get_splice_seq($fusion_info_ref->{$cluster_id}{0});
		my $splice_seq2 = get_splice_seq($fusion_info_ref->{$cluster_id}{1});
		
		my $splice_seqf = $splice_seq1.revcomp($splice_seq2);
		my $splice_seqr = $splice_seq2.revcomp($splice_seq1);
		
		my @splice_scores;
		push @splice_scores, calc_edit_dist("GTAG", $splice_seqf);
		push @splice_scores, calc_edit_dist("GTAG", $splice_seqr);
		push @splice_scores, calc_edit_dist("ATAC", $splice_seqf);
		push @splice_scores, calc_edit_dist("ATAC", $splice_seqr);
		
		$fusion_info_ref->{$cluster_id}{splice_score} = 4 - min(@splice_scores);
	}
}


#
# Repeat Proportion Feature
#

sub get_repeat_proportion
{
	my $chromosome = shift;
	my $regions = shift;
	my $repeats_ref = shift;

	my @repeat_overlaps;	
	foreach my $region (@{$regions})
	{
		foreach my $bin (get_bins($region->[0],$region->[1],$repeats_ref->{bin_spacing}))
		{
			foreach my $repeat (@{$repeats_ref->{binned_repeats}{$chromosome}{$bin}})
			{
				if (overlap($repeat,$region))
				{
					push @repeat_overlaps, overlapsize($repeat,$region);
				}
			}
		}
	}
	
	my $repeat_proportion = max(0,@repeat_overlaps) / regions_length(@{$regions});
	
	return $repeat_proportion;
}

sub get_repeat_list
{
	my $chromosome = shift;
	my $regions = shift;
	my $repeats_ref = shift;

	my @repeat_list;
	foreach my $region (@{$regions})
	{
		foreach my $bin (get_bins($region->[0],$region->[1],$repeats_ref->{bin_spacing}))
		{
			foreach my $repeat (@{$repeats_ref->{binned_repeats}{$chromosome}{$bin}})
			{
				if (overlap($repeat,$region))
				{
					push @repeat_list, $repeat->[2];
				}
			}
		}
	}
	
	return @repeat_list;
}
sub calculate_repeat_proportion_feature
{
	my $fusion_info_ref = shift;
	
	my %repeats;
	read_repeats($repeats_regions, \%repeats);
	
	foreach my $cluster_id (keys %{$fusion_info_ref})
	{
		foreach my $cluster_end ("0","1")
		{
			my $chromosome = $fusion_info_ref->{$cluster_id}{$cluster_end}{chromosome};
			my $regions = $fusion_info_ref->{$cluster_id}{$cluster_end}{genomic_regions};

			$fusion_info_ref->{$cluster_id}{$cluster_end}{repeat_proportion} = get_repeat_proportion($chromosome, $regions, \%repeats);
			$fusion_info_ref->{$cluster_id}{$cluster_end}{repeat_list} = join ",", get_repeat_list($chromosome, $regions, \%repeats);
		}
	}
}

