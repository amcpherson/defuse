#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use List::Util qw[min max];

use lib dirname($0);
use configdata;
use cmdrunner;

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
my $genome_fasta			= $config->get_value("genome_fasta");
my $gene_fasta				= $config->get_value("gene_fasta");
my $cdna_fasta				= $config->get_value("cdna_fasta");
my $rrna_fasta				= $config->get_value("rrna_fasta");
my $est_fasta				= $config->get_value("est_fasta");
my $est_alignments			= $config->get_value("est_alignments");
my $exons_fasta				= $config->get_value("exons_fasta");
my $exons_regions			= $config->get_value("exons_regions");
my $cds_fasta				= $config->get_value("cds_fasta");
my $utr5p_fasta				= $config->get_value("utr5p_fasta");
my $utr3p_fasta				= $config->get_value("utr3p_fasta");
my $tss_regions				= $config->get_value("tss_regions");
my $tts_regions				= $config->get_value("tts_regions");
my $intronic_fasta			= $config->get_value("intronic_fasta");
my $upstream_fasta			= $config->get_value("upstream_fasta");
my $downstream_fasta		= $config->get_value("downstream_fasta");
my $gene_info_list			= $config->get_value("gene_info_list");
my $gene_adj_list			= $config->get_value("gene_adj_list");
my $cdna_gene_fasta			= $config->get_value("cdna_gene_fasta");
my $cdna_gene_regions		= $config->get_value("cdna_gene_regions");
my $cdna_regions			= $config->get_value("cdna_regions");
my $gene_tran_list			= $config->get_value("gene_tran_list");
my $splice_bias             = $config->get_value("splice_bias");
my $denovo_assembly			= $config->get_value("denovo_assembly");
my $tools_directory			= $config->get_value("tools_directory");
my $scripts_directory		= $config->get_value("scripts_directory");
my $blat_bin				= $config->get_value("blat_bin");
my $samtools_bin			= $config->get_value("samtools_bin");

my $genome_max_ins = 2000;
my $est_max_ins = 10000;
my $cdna_max_ins = 10000000;

my $entropy_breakpoint_adjacent_size = 40;

my $estislandsbin = $tools_directory."/estislands";
my $break_concordant_script = $scripts_directory."/calc_break_concordant.pl";
my $interrupted_script = $scripts_directory."/calc_interrupted.pl";

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/annotate";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("annotate");
$runner->prefix($log_prefix);

# Read gene info
my %gene_info;
read_gene_info($gene_info_list, \%gene_info);

my $splitr_break_filename = $output_directory."/splitr.break";
my $denovo_break_filename = $output_directory."/denovo.break";
my $splitr_seq_filename = $output_directory."/splitr.seq";
my $denovo_seq_filename = $output_directory."/denovo.seq";
my $splitr_span_pval_filename = $output_directory."/splitr.span.pval";
my $denovo_span_pval_filename = $output_directory."/denovo.span.pval";

my %splitr_break;
my %denovo_break;
my %splitr_seq;
my %denovo_seq;
my %splitr_span_pval;
my %denovo_span_pval;

read_breaks($splitr_break_filename, \%splitr_break);
read_breaks($denovo_break_filename, \%denovo_break);
read_splitr_seq($splitr_seq_filename, \%splitr_seq);
read_denovo_seq($denovo_seq_filename, \%denovo_seq);
read_span_pval($splitr_span_pval_filename, \%splitr_span_pval);
read_span_pval($denovo_span_pval_filename, \%denovo_span_pval);

# Read in the read stats
my $read_stats = $output_directory."/concordant.read.stats";
my %read_stat_values;
get_stats($read_stats, \%read_stat_values);

# Read in the cdna gene regions
my %cdna_gene_reg;
read_regions($cdna_gene_regions, \%cdna_gene_reg);

# Read in the exons regions
my %exons_reg;
read_regions($exons_regions, \%exons_reg);

# Read in expression
my %expression;
my $expression_filename = $output_directory."/expression.txt";
read_expression($expression_filename, \%expression);

sub get_expression
{
	my $gene = shift;
	return 0 if not defined $expression{$gene};
	return $expression{$gene}{expression};
}

my %fusion_gene_lookup;
my %fusion_gene1;
my %fusion_gene2;
my %fusion_ref_name1;
my %fusion_ref_name2;
my %fusion_strand1;
my %fusion_strand2;
my %genomic_breakpos1;
my %genomic_breakpos2;
my %fusion_break_predict;
my %break_adj_entropy1;
my %break_adj_entropy2;
my %fusion_seq_length;
my %fusion_seq1_length;
my %fusion_seq2_length;

sub calcentropy
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

my $breakpoint_sequences_temp_filename = $output_directory."/breakpoint.sequences.fa.tmp";
my $breakpoint_sequences_filename = $output_directory."/breakpoint.sequences.fa";

# Find fusion breakpoint information
# Create breakpoint sequence fasta
open SEQ, ">".$breakpoint_sequences_temp_filename or die "Error: Unable to open file $breakpoint_sequences_temp_filename\n";
my %cluster_ids = (%splitr_seq, %denovo_seq);
foreach my $cluster_id (keys %cluster_ids)
{
	my $ref_name1 = $splitr_break{$cluster_id}{breakpos}->[0]->[0];
	my $ref_name2 = $splitr_break{$cluster_id}{breakpos}->[1]->[0];
	
	$ref_name1 =~ /(ENSG\d+)/ or die "Error: Unable to interpret $ref_name1\n"; my $gene1 = $1;
	$ref_name2 =~ /(ENSG\d+)/ or die "Error: Unable to interpret $ref_name2\n"; my $gene2 = $1;

	my $strand1 = $splitr_break{$cluster_id}{breakpos}->[0]->[1];
	my $strand2 = $splitr_break{$cluster_id}{breakpos}->[1]->[1];
	
	$genomic_breakpos1{$cluster_id} = calc_genomic_position($splitr_break{$cluster_id}{breakpos}->[0]->[2], $cdna_gene_reg{$ref_name1});
	$genomic_breakpos2{$cluster_id} = calc_genomic_position($splitr_break{$cluster_id}{breakpos}->[1]->[2], $cdna_gene_reg{$ref_name2});
	
	# Select the breakpoint prediction with the highest spanning pvalue
	my $break_sequence = "";
	my $break_predict = "";
	if (lc($denovo_assembly) eq "yes" and $denovo_span_pval{$cluster_id}{pvalue} > $splitr_span_pval{$cluster_id}{pvalue})
	{
		$break_sequence = $denovo_seq{$cluster_id}{sequence};
		$break_predict = "denovo";
	}
	else
	{
		$break_sequence = $splitr_seq{$cluster_id}{sequence};
		$break_predict = "splitr";
	}
	
	# Skip if no breakpoint sequence was predicted
	delete $cluster_ids{$cluster_id} and next if $break_sequence eq "";

	# Find position of breakpoint in break sequence	
	my $break_in_seq = index $break_sequence, "|";

	# Skip if we dont find the break marker
	next if $break_in_seq < 0;
	
	# Remove breakpoint marker
	$break_sequence =~ s/\|//g;
	
	die "Error: Unable to find break infor for $cluster_id\n" if scalar @{$splitr_break{$cluster_id}{breakpos}} != 2;
	
	print SEQ ">$cluster_id\n".$break_sequence."\n";
	
	$fusion_gene_lookup{$cluster_id}{$gene1} = 1;
	$fusion_gene_lookup{$cluster_id}{$gene2} = 1;

	$fusion_gene1{$cluster_id} = $gene1;
	$fusion_gene2{$cluster_id} = $gene2;

	$fusion_ref_name1{$cluster_id} = $ref_name1;
	$fusion_ref_name2{$cluster_id} = $ref_name2;

	$fusion_strand1{$cluster_id} = $strand1;
	$fusion_strand2{$cluster_id} = $strand2;
	
	$fusion_break_predict{$cluster_id} = $break_predict;

	my $break_adj_seq1 = substr $break_sequence, max(0, $break_in_seq - $entropy_breakpoint_adjacent_size), min($break_in_seq, $entropy_breakpoint_adjacent_size);
	my $break_adj_seq2 = substr $break_sequence, $break_in_seq, min(length($break_sequence) - $break_in_seq, $entropy_breakpoint_adjacent_size);
	
	$break_adj_entropy1{$cluster_id} = calcentropy($break_adj_seq1);
	$break_adj_entropy2{$cluster_id} = calcentropy($break_adj_seq2);
	
	$fusion_seq_length{$cluster_id} = length($break_sequence);
	$fusion_seq1_length{$cluster_id} = $break_in_seq;
	$fusion_seq2_length{$cluster_id} = length($break_sequence) - $break_in_seq;
}
close SEQ;

# Replace old file if we generated a different one
$runner->replaceifdifferent($breakpoint_sequences_temp_filename, $breakpoint_sequences_filename);

my %clusters;
my $clusters_sam = $output_directory."/clusters.sc.sam";
read_sam_clusters($clusters_sam, \%clusters);

my $read_sequences_temp_filename = $output_directory."/read.sequences.fa.tmp";
my $read_sequences_filename = $output_directory."/read.sequences.fa";

my %fragment_ids;
my %fusion_fragments;
my %fusion_span_count;

# Find cluster fragment info
# Create a read sequence fasta
open RFA, ">".$read_sequences_temp_filename or die "Error: Unable to open file $read_sequences_temp_filename\n";
foreach my $cluster_id (keys %cluster_ids)
{
	die "Error: fusion $cluster_id does not refer to 2 reference sequences\n" if scalar keys %{$clusters{$cluster_id}} != 2;
	foreach my $ref_name (keys %{$clusters{$cluster_id}})
	{
		foreach my $fragment_id (keys %{$clusters{$cluster_id}{$ref_name}})
		{
			print RFA ">".$clusters{$cluster_id}{$ref_name}{$fragment_id}{read_id}."\n";
			print RFA $clusters{$cluster_id}{$ref_name}{$fragment_id}{sequence}."\n";
					        
			$fragment_ids{$fragment_id} = 1;
			$fusion_fragments{$cluster_id}{$fragment_id} = 1;
		}
		
		$fusion_span_count{$cluster_id} = scalar keys %{$clusters{$cluster_id}{$ref_name}};
	}
}
close RFA;

# Replace old file if we generated a different one
$runner->replaceifdifferent($read_sequences_temp_filename, $read_sequences_filename);

# Calculate break concordant counts
my %splitr_break_concordant;
my %denovo_break_concordant;
calculate_break_concordant($splitr_break_filename, \%splitr_break_concordant);
calculate_break_concordant($denovo_break_filename, \%denovo_break_concordant);

# Calculate interrupted
my %splitr_interruption_info;
my %denovo_interruption_info;
calculate_interrupted($splitr_break_filename, \%splitr_interruption_info);
calculate_interrupted($denovo_break_filename, \%denovo_interruption_info);

# Calculate splicing index
my %fusion_splicing_index1;
my %fusion_splicing_index2;
foreach my $cluster_id (keys %cluster_ids)
{
	if ($fusion_break_predict{$cluster_id} eq "splitr")
	{
		$fusion_splicing_index1{$cluster_id} = $splitr_break_concordant{$cluster_id}{$fusion_gene1{$cluster_id}}{break_concordant} / $fusion_span_count{$cluster_id};
		$fusion_splicing_index2{$cluster_id} = $splitr_break_concordant{$cluster_id}{$fusion_gene2{$cluster_id}}{break_concordant} / $fusion_span_count{$cluster_id};
	}
	elsif ($fusion_break_predict{$cluster_id} eq "denovo")
	{
		$fusion_splicing_index1{$cluster_id} = $denovo_break_concordant{$cluster_id}{$fusion_gene1{$cluster_id}}{break_concordant} / $fusion_span_count{$cluster_id};
		$fusion_splicing_index2{$cluster_id} = $denovo_break_concordant{$cluster_id}{$fusion_gene2{$cluster_id}}{break_concordant} / $fusion_span_count{$cluster_id};
	}
}

sub calculate_interrupted_index
{
	my $interrupted_ref = shift;

	return "-" if not defined $interrupted_ref->{count_before};
	
	my $expression_before = $interrupted_ref->{count_before} / ($interrupted_ref->{size_before} + 1) + 1;
	my $expression_after = $interrupted_ref->{count_after} / ($interrupted_ref->{size_after} + 1) + 1;
	
	return $expression_before / $expression_after;
}

# Calculate interrupted index
my %fusion_interrupted_index1;
my %fusion_interrupted_index2;
foreach my $cluster_id (keys %cluster_ids)
{
	if ($fusion_break_predict{$cluster_id} eq "splitr")
	{
		$fusion_interrupted_index1{$cluster_id} = calculate_interrupted_index($splitr_interruption_info{$cluster_id}{$fusion_gene1{$cluster_id}});
		$fusion_interrupted_index2{$cluster_id} = calculate_interrupted_index($splitr_interruption_info{$cluster_id}{$fusion_gene2{$cluster_id}});
	}
	elsif ($fusion_break_predict{$cluster_id} eq "denovo")
	{
		$fusion_interrupted_index1{$cluster_id} = calculate_interrupted_index($splitr_interruption_info{$cluster_id}{$fusion_gene1{$cluster_id}});
		$fusion_interrupted_index2{$cluster_id} = calculate_interrupted_index($splitr_interruption_info{$cluster_id}{$fusion_gene2{$cluster_id}});
	}
}

my %gene_tss;
open TSS, $tss_regions or die "Error: Unable to open $tss_regions\n";
while (<TSS>)
{
	chomp;
	my @fields = split /\t/;
	
	my $gene_id = $fields[1];
	my $tss_start = $fields[3];
	my $tss_end = $fields[4];
	
	push @{$gene_tss{$gene_id}}, [$tss_start, $tss_end];
}
close TSS;

my %gene_tts;
open TTS, $tts_regions or die "Error: Unable to open $tts_regions\n";
while (<TTS>)
{
	chomp;
	my @fields = split /\t/;
	
	my $gene_id = $fields[1];
	my $tts_start = $fields[3];
	my $tts_end = $fields[4];
	
	push @{$gene_tts{$gene_id}}, [$tts_start, $tts_end];
}
close TTS;

my %gene_adjacency;
open GAJ, $gene_adj_list or die "Error: Unable to open $gene_adj_list\n";
while (<GAJ>)
{
	chomp;
	my @fields = split /\t/;
	
	my $gene_id1 = $fields[0];
	my $gene_id2 = $fields[1];

	$gene_adjacency{$gene_id1}{$gene_id2} = 1;
}
close GAJ;

sub find_alignregion
{
	my $sequences_fasta = $_[0];
	my $blat_out = $_[1];
	my $align_strand = $_[2];
	my $query_region = $_[3];
	my $target_region = $_[4];
	
	$runner->run("$blat_bin -noHead $sequences_fasta #<1 #>1", [$breakpoint_sequences_filename], [$blat_out]);

	open PSL, $blat_out or die "Error: Unable to open $blat_out\n";
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
	
		$target_seq_name =~ /(ENSG\d+)/;
		my $gene = $1;
	
		next unless $fusion_gene_lookup{$cluster_id}{$gene};
		
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
			
			push @{$align_strand->{$cluster_id}{$gene}}, $strand;
			push @{$query_region->{$cluster_id}{$gene}}, [$query_start, $query_end, $query_size];
			push @{$target_region->{$cluster_id}{$gene}}, [$target_start, $target_end, $target_size, $target_seq_name];
		}
	}
	close PSL;
}

sub find_overlap
{
	my $sequences_fasta = $_[0];
	my $blat_out = $_[1];
	my $overlapping = $_[2];

	$runner->run("$blat_bin -noHead $sequences_fasta #<1 #>1", [$breakpoint_sequences_filename], [$blat_out]);

	open PSL, $blat_out or die "Error: Unable to open $blat_out\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
		
		my $num_matches = $psl_fields[0];
		my $num_target_bases_inserted = $psl_fields[7];
		my $strand = $psl_fields[8];
		my $cluster_id = $psl_fields[9];
		my $breakpoint_seq_length = $psl_fields[10];
		my $target_seq_name = $psl_fields[13];
		my $query_starts = $psl_fields[19];
		my $target_starts = $psl_fields[20];
		
		$target_seq_name =~ /(ENSG\d+)/;
		my $gene = $1;
		
		next unless $fusion_gene_lookup{$cluster_id}{$gene};
		
		$overlapping->{$cluster_id}{$gene} = 1;
	}
	close PSL;
}

sub contains
{
	my $region = $_[0];
	my $testregion = $_[1];
	
	if ($region->[0] <= $testregion->[0] and $region->[1] >= $testregion->[1])
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

sub contains_array
{
	my $region = $_[0];
	my $testregions = $_[1];
	
	foreach my $testregion (@{$testregions})
	{
		if (contains($region, $testregion))
		{
			return 1;
		}
	}
	
	return 0;
}

my %exon_align_strand;
my %exon_query_region;
my %exon_target_region;

my %cds_align_strand;
my %cds_query_region;
my %cds_target_region;

my %gene_align_strand;
my %gene_query_region;
my %gene_target_region;

my $exon_align_output = $output_directory."/exon.psl";
my $cds_align_output = $output_directory."/cds.psl";
my $gene_align_output = $output_directory."/gene.psl";

find_alignregion($exons_fasta, $exon_align_output, \%exon_align_strand, \%exon_query_region, \%exon_target_region);
find_alignregion($cds_fasta, $cds_align_output, \%cds_align_strand, \%cds_query_region, \%cds_target_region);
find_alignregion($gene_fasta, $gene_align_output, \%gene_align_strand, \%gene_query_region, \%gene_target_region);

my $cdna_overlap_output = $output_directory."/cdna_overlap.psl";
my $cds_overlap_output = $output_directory."/cds_overlap.psl";
my $utr5p_overlap_output = $output_directory."/utr5p_overlap.psl";
my $utr3p_overlap_output = $output_directory."/utr3p_overlap.psl";
my $intronic_overlap_output = $output_directory."/intronic_overlap.psl";
my $upstream_overlap_output = $output_directory."/upstream_overlap.psl";
my $downstream_overlap_output = $output_directory."/downstream_overlap.psl";

my %cdna_overlap;
my %cds_overlap;
my %utr5p_overlap;
my %utr3p_overlap;
my %intronic_overlap;
my %upstream_overlap;
my %downstream_overlap;

find_overlap($cdna_fasta, $cdna_overlap_output, \%cdna_overlap);
find_overlap($cds_fasta, $cds_overlap_output, \%cds_overlap);
find_overlap($utr5p_fasta, $utr5p_overlap_output, \%utr5p_overlap);
find_overlap($utr3p_fasta, $utr3p_overlap_output, \%utr3p_overlap);
find_overlap($intronic_fasta, $intronic_overlap_output, \%intronic_overlap);
find_overlap($upstream_fasta, $upstream_overlap_output, \%upstream_overlap);
find_overlap($downstream_fasta, $downstream_overlap_output, \%downstream_overlap);

my %concordant_reads;
find_concordant_reads($genome_fasta, $read_sequences_filename, \%concordant_reads, $genome_max_ins);
find_concordant_reads($cdna_fasta, $read_sequences_filename, \%concordant_reads, $cdna_max_ins);
find_anchored_concordant_reads($rrna_fasta, $read_sequences_filename, \%concordant_reads);

my %concordant_ratio;
foreach my $cluster_id (keys %fusion_span_count)
{
	my $span_count = $fusion_span_count{$cluster_id};
	
	my $concordant_fragments = 0;
	foreach my $fragment_id (keys %{$fusion_fragments{$cluster_id}})
	{
		$concordant_fragments++ if $concordant_reads{$fragment_id};
	}

	$concordant_ratio{$cluster_id} = $concordant_fragments / $span_count;
}

my $minimum_coverage = $read_stat_values{"fraglength_mean"} - $read_stat_values{"readlength_min"};

my %span_coverage;
foreach my $cluster_id (keys %clusters)
{
	foreach my $ref_name (keys %{$clusters{$cluster_id}})
	{
		my %covered;
		foreach my $fragment_id (keys %{$clusters{$cluster_id}{$ref_name}})
		{
			foreach my $pos ($clusters{$cluster_id}{$ref_name}{$fragment_id}{start} .. $clusters{$cluster_id}{$ref_name}{$fragment_id}{end})
			{
				$covered{$pos} = 1;
			}
		}
		
		my $num_covered = scalar keys %covered;
		my $normalized_coverage = $num_covered / $minimum_coverage;
		
		$span_coverage{$cluster_id}{$ref_name} = $normalized_coverage;
	}
}

my %breakseqs_percident;
find_breakseqs_percident($genome_fasta, $breakpoint_sequences_filename, \%{$breakseqs_percident{genome}}, $genome_max_ins);
find_breakseqs_percident($cdna_fasta, $breakpoint_sequences_filename, \%{$breakseqs_percident{cdna}}, $cdna_max_ins);
find_breakseqs_percident($est_fasta, $breakpoint_sequences_filename, \%{$breakseqs_percident{est}}, $est_max_ins);
find_breakseqs_estislands_percident($genome_fasta, $breakpoint_sequences_filename, \%{$breakseqs_percident{estisland}});

my %breakpoint_matches;
find_breakseqs_matches($cdna_gene_fasta, $breakpoint_sequences_filename, \%breakpoint_matches);

# Scale all the percent identities
# Calculate breakpoint homologies
my %breakpoint_homology;
foreach my $cluster_id (keys %cluster_ids)
{
	foreach my $ref_type ("genome", "cdna", "est", "estisland")
	{
		next if not defined $breakseqs_percident{$ref_type}{$cluster_id};
		
		my $mismatches = (1 - $breakseqs_percident{$ref_type}{$cluster_id}) * $fusion_seq_length{$cluster_id};
		my $adjusted_percident = 1 - ($mismatches / min($fusion_seq1_length{$cluster_id}, $fusion_seq2_length{$cluster_id}));
		$adjusted_percident = max(0,$adjusted_percident);
	
		$breakseqs_percident{$ref_type}{$cluster_id} = $adjusted_percident;
	}
	
	$breakpoint_matches{$cluster_id}{$fusion_gene1{$cluster_id}} = $fusion_seq1_length{$cluster_id} if not defined $breakpoint_matches{$cluster_id}{$fusion_gene1{$cluster_id}};
	$breakpoint_matches{$cluster_id}{$fusion_gene2{$cluster_id}} = $fusion_seq2_length{$cluster_id} if not defined $breakpoint_matches{$cluster_id}{$fusion_gene2{$cluster_id}};

	my $matches1 = $breakpoint_matches{$cluster_id}{$fusion_gene1{$cluster_id}} - $fusion_seq1_length{$cluster_id};
	my $matches2 = $breakpoint_matches{$cluster_id}{$fusion_gene2{$cluster_id}} - $fusion_seq2_length{$cluster_id};
	
	$breakpoint_homology{$cluster_id} = $matches1 + $matches2;
}

foreach my $cluster_id (sort {$a <=> $b} keys %cluster_ids)
{
	my $gene1 = $fusion_gene1{$cluster_id};
	my $gene2 = $fusion_gene2{$cluster_id};
	
	my $ref_name1 = $fusion_ref_name1{$cluster_id};
	my $ref_name2 = $fusion_ref_name2{$cluster_id};

	my $upstream1 = "N";
	$upstream1 = "Y" if defined $upstream_overlap{$cluster_id}{$gene1};

	my $upstream2 = "N";
	$upstream2 = "Y" if defined $upstream_overlap{$cluster_id}{$gene2};

	my $exonic1 = "N";
	$exonic1 = "Y" if defined $cdna_overlap{$cluster_id}{$gene1};

	my $exonic2 = "N";
	$exonic2 = "Y" if defined $cdna_overlap{$cluster_id}{$gene2};

	my $intronic1 = "N";
	$intronic1 = "Y" if defined $intronic_overlap{$cluster_id}{$gene1};

	my $intronic2 = "N";
	$intronic2 = "Y" if defined $intronic_overlap{$cluster_id}{$gene2};

	my $downstream1 = "N";
	$downstream1 = "Y" if defined $downstream_overlap{$cluster_id}{$gene1};

	my $downstream2 = "N";
	$downstream2 = "Y" if defined $downstream_overlap{$cluster_id}{$gene2};

	my $utr5p1 = "N";
	$utr5p1 = "Y" if defined $utr5p_overlap{$cluster_id}{$gene1};

	my $utr5p2 = "N";
	$utr5p2 = "Y" if defined $utr5p_overlap{$cluster_id}{$gene2};

	my $coding1 = "N";
	$coding1 = "Y" if defined $cds_overlap{$cluster_id}{$gene1};

	my $coding2 = "N";
	$coding2 = "Y" if defined $cds_overlap{$cluster_id}{$gene2};

	my $utr3p1 = "N";
	$utr3p1 = "Y" if defined $utr3p_overlap{$cluster_id}{$gene1};

	my $utr3p2 = "N";
	$utr3p2 = "Y" if defined $utr3p_overlap{$cluster_id}{$gene2};

	my $orf = "N";
	foreach my $start_index1 (0..$#{$cds_align_strand{$cluster_id}{$gene1}})
	{
		my $strand1 = $cds_align_strand{$cluster_id}{$gene1}->[$start_index1];
		my $query_region1 = $cds_query_region{$cluster_id}{$gene1}->[$start_index1];
		my $target_region1 = $cds_target_region{$cluster_id}{$gene1}->[$start_index1];
	
		foreach my $start_index2 (0..$#{$cds_align_strand{$cluster_id}{$gene2}})
		{
			my $strand2 = $cds_align_strand{$cluster_id}{$gene2}->[$start_index2];
			my $query_region2 = $cds_query_region{$cluster_id}{$gene2}->[$start_index2];
			my $target_region2 = $cds_target_region{$cluster_id}{$gene2}->[$start_index2];

			next unless $strand1 eq $strand2;

			my $query_phase = ($query_region1->[0] - $query_region2->[0]) % 3;
			my $target_phase = ($target_region1->[0] - $target_region2->[0]) % 3;

			$orf = "Y" if $query_phase == $target_phase;
		}
	}

	my $utr5pexchange = "N";
	my $utr3pexchange = "N";
	foreach my $start_index1 (0..$#{$gene_align_strand{$cluster_id}{$gene1}})
	{
		my $strand1 = $gene_align_strand{$cluster_id}{$gene1}->[$start_index1];
		my $query_region1 = $gene_query_region{$cluster_id}{$gene1}->[$start_index1];
		my $target_region1 = $gene_target_region{$cluster_id}{$gene1}->[$start_index1];

		foreach my $start_index2 (0..$#{$gene_align_strand{$cluster_id}{$gene2}})
		{
			my $strand2 = $gene_align_strand{$cluster_id}{$gene2}->[$start_index2];
			my $query_region2 = $gene_query_region{$cluster_id}{$gene2}->[$start_index2];
			my $target_region2 = $gene_target_region{$cluster_id}{$gene2}->[$start_index2];

			next unless $strand1 eq $strand2;

			$utr5pexchange = "Y" if $utr5p1 eq "Y" and $utr5p2 eq "Y";
			$utr3pexchange = "Y" if $utr3p1 eq "Y" and $utr3p2 eq "Y";
		}
	}

	my $exonboundaries = "N";
	foreach my $start_index1 (0..$#{$exon_align_strand{$cluster_id}{$gene1}})
	{
		my $strand1 = $exon_align_strand{$cluster_id}{$gene1}->[$start_index1];
		my $query_region1 = $exon_query_region{$cluster_id}{$gene1}->[$start_index1];
		my $target_region1 = $exon_target_region{$cluster_id}{$gene1}->[$start_index1];

		my $target_size1 = $target_region1->[2];
		my $target_name1 = $target_region1->[3];

		foreach my $start_index2 (0..$#{$exon_align_strand{$cluster_id}{$gene2}})
		{	
			my $strand2 = $exon_align_strand{$cluster_id}{$gene2}->[$start_index2];
			my $query_region2 = $exon_query_region{$cluster_id}{$gene2}->[$start_index2];
			my $target_region2 = $exon_target_region{$cluster_id}{$gene2}->[$start_index2];

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
					$exonboundaries = "Y";
					
					$genomic_breakpos1{$cluster_id} = calc_genomic_position($query_end1_targetpos, $exons_reg{$target_name1});
					$genomic_breakpos2{$cluster_id} = calc_genomic_position($query_start2_targetpos, $exons_reg{$target_name2});
					
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
					$exonboundaries = "Y";
					
					$genomic_breakpos1{$cluster_id} = calc_genomic_position($query_start1_targetpos, $exons_reg{$target_name1});
					$genomic_breakpos2{$cluster_id} = calc_genomic_position($query_end2_targetpos, $exons_reg{$target_name2});
					
					last;
				}
			}
		}
		
		last if $exonboundaries eq "Y";
	}

	my $adjacent = "N";
	$adjacent = "Y" if defined $gene_adjacency{$gene1}{$gene2};

	my $genome_strand1 = ($gene_info{$gene1}{strand} eq $fusion_strand1{$cluster_id}) ? "+" : "-";
	my $genome_strand2 = ($gene_info{$gene2}{strand} eq $fusion_strand2{$cluster_id}) ? "+" : "-";

	my $interchromosomal = "N";
	$interchromosomal = "Y" if $gene_info{$gene1}{chromosome} ne $gene_info{$gene2}{chromosome};

	my $inversion = "N";
	$inversion = "Y" if $interchromosomal eq "N" and $genome_strand1 eq $genome_strand2;

	my $eversion = "N";
	$eversion = "Y" if $interchromosomal eq "N" and $gene_info{$gene1}{start} < $gene_info{$gene2}{start} and $genome_strand1 eq "-" and $genome_strand2 eq "+";
	$eversion = "Y" if $interchromosomal eq "N" and $gene_info{$gene1}{start} > $gene_info{$gene2}{start} and $genome_strand1 eq "+" and $genome_strand2 eq "-";
	
	my $deletion = "N";
	$deletion = "Y" if $interchromosomal eq "N" and $inversion eq "N" and $eversion eq "N";

	my $read_through = "N";
	$read_through = "Y" if $deletion eq "Y" and $adjacent eq "Y";

	$breakseqs_percident{genome}{$cluster_id} = 0 if not defined $breakseqs_percident{genome}{$cluster_id};
	$breakseqs_percident{cdna}{$cluster_id} = 0 if not defined $breakseqs_percident{cdna}{$cluster_id};
	$breakseqs_percident{est}{$cluster_id} = 0 if not defined $breakseqs_percident{est}{$cluster_id};
	$breakseqs_percident{estisland}{$cluster_id} = 0 if not defined $breakseqs_percident{estisland}{$cluster_id};

	print $cluster_id."\tlibrary_name\t".$library_name."\n";

	print $cluster_id."\tgene1\t".$gene1."\n";
	print $cluster_id."\tgene_name1\t".$gene_info{$gene1}{name}."\n";
	print $cluster_id."\tgene_chromosome1\t".$gene_info{$gene1}{chromosome}."\n";
	print $cluster_id."\tgene_strand1\t".$gene_info{$gene1}{strand}."\n";
	print $cluster_id."\tgene_start1\t".$gene_info{$gene1}{start}."\n";
	print $cluster_id."\tgene_end1\t".$gene_info{$gene1}{end}."\n";

	print $cluster_id."\tgene2\t".$gene2."\n";
	print $cluster_id."\tgene_name2\t".$gene_info{$gene2}{name}."\n";
	print $cluster_id."\tgene_chromosome2\t".$gene_info{$gene2}{chromosome}."\n";
	print $cluster_id."\tgene_strand2\t".$gene_info{$gene2}{strand}."\n";
	print $cluster_id."\tgene_start2\t".$gene_info{$gene2}{start}."\n";
	print $cluster_id."\tgene_end2\t".$gene_info{$gene2}{end}."\n";
	
	print $cluster_id."\tgene_align_strand1\t".$fusion_strand1{$cluster_id}."\n";
	print $cluster_id."\tgene_align_strand2\t".$fusion_strand2{$cluster_id}."\n";
	
	print $cluster_id."\tgenomic_break_pos1\t".$genomic_breakpos1{$cluster_id}."\n";
	print $cluster_id."\tgenomic_break_pos2\t".$genomic_breakpos2{$cluster_id}."\n";	
	print $cluster_id."\tgenomic_strand1\t".$genome_strand1."\n";
	print $cluster_id."\tgenomic_strand2\t".$genome_strand2."\n";	

	print $cluster_id."\tsplicing_index1\t".$fusion_splicing_index1{$cluster_id}."\n";
	print $cluster_id."\tsplicing_index2\t".$fusion_splicing_index2{$cluster_id}."\n";
	print $cluster_id."\tinterrupted_index1\t".$fusion_interrupted_index1{$cluster_id}."\n";
	print $cluster_id."\tinterrupted_index2\t".$fusion_interrupted_index2{$cluster_id}."\n";
	print $cluster_id."\tspan_coverage1\t".$span_coverage{$cluster_id}{$ref_name1}."\n";
	print $cluster_id."\tspan_coverage2\t".$span_coverage{$cluster_id}{$ref_name2}."\n";
	print $cluster_id."\texpression1\t".get_expression($gene1)."\n";
	print $cluster_id."\texpression2\t".get_expression($gene2)."\n";
	
	print $cluster_id."\tupstream1\t".$upstream1."\n";
	print $cluster_id."\tupstream2\t".$upstream2."\n";
	print $cluster_id."\texonic1\t".$exonic1."\n";
	print $cluster_id."\texonic2\t".$exonic2."\n";
	print $cluster_id."\tintronic1\t".$intronic1."\n";
	print $cluster_id."\tintronic2\t".$intronic2."\n";
	print $cluster_id."\tdownstream1\t".$downstream1."\n";
	print $cluster_id."\tdownstream2\t".$downstream2."\n";
	print $cluster_id."\tutr5p1\t".$utr5p1."\n";
	print $cluster_id."\tutr5p2\t".$utr5p2."\n";
	print $cluster_id."\tcoding1\t".$coding1."\n";
	print $cluster_id."\tcoding2\t".$coding2."\n";
	print $cluster_id."\tutr3p1\t".$utr3p1."\n";
	print $cluster_id."\tutr3p2\t".$utr3p2."\n";
	print $cluster_id."\torf\t".$orf."\n";
	print $cluster_id."\tutr5pexchange\t".$utr5pexchange."\n";
	print $cluster_id."\tutr3pexchange\t".$utr3pexchange."\n";
	print $cluster_id."\texonboundaries\t".$exonboundaries."\n";
	print $cluster_id."\tadjacent\t".$adjacent."\n";
	print $cluster_id."\tinterchromosomal\t".$interchromosomal."\n";
	print $cluster_id."\tinversion\t".$inversion."\n";
	print $cluster_id."\teversion\t".$eversion."\n";
	print $cluster_id."\tdeletion\t".$deletion."\n";
	print $cluster_id."\tread_through\t".$read_through."\n";
	
	print $cluster_id."\tspan_count\t".$fusion_span_count{$cluster_id}."\n";
	print $cluster_id."\tconcordant_ratio\t".$concordant_ratio{$cluster_id}."\n";
	print $cluster_id."\tgenome_breakseqs_percident\t".$breakseqs_percident{genome}{$cluster_id}."\n";
	print $cluster_id."\tcdna_breakseqs_percident\t".$breakseqs_percident{cdna}{$cluster_id}."\n";
	print $cluster_id."\test_breakseqs_percident\t".$breakseqs_percident{est}{$cluster_id}."\n";
	print $cluster_id."\tbreakseqs_estislands_percident\t".$breakseqs_percident{estisland}{$cluster_id}."\n";
	print $cluster_id."\tbreak_predict\t".$fusion_break_predict{$cluster_id}."\n";

	print $cluster_id."\tbreak_adj_entropy1\t".$break_adj_entropy1{$cluster_id}."\n";
	print $cluster_id."\tbreak_adj_entropy2\t".$break_adj_entropy2{$cluster_id}."\n";
	print $cluster_id."\tbreakpoint_homology\t".$breakpoint_homology{$cluster_id}."\n";
	
	print $cluster_id."\tbreak_adj_entropy_min\t".min($break_adj_entropy1{$cluster_id},$break_adj_entropy2{$cluster_id})."\n";
	print $cluster_id."\tspan_coverage_min\t".min($span_coverage{$cluster_id}{$ref_name1},$span_coverage{$cluster_id}{$ref_name2})."\n";
	print $cluster_id."\tspan_coverage_max\t".max($span_coverage{$cluster_id}{$ref_name1},$span_coverage{$cluster_id}{$ref_name2})."\n";
}

sub read_splitr_seq
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
		$seqs_hash_ref->{$cluster_id}{inter_length} = $fields[2];
		$seqs_hash_ref->{$cluster_id}{split_count} = $fields[3];
		$seqs_hash_ref->{$cluster_id}{split_pos_average} = $fields[4];
		$seqs_hash_ref->{$cluster_id}{split_min_average} = $fields[5];
	}
	close SEQ;
}

sub read_denovo_seq
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
		$seqs_hash_ref->{$cluster_id}{inter_length} = $fields[2];
		$seqs_hash_ref->{$cluster_id}{min_count} = $fields[3];
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
		my $reference = $fields[1];
		my $strand = $fields[2];
		my $breakpos = $fields[3];
		
		push @{$breaks_hash_ref->{$cluster_id}{breakpos}}, [$reference,$strand,$breakpos];
	}
	close BR;
}

sub read_span_pval
{
	my $span_pval_filename = shift;
	my $span_pval_hash_ref = shift;
	
	open SPP, $span_pval_filename or die "Error: Unable to find $span_pval_filename: $!\n";
	while (<SPP>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		
		$span_pval_hash_ref->{$cluster_id}{pvalue} = $fields[2];
	}
	close SPP;
}

sub read_sam_clusters
{
	my $clusters_filename = shift;
	my $clusters_hash_ref = shift;
	
	open CLU, $clusters_filename or die "Error: Unable to find $clusters_filename: $!\n";
	while (<CLU>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $read_id = $fields[1];
		my $flag = $fields[2];
		my $ref_name = $fields[3];
		my $start = $fields[4];
		my $end = $start + length($fields[10]) - 1;
		my $sequence = $fields[10];
	
		$read_id =~ /(.*)\/([12])/;
		my $fragment_id = $1;
		my $read_end = $2;
	
		my $strand;
		if ($flag & hex('0x0010'))
		{
			$strand = "-";
		}
		else
		{
			$strand = "+";
		}
		
		if ($strand eq "-")
		{
			$sequence = reverse($sequence);
			$sequence =~ tr/ACGTacgt/TGCAtgca/;
		}
	
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{read_id} = $read_id;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{read_end} = $read_end;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{strand} = $strand;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{sequence} = $sequence;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{start} = $start;
		$clusters_hash_ref->{$cluster_id}{$ref_name}{$fragment_id}{end} = $end;
	}
	close CLU;
}

sub find_concordant_reads
{
	my $reference_fasta = shift;
	my $reads_fasta = shift;
	my $concordant_ref = shift;
	my $max_insert_size = shift;
	
	my $reads_psl = $reads_fasta.".".basename($reference_fasta).".psl";

	$runner->run("$blat_bin -noHead $reference_fasta #<1 #>1", [$reads_fasta], [$reads_psl]);
	
	my %read_align;
		
	open PSL, $reads_psl or die "Error: Unable to open $reads_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
		
		my $matches = $psl_fields[0];
		my $read_id = $psl_fields[9];
		my $qsize = $psl_fields[10];
		my $align_chr = $psl_fields[13];
		my $align_start = $psl_fields[15];
		my $align_end = $psl_fields[16];
		my $align_strand = $psl_fields[8];
	
		$read_id =~ /(.*)\/([12])/;
		my $fragment_id = $1;
		my $read_end = $2;

		my $percent_identity = $matches / $qsize;
		next if $percent_identity < 0.90;

		$align_chr =~ s/ENST\d+//g;

		push @{$read_align{$fragment_id}{$read_end}}, [$align_chr,$align_strand,$align_start,$align_end];
	}
	close PSL;
	
	foreach my $fragment_id (keys %fragment_ids)
	{
		foreach my $aligninfo1 (@{$read_align{$fragment_id}{'1'}})
		{
			foreach my $aligninfo2 (@{$read_align{$fragment_id}{'2'}})
			{
				next unless $aligninfo1->[0] eq $aligninfo2->[0];

				if (abs($aligninfo2->[2] - $aligninfo1->[3]) < $max_insert_size)
				{
					$concordant_ref->{$fragment_id} = 1;
				}
				elsif (abs($aligninfo1->[2] - $aligninfo2->[3]) < $max_insert_size)
				{
					$concordant_ref->{$fragment_id} = 1;
				}
				
				last if defined $concordant_ref->{$fragment_id};
			}
			
			last if defined $concordant_ref->{$fragment_id};
		}
	}
}

sub find_multimapped_reads
{
	my $reference_fasta = shift;
	my $reads_fasta = shift;
	my $multimap_ref = shift;
	
	my $reads_psl = $reads_fasta.".".basename($reference_fasta).".psl";

	$runner->run("$blat_bin -noHead $reference_fasta #<1 #>1", [$reads_fasta], [$reads_psl]);
	
	my %read_align;
		
	open PSL, $reads_psl or die "Error: Unable to open $reads_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
		
		my $matches = $psl_fields[0];
		my $read_id = $psl_fields[9];
		my $qsize = $psl_fields[10];
		my $align_chr = $psl_fields[13];
		my $align_start = $psl_fields[15];
		my $align_end = $psl_fields[16];
		my $align_strand = $psl_fields[8];
	
		$read_id =~ /(.*)\/([12])/;
		my $fragment_id = $1;
		my $read_end = $2;

		my $percent_identity = $matches / $qsize;
		next if $percent_identity < 0.90;

		push @{$read_align{$fragment_id}{$read_end}}, [$align_chr,$align_strand,$align_start,$align_end];
	}
	close PSL;

	foreach my $fragment_id (keys %read_align)
	{
		foreach my $read_end (keys %{$read_align{$fragment_id}})
		{
			$multimap_ref->{$fragment_id}{$read_end} = scalar @{$read_align{$fragment_id}{$read_end}};
		}
	}
}

sub alignoverlap
{
	my $aligninfo1 = shift;
	my $aligninfo2 = shift;

	if ($aligninfo1->[0] ne $aligninfo2->[0])
	{
		return 0;
	}

	if ($aligninfo1->[2] > $aligninfo2->[3] or $aligninfo2->[2] > $aligninfo1->[3])
	{
		return 0;
	}

	return 1;
}

sub find_anchored_concordant_reads
{
	my $reference_fasta = shift;
	my $reads_fasta = shift;
	my $concordant_ref = shift;

	my $reads_psl = $reads_fasta.".".basename($reference_fasta).".psl";

	$runner->run("$blat_bin -noHead $reference_fasta #<1 #>1", [$reads_fasta], [$reads_psl]);
		
	open PSL, $reads_psl or die "Error: Unable to open $reads_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;

		my $matches = $psl_fields[0];
		my $qsize = $psl_fields[10];
		my $read_id = $psl_fields[9];

		my $percent_identity = $matches / $qsize;	
		next if $percent_identity < 0.9;

		$read_id =~ /(.*)\/([12])/;
		my $fragment_id = $1;
		my $read_end = $2;

		$concordant_ref->{$fragment_id} = 1;
	}
	close PSL;
}

sub find_breakseqs_percident
{
	my $reference_fasta = shift;
	my $breakpoints_fasta = shift;
	my $max_percident_ref = shift;
	my $max_ins = shift;

	my $breakpoints_psl = $breakpoints_fasta.".".basename($reference_fasta).".psl";
	
	$runner->run("$blat_bin -noHead $reference_fasta #<1 #>1", [$breakpoints_fasta], [$breakpoints_psl]);

	open PSL, $breakpoints_psl or die "Error: Unable to open $breakpoints_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
	
		my $num_matches = $psl_fields[0];
		my $num_target_bases_inserted = $psl_fields[7];
		my $cluster_id = $psl_fields[9];
		my $breakpoint_seq_length = $psl_fields[10];
	
		next if $num_target_bases_inserted > $max_ins;
	
		my $percent_identity = $num_matches / $breakpoint_seq_length;
		
		$max_percident_ref->{$cluster_id} = 0 if not defined $max_percident_ref->{$cluster_id};
		$max_percident_ref->{$cluster_id} = max($max_percident_ref->{$cluster_id}, $percent_identity);
	}
	close PSL;
}

sub find_breakseqs_matches
{
	my $reference_fasta = shift;
	my $breakpoints_fasta = shift;
	my $max_matches_ref = shift;

	my $breakpoints_psl = $breakpoints_fasta.".".basename($reference_fasta).".psl";
	
	$runner->run("$blat_bin -noHead $reference_fasta #<1 #>1", [$breakpoints_fasta], [$breakpoints_psl]);

	open PSL, $breakpoints_psl or die "Error: Unable to open $breakpoints_psl\n";
	while (<PSL>)
	{
		chomp;
		my @psl_fields = split /\t/;
	
		my $num_matches = $psl_fields[0];
		my $num_target_bases_inserted = $psl_fields[7];
		my $cluster_id = $psl_fields[9];
		my $breakpoint_seq_length = $psl_fields[10];
		my $target_seq_name = $psl_fields[13];

		$target_seq_name =~ /(ENSG\d+)/;
		my $gene = $1;

		next unless $fusion_gene_lookup{$cluster_id}{$gene};

		$max_matches_ref->{$cluster_id}{$gene} = 0 if not defined $max_matches_ref->{$cluster_id}{$gene};
		$max_matches_ref->{$cluster_id}{$gene} = max($max_matches_ref->{$cluster_id}{$gene}, $num_matches);
	}
	close PSL;
}

sub find_breakseqs_estislands_percident
{
	my $genome_fasta = shift;
	my $breakpoints_fasta = shift;
	my $max_percident_ref = shift;
	
	my $genome_psl = $breakpoints_fasta.".".basename($genome_fasta).".psl";
	my $estislands_psl = $breakpoints_fasta.".".basename($genome_fasta).".estisl.psl";
	
	$runner->run("$blat_bin -noHead $genome_fasta #<1 #>1", [$breakpoints_fasta], [$genome_psl]);
	$runner->run("$estislandsbin -b #<1 -e $est_alignments -o #>1", [$genome_psl], [$estislands_psl]);
	
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

sub read_gene_info
{
	my $gene_info_list = shift;
	my $gene_info_ref = shift;
	
	# Read in gene information
	open GI, $gene_info_list or die "Error: Unable to open $gene_info_list: $!\n";
	while (<GI>)
	{
		chomp;
		my @fields = split /\t/;
	
		my $ensgene = $fields[0];

		$gene_info_ref->{$ensgene}{name} = $fields[1];
		$gene_info_ref->{$ensgene}{chromosome} = $fields[2];
		$gene_info_ref->{$ensgene}{strand} = $fields[3];
		$gene_info_ref->{$ensgene}{start} = $fields[4];
		$gene_info_ref->{$ensgene}{end} = $fields[5];
	}
	close GI;
}

sub read_gene_transcript
{
	my $gene_tran_filename = shift;
	my $gene_tran_ref = shift;
	
	# Read in gene transcript mapping
	open GT, $gene_tran_filename or die "Error: Unable to open $gene_tran_filename: $!\n";
	while (<GT>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $ensgene = $fields[0];
		my $enstran = $fields[1];
		
		push @{$gene_tran_ref->{$ensgene}}, $enstran;
	}
	close GT;
}

sub read_regions
{
	my $regions_filename = shift;
	my $regions_hash_ref = shift;

	open REG, $regions_filename or die;
	while (<REG>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $gene = $fields[0];
		my $chromosome = $fields[1];
		my $strand = $fields[2];
	
		my @exons;
		my $fieldindex = 4;
		while ($fieldindex <= $#fields)
		{
			push @exons, [$fields[$fieldindex-1],$fields[$fieldindex]];
			$fieldindex += 2;
		}
		
		$regions_hash_ref->{$gene}{chromosome} = $chromosome;
		$regions_hash_ref->{$gene}{strand} = $strand;
		$regions_hash_ref->{$gene}{exons} = [@exons];
	}
	close REG;
}

sub calculate_break_concordant
{
	my $breaks_filename = shift;
	my $break_concordant_ref = shift;

	my $break_concordant_filename = $breaks_filename.".concordant.counts";

	$runner->run("$break_concordant_script -c $config_filename -o $output_directory -b #<1 > #>1", [$breaks_filename], [$break_concordant_filename]);

	open BRC, $break_concordant_filename or die "Error: Unable to open $break_concordant_filename: $!\n";
	while (<BRC>)
	{
		chomp;
		my ($cluster_id, $gene, $break_concordant) = split /\t/;
		$break_concordant_ref->{$cluster_id}{$gene}{break_concordant} = $break_concordant;
	}
	close BRC;
}

sub calculate_interrupted
{
	my $breaks_filename = shift;
	my $interrupted_ref = shift;

	my $interrupted_filename = $breaks_filename.".interrupted.counts";

	$runner->run("$interrupted_script -c $config_filename -o $output_directory -b #<1 > #>1", [$breaks_filename], [$interrupted_filename]);

	open BIN, $interrupted_filename or die "Error: Unable to open $interrupted_filename: $!\n";
	while (<BIN>)
	{
		chomp;
		my ($cluster_id, $gene, $size_before, $size_after, $count_before, $count_after) = split /\t/;
		$interrupted_ref->{$cluster_id}{$gene}{size_before} = $size_before;
		$interrupted_ref->{$cluster_id}{$gene}{size_after} = $size_after;
		$interrupted_ref->{$cluster_id}{$gene}{count_before} = $count_before;
		$interrupted_ref->{$cluster_id}{$gene}{count_after} = $count_after;
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

# Find position in genome given a position and the strand and exons of the transcript
sub calc_genomic_position
{
	my $position = shift;
	my $transcript_ref = shift;
	
	my $strand = $transcript_ref->{strand};
	my $exons = $transcript_ref->{exons};
	
	if ($strand eq "-")
	{
		$position = regions_length(@{$exons}) - $position + 1;
	}
	
	if ($position < 1)
	{
		return $exons->[0]->[0] + $position - 1;
	}
	
	my $local_offset = 0;
	foreach my $exon (@{$exons})
	{
		my $exonsize = $exon->[1] - $exon->[0] + 1;
			
		if ($position <= $local_offset + $exonsize)
		{
			return $position - $local_offset - 1 + $exon->[0];
		}
				
		$local_offset += $exonsize;
	}
	
	return $position - $local_offset + $exons->[$#{$exons}]->[1];
}

# Find position in a transcript given a genomic position and strand and exons of the transcript
# This version returns the position of the beginning of the next exon if the genomic position is intronic
sub calc_transcript_position
{
	my $position = shift;
	my $transcript_ref = shift;

	my $strand = $transcript_ref->{strand};
	my $exons = $transcript_ref->{exons};
	
	my $local_offset = 0;
	my $transcript_position;
	foreach my $exon (@{$exons})
	{
		my $exonsize = $exon->[1] - $exon->[0] + 1;

		if ($position <= $exon->[1])
		{
			if ($position < $exon->[0])
			{
				$transcript_position = $local_offset + 1;
				last;
			}
			else
			{
				$transcript_position = $local_offset + $position - $exon->[0] + 1;
				last;
			}
		}
				
		$local_offset += $exonsize;
	}

	$transcript_position = regions_length(@{$exons}) if not defined $transcript_position;
	
	if ($strand eq "-")
	{
		$transcript_position = regions_length(@{$exons}) - $transcript_position + 1;
	}
	
	return $transcript_position;
}


