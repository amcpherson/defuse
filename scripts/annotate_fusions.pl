#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];
use Cwd qw[abs_path];

$| = 1;

use FindBin;
use lib "$FindBin::RealBin";
use configdata;
use cmdrunner;
use gene_models;

use lib "$FindBin::RealBin/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --name      Library Name\n";

my $help;
my $config_filename;
my $dataset_directory;
my $output_directory;
my $library_name;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'output=s'    => \$output_directory,
	'name=s'      => \$library_name,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $dataset_directory or die @usage;
defined $output_directory or die @usage;
defined $library_name or die @usage;

my $source_directory = abs_path("$FindBin::RealBin/../");

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

# Config values
my $reference_fasta			= $config->get_value("reference_fasta");
my $gene_models_filename	= $config->get_value("gene_models");
my $genome_fasta			= $config->get_value("genome_fasta");
my $cdna_fasta				= $config->get_value("cdna_fasta");
my $exons_fasta				= $config->get_value("exons_fasta");
my $cds_fasta				= $config->get_value("cds_fasta");
my $est_fasta				= $config->get_value("est_fasta");
my $est_alignments			= $config->get_value("est_alignments");
my $repeats_regions			= $config->get_value("repeats_regions");
my $splice_bias				= $config->get_value("splice_bias");
my $tools_directory			= $config->get_value("tools_directory");
my $scripts_directory		= $config->get_value("scripts_directory");
my $samtools_bin			= $config->get_value("samtools_bin");
my $percident_threshold		= $config->get_value("percent_identity_threshold");
my $calc_extra_anno			= $config->get_value("calculate_extra_annotations");

# Get samtools major version (0.1.x or 1.x)
my $samtools_version_major	= $1 if `$samtools_bin 2>&1` =~ /^Version:\s(\d+\.\d+).*$/m;

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

sub get_splice_seq
{
	my $chromosome = shift;
	my $position = shift;
	my $strand = shift;
	
	my $seq;
	if ($strand eq "+")
	{
		$seq = $genome_db->seq($chromosome, $position + 1, $position + 2);
	}
	elsif ($strand eq "-")
	{
		$seq = revcomp($genome_db->seq($chromosome, $position - 2, $position - 1));
	}
	defined $seq or die "Error: unable to obtain sequence at $chromosome:$position\n";
	
	return $seq;
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

my $breakpoints_genome_psl = $output_directory."/breakpoints.genome.psl";
my $breakpoints_genome_nointron_psl = $output_directory."/breakpoints.genome.nointron.psl";
my $breakpoints_cdna_psl = $output_directory."/breakpoints.cdna.psl";
my $breakpoints_est_psl = $output_directory."/breakpoints.est.psl";
my $breakpoints_exons_psl = $output_directory."/breakpoints.exons.psl";
my $breakpoints_cds_psl = $output_directory."/breakpoints.cds.psl";

my $break_filename = $output_directory."/splitreads.break";
my $seq_filename = $output_directory."/splitreads.seq";

my %break;
my %seq;

read_breaks($break_filename, \%break);
read_seq($seq_filename, \%seq);

# Read in the read stats
my $read_stats = $output_directory."/concordant.read.stats";
my %read_stat_values;
get_stats($read_stats, \%read_stat_values);

# Read in the gene models
my $gene_models = gene_models->new($gene_models_filename);

# Read in expression
my %expression;
my $expression_filename = $output_directory."/expression.txt";
read_expression($expression_filename, \%expression);

# Read in repeats list
my %repeats;
read_repeats($repeats_regions, \%repeats);

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
my %genomic_strand1;
my %genomic_strand2;
my %genomic_starts1;
my %genomic_starts2;
my %genomic_ends1;
my %genomic_ends2;
my %gene_location1;
my %gene_location2;
my %break_adj_entropy1;
my %break_adj_entropy2;
my %fusion_seq_length;
my %fusion_seq1_length;
my %fusion_seq2_length;
my %fusion_repeat_proportion1;
my %fusion_repeat_proportion2;
my %fusion_splice_score;
my %fusion_splice_variants;

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

sub get_start_end
{
	my $regions = shift;
	my $idx = shift;

	my @positions;
	foreach my $region (@{$regions})
	{
		push @positions, $region->[$idx];
	}

	return \@positions;
}

my %clusters;
my $clusters_sc = $output_directory."/clusters.sc";
read_clusters($clusters_sc, \%clusters);

foreach my $cluster_id (keys %clusters)
{
	if (not defined $break{$cluster_id})
	{
		delete $clusters{$cluster_id};
	}
}

# Calculate fusion align regions
my %fusion_region;
my %fusion_align_region;
foreach my $cluster_id (keys %clusters)
{
	foreach my $cluster_end (keys %{$clusters{$cluster_id}})
	{
		foreach my $fragment_id (keys %{$clusters{$cluster_id}{$cluster_end}})
		{
			my $align_start = $clusters{$cluster_id}{$cluster_end}{$fragment_id}{start};
			my $align_end = $clusters{$cluster_id}{$cluster_end}{$fragment_id}{end};
			
			$fusion_align_region{$cluster_id}{$cluster_end} = [$align_start, $align_end] if not defined $fusion_align_region{$cluster_id}{$cluster_end};
			
			$fusion_align_region{$cluster_id}{$cluster_end}->[0] = min($fusion_align_region{$cluster_id}{$cluster_end}->[0], $align_start);
			$fusion_align_region{$cluster_id}{$cluster_end}->[1] = max($fusion_align_region{$cluster_id}{$cluster_end}->[1], $align_end);
		}

		my $break_pos = $break{$cluster_id}{$cluster_end}{breakpos};
		my $strand = $break{$cluster_id}{$cluster_end}{strand};

		@{$fusion_region{$cluster_id}{$cluster_end}} = @{$fusion_align_region{$cluster_id}{$cluster_end}};

		if ($strand eq "+")
		{
			$fusion_region{$cluster_id}{$cluster_end}->[1] = $break_pos;
			$fusion_region{$cluster_id}{$cluster_end}->[0] = min($fusion_region{$cluster_id}{$cluster_end}->[0], $break_pos);
		}
		elsif ($strand eq "-")
		{
			$fusion_region{$cluster_id}{$cluster_end}->[0] = $break_pos;
			$fusion_region{$cluster_id}{$cluster_end}->[1] = max($fusion_region{$cluster_id}{$cluster_end}->[1], $break_pos);
		}
		else
		{
			die;
		}
	}
}

# Find fusion breakpoint information
my %cluster_ids;
foreach my $cluster_id (keys %break)
{
	my $ref_name1 = $break{$cluster_id}{"0"}{reference};
	my $ref_name2 = $break{$cluster_id}{"1"}{reference};
	
	my $strand1 = $break{$cluster_id}{"0"}{strand};
	my $strand2 = $break{$cluster_id}{"1"}{strand};
	
	my $break_pos1 = $break{$cluster_id}{"0"}{breakpos};
	my $break_pos2 = $break{$cluster_id}{"1"}{breakpos};
	
	my $gene1 = $gene_models->calc_gene($ref_name1, $break_pos1);
	my $gene2 = $gene_models->calc_gene($ref_name2, $break_pos2);
	
	my $genomic_breakpos1 = $gene_models->calc_genomic_position($ref_name1, $break_pos1);
	my $genomic_breakpos2 = $gene_models->calc_genomic_position($ref_name2, $break_pos2);
	
	my $genomic_strand1 = $gene_models->calc_genomic_strand($ref_name1, $strand1);
	my $genomic_strand2 = $gene_models->calc_genomic_strand($ref_name2, $strand2);
	
	my $gene_location1 = $gene_models->calc_gene_location($gene1, $genomic_breakpos1);
	my $gene_location2 = $gene_models->calc_gene_location($gene2, $genomic_breakpos2);
	
	my @genomic_align_regions1 = $gene_models->calc_genomic_regions($ref_name1, $fusion_align_region{$cluster_id}{"0"});
	my @genomic_align_regions2 = $gene_models->calc_genomic_regions($ref_name2, $fusion_align_region{$cluster_id}{"1"});

	my @genomic_regions1 = $gene_models->calc_genomic_regions($ref_name1, $fusion_region{$cluster_id}{"0"});
	my @genomic_regions2 = $gene_models->calc_genomic_regions($ref_name2, $fusion_region{$cluster_id}{"1"});
	
	my $chromosome1 = $gene_models->calc_genomic_chromosome($ref_name1);
	my $chromosome2 = $gene_models->calc_genomic_chromosome($ref_name2);
	
	$fusion_repeat_proportion1{$cluster_id} = get_repeat_proportion($chromosome1, \@genomic_align_regions1, \%repeats);
	$fusion_repeat_proportion2{$cluster_id} = get_repeat_proportion($chromosome2, \@genomic_align_regions2, \%repeats);
	
	my $gene_strand_a = ($gene1 lt $gene2) ? $gene1.$strand1 : $gene2.$strand2;
	my $gene_strand_b = ($gene1 lt $gene2) ? $gene2.$strand2 : $gene1.$strand1;
	my $breakpos_a = ($gene1 lt $gene2) ? $genomic_breakpos1 : $genomic_breakpos2;
	my $breakpos_b = ($gene1 lt $gene2) ? $genomic_breakpos2 : $genomic_breakpos1;

	$fusion_splice_variants{$gene_strand_a}{$gene_strand_b}{$breakpos_a."-".$breakpos_b} = 1;
	
	# Select the breakpoint prediction with the highest spanning pvalue
	my $break_sequence = $seq{$cluster_id}{sequence};
	
	# Skip if no breakpoint sequence was predicted
	next if $break_sequence eq "N" or $break_sequence eq "";
	
	# Find position of breakpoint in break sequence	
	my $break_in_seq = index $break_sequence, "|";
	
	# Skip if we dont find the break marker
	# TODO: fix this for denovo breakpoints
	next if $break_in_seq < 0;
	
	# Remove breakpoint marker
	$break_sequence =~ s/\|//g;
	
	$genomic_breakpos1{$cluster_id} = $genomic_breakpos1;
	$genomic_breakpos2{$cluster_id} = $genomic_breakpos2;
	
	$genomic_strand1{$cluster_id} = $genomic_strand1;
	$genomic_strand2{$cluster_id} = $genomic_strand2;

	$genomic_starts1{$cluster_id} = join ",", @{get_start_end(\@genomic_regions1, 0)};
	$genomic_starts2{$cluster_id} = join ",", @{get_start_end(\@genomic_regions2, 0)};

	$genomic_ends1{$cluster_id} = join ",", @{get_start_end(\@genomic_regions1, 1)};
	$genomic_ends2{$cluster_id} = join ",", @{get_start_end(\@genomic_regions2, 1)};

	$genomic_starts1{$cluster_id} = ($genomic_starts1{$cluster_id} eq "") ? "NA" : $genomic_starts1{$cluster_id};
	$genomic_starts2{$cluster_id} = ($genomic_starts2{$cluster_id} eq "") ? "NA" : $genomic_starts2{$cluster_id};

	$genomic_ends1{$cluster_id} = ($genomic_ends1{$cluster_id} eq "") ? "NA" : $genomic_ends1{$cluster_id};
	$genomic_ends2{$cluster_id} = ($genomic_ends2{$cluster_id} eq "") ? "NA" : $genomic_ends2{$cluster_id};

	$gene_location1{$cluster_id} = $gene_location1;
	$gene_location2{$cluster_id} = $gene_location2;
	
	$fusion_gene_lookup{$cluster_id}{$gene1} = 1;
	$fusion_gene_lookup{$cluster_id}{$gene2} = 1;

	$fusion_gene1{$cluster_id} = $gene1;
	$fusion_gene2{$cluster_id} = $gene2;

	$fusion_ref_name1{$cluster_id} = $ref_name1;
	$fusion_ref_name2{$cluster_id} = $ref_name2;

	$fusion_strand1{$cluster_id} = $strand1;
	$fusion_strand2{$cluster_id} = $strand2;
	
	my $break_adj_seq1 = substr $break_sequence, max(0, $break_in_seq - $entropy_breakpoint_adjacent_size), min($break_in_seq, $entropy_breakpoint_adjacent_size);
	my $break_adj_seq2 = substr $break_sequence, $break_in_seq, min(length($break_sequence) - $break_in_seq, $entropy_breakpoint_adjacent_size);
	
	$break_adj_entropy1{$cluster_id} = calcentropy($break_adj_seq1);
	$break_adj_entropy2{$cluster_id} = calcentropy($break_adj_seq2);
	
	$fusion_seq_length{$cluster_id} = length($break_sequence);
	$fusion_seq1_length{$cluster_id} = $break_in_seq;
	$fusion_seq2_length{$cluster_id} = length($break_sequence) - $break_in_seq;
	
	$cluster_ids{$cluster_id} = 1;
}
close SEQ;

my %fragment_ids;
my %fusion_fragments;
my %fusion_span_count;

# Find cluster fragment info
foreach my $cluster_id (keys %cluster_ids)
{
	foreach my $cluster_end (keys %{$clusters{$cluster_id}})
	{
		foreach my $fragment_id (keys %{$clusters{$cluster_id}{$cluster_end}})
		{
			$fragment_ids{$fragment_id} = 1;
			$fusion_fragments{$cluster_id}{$fragment_id} = 1;
		}
		
		$fusion_span_count{$cluster_id} = scalar keys %{$clusters{$cluster_id}{$cluster_end}};
	}
}

my $cdna_fasta_index = $cdna_fasta.".fai";
my $cdna_pair_sam = $output_directory."/cdna.pair.sam";
my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
my $cdna_pair_bam_bai = $cdna_pair_bam.".bai";
my $cdna_pair_bam_prefix = $cdna_pair_bam.".sort";
if ($calc_extra_anno eq "yes")
{
	my $samtools_sort_cmd = "$samtools_bin view -bt $cdna_fasta_index #<1 | $samtools_bin sort ".(($samtools_version_major < 1) ? "-o - $cdna_pair_bam_prefix > #>1" : "-T $cdna_pair_bam_prefix -O bam - > #>1");
	$runner->run($samtools_sort_cmd, [$cdna_pair_sam], [$cdna_pair_bam]);
	$runner->run("$samtools_bin index #<1", [$cdna_pair_bam], [$cdna_pair_bam_bai]);
}

# Calculate mapping statistics
my %mapping_stats;
calculate_mapping_stats(\%mapping_stats);

# Calculate break concordant counts
my %break_concordant;
if ($calc_extra_anno eq "yes")
{
	calculate_break_concordant($break_filename, \%break_concordant);
}

# Calculate interrupted
my %interruption_info;
if ($calc_extra_anno eq "yes")
{
	calculate_interrupted($break_filename, \%interruption_info);
}

sub calculate_splicing_index
{
	my $break_concordant_ref = shift;
	my $cluster_id = shift;
	my $cluster_end = shift;
	
	return "-" if not defined $break_concordant_ref->{$cluster_id} or not defined $break_concordant_ref->{$cluster_id}{$cluster_end};
	
	return $break_concordant_ref->{$cluster_id}{$cluster_end}{break_concordant} / $fusion_span_count{$cluster_id};
}

# Calculate splicing index
my %fusion_splicing_index1;
my %fusion_splicing_index2;
foreach my $cluster_id (keys %cluster_ids)
{
	$fusion_splicing_index1{$cluster_id} = calculate_splicing_index(\%break_concordant, $cluster_id, "0");
	$fusion_splicing_index2{$cluster_id} = calculate_splicing_index(\%break_concordant, $cluster_id, "1");
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
	$fusion_interrupted_index1{$cluster_id} = calculate_interrupted_index($interruption_info{$cluster_id}{"0"});
	$fusion_interrupted_index2{$cluster_id} = calculate_interrupted_index($interruption_info{$cluster_id}{"1"});
}

sub find_alignregion
{
	my $breakpoints_psl = shift;
	my $align_strand = shift;
	my $query_region = shift;
	my $target_region = shift;
	
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

print $target_seq_name."\n" if not defined $gene;	
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

find_alignregion($breakpoints_exons_psl, \%exon_align_strand, \%exon_query_region, \%exon_target_region);
find_alignregion($breakpoints_cds_psl, \%cds_align_strand, \%cds_query_region, \%cds_target_region);

my $minimum_coverage = $read_stat_values{"fraglength_mean"} - $read_stat_values{"readlength_min"};

my %span_coverage;
foreach my $cluster_id (keys %clusters)
{
	foreach my $cluster_end (keys %{$clusters{$cluster_id}})
	{
		my %covered;
		foreach my $fragment_id (keys %{$clusters{$cluster_id}{$cluster_end}})
		{
			foreach my $pos ($clusters{$cluster_id}{$cluster_end}{$fragment_id}{start} .. $clusters{$cluster_id}{$cluster_end}{$fragment_id}{end})
			{
				$covered{$pos} = 1;
			}
		}
		
		my $num_covered = scalar keys %covered;
		my $normalized_coverage = $num_covered / $minimum_coverage;
		
		$span_coverage{$cluster_id}{$cluster_end} = $normalized_coverage;
	}
}

my %breakseqs_percident;
find_breakseqs_percident($breakpoints_genome_psl, \%{$breakseqs_percident{genome}}, $genome_max_ins);
find_breakseqs_percident($breakpoints_cdna_psl, \%{$breakseqs_percident{cdna}}, $cdna_max_ins);
find_breakseqs_percident($breakpoints_est_psl, \%{$breakseqs_percident{est}}, $est_max_ins);
find_breakseqs_estislands_percident($breakpoints_genome_psl, \%{$breakseqs_percident{estisland}});

my %max_left_end;
my %min_right_start;
find_breakseqs_overlap($breakpoints_genome_nointron_psl, $percident_threshold, \%max_left_end, \%min_right_start);
find_breakseqs_overlap($breakpoints_cdna_psl, $percident_threshold, \%max_left_end, \%min_right_start);

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
	
	$breakpoint_homology{$cluster_id} = 0;
	if (defined $max_left_end{$cluster_id} and defined $min_right_start{$cluster_id})
	{
		$breakpoint_homology{$cluster_id} = max(0, $max_left_end{$cluster_id} - $min_right_start{$cluster_id} + 1);
	}
}

my %genome_breakseqs_percident;
find_breakseqs_percident($breakpoints_genome_psl, \%genome_breakseqs_percident);

my %altsplice;
foreach my $cluster_id (keys %cluster_ids)
{
	if (defined $genome_breakseqs_percident{$cluster_id} and $genome_breakseqs_percident{$cluster_id} > $percident_threshold)
	{
		$altsplice{$cluster_id} = "Y";
	}
	else
	{
		$altsplice{$cluster_id} = "N";
	}
}

foreach my $cluster_id (sort {$a <=> $b} keys %cluster_ids)
{
	my $gene1 = $fusion_gene1{$cluster_id};
	my $gene2 = $fusion_gene2{$cluster_id};
	
	my $ref_name1 = $fusion_ref_name1{$cluster_id};
	my $ref_name2 = $fusion_ref_name2{$cluster_id};
	
	my $genome_strand1 = $genomic_strand1{$cluster_id};
	my $genome_strand2 = $genomic_strand2{$cluster_id};

	my $transcript1 = (defined $gene_models->{transcripts}{$ref_name1}) ? $ref_name1 : 'NA';
	my $transcript2 = (defined $gene_models->{transcripts}{$ref_name2}) ? $ref_name2 : 'NA';

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
			
			my $target_phase;
			if ($strand1 eq "+")
			{
				$target_phase = ($target_region1->[0] - $target_region2->[0]) % 3;
			}
			else
			{
				$target_phase = ($target_region1->[1] - $target_region2->[1]) % 3;
			}

			$orf = "Y" if $query_phase == $target_phase;
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
					
					$genomic_breakpos1{$cluster_id} = $gene_models->exon_to_genome($target_name1, $query_end1_targetpos);
					$genomic_breakpos2{$cluster_id} = $gene_models->exon_to_genome($target_name2, $query_start2_targetpos);
					
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
					
					$genomic_breakpos1{$cluster_id} = $gene_models->exon_to_genome($target_name1, $query_start1_targetpos);
					$genomic_breakpos2{$cluster_id} = $gene_models->exon_to_genome($target_name2, $query_end2_targetpos);
					
					last;
				}
			}
		}
		
		last if $exonboundaries eq "Y";
	}

	my $splice_seq1 = get_splice_seq($gene_models->{genes}{$gene1}{chromosome}, $genomic_breakpos1{$cluster_id}, $genome_strand1);
	my $splice_seq2 = get_splice_seq($gene_models->{genes}{$gene2}{chromosome}, $genomic_breakpos2{$cluster_id}, $genome_strand2);

	my $splice_seqf = $splice_seq1.revcomp($splice_seq2);
	my $splice_seqr = $splice_seq2.revcomp($splice_seq1);

	my @splice_scores;
	push @splice_scores, calc_edit_dist("GTAG", $splice_seqf);
	push @splice_scores, calc_edit_dist("GTAG", $splice_seqr);
	push @splice_scores, calc_edit_dist("ATAC", $splice_seqf);
	push @splice_scores, calc_edit_dist("ATAC", $splice_seqr);
	
	$fusion_splice_score{$cluster_id} = 4 - min(@splice_scores);

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

	$breakseqs_percident{genome}{$cluster_id} = 0 if not defined $breakseqs_percident{genome}{$cluster_id};
	$breakseqs_percident{cdna}{$cluster_id} = 0 if not defined $breakseqs_percident{cdna}{$cluster_id};
	$breakseqs_percident{est}{$cluster_id} = 0 if not defined $breakseqs_percident{est}{$cluster_id};
	$breakseqs_percident{estisland}{$cluster_id} = 0 if not defined $breakseqs_percident{estisland}{$cluster_id};
	
	my $gene_strand_a = ($gene1 lt $gene2) ? $gene1.$fusion_strand1{$cluster_id} : $gene2.$fusion_strand2{$cluster_id};
	my $gene_strand_b = ($gene1 lt $gene2) ? $gene2.$fusion_strand2{$cluster_id} : $gene1.$fusion_strand1{$cluster_id};
	
	my $num_splice_variants = scalar keys %{$fusion_splice_variants{$gene_strand_a}{$gene_strand_b}};
	
	my $gene_align_strand1 = ($genome_strand1 eq $gene_models->{genes}{$gene1}{strand}) ? "+" : "-";
	my $gene_align_strand2 = ($genome_strand2 eq $gene_models->{genes}{$gene2}{strand}) ? "+" : "-";
	
	print $cluster_id."\tlibrary_name\t".$library_name."\n";

	print $cluster_id."\tgene1\t".$gene1."\n";
	print $cluster_id."\ttranscript1\t".$transcript1."\n";

	print $cluster_id."\tgene_name1\t".$gene_models->{genes}{$gene1}{name}."\n";
	print $cluster_id."\tgene_chromosome1\t".$gene_models->{genes}{$gene1}{chromosome}."\n";
	print $cluster_id."\tgene_strand1\t".$gene_models->{genes}{$gene1}{strand}."\n";
	print $cluster_id."\tgene_start1\t".$gene_models->{genes}{$gene1}{region}->[0]."\n";
	print $cluster_id."\tgene_end1\t".$gene_models->{genes}{$gene1}{region}->[1]."\n";

	print $cluster_id."\tgene2\t".$gene2."\n";
	print $cluster_id."\ttranscript2\t".$transcript2."\n";

	print $cluster_id."\tgene_name2\t".$gene_models->{genes}{$gene2}{name}."\n";
	print $cluster_id."\tgene_chromosome2\t".$gene_models->{genes}{$gene2}{chromosome}."\n";
	print $cluster_id."\tgene_strand2\t".$gene_models->{genes}{$gene2}{strand}."\n";
	print $cluster_id."\tgene_start2\t".$gene_models->{genes}{$gene2}{region}->[0]."\n";
	print $cluster_id."\tgene_end2\t".$gene_models->{genes}{$gene2}{region}->[1]."\n";
	
	print $cluster_id."\tgene_align_strand1\t".$gene_align_strand1."\n";
	print $cluster_id."\tgene_align_strand2\t".$gene_align_strand2."\n";
	
	print $cluster_id."\tgenomic_break_pos1\t".$genomic_breakpos1{$cluster_id}."\n";
	print $cluster_id."\tgenomic_break_pos2\t".$genomic_breakpos2{$cluster_id}."\n";	
	print $cluster_id."\tgenomic_strand1\t".$genome_strand1."\n";
	print $cluster_id."\tgenomic_strand2\t".$genome_strand2."\n";

	print $cluster_id."\tgenomic_starts1\t".$genomic_starts1{$cluster_id}."\n";	
	print $cluster_id."\tgenomic_starts2\t".$genomic_starts2{$cluster_id}."\n";	
	print $cluster_id."\tgenomic_ends1\t".$genomic_ends1{$cluster_id}."\n";	
	print $cluster_id."\tgenomic_ends2\t".$genomic_ends2{$cluster_id}."\n";	

	print $cluster_id."\tsplicing_index1\t".$fusion_splicing_index1{$cluster_id}."\n";
	print $cluster_id."\tsplicing_index2\t".$fusion_splicing_index2{$cluster_id}."\n";
	print $cluster_id."\tinterrupted_index1\t".$fusion_interrupted_index1{$cluster_id}."\n";
	print $cluster_id."\tinterrupted_index2\t".$fusion_interrupted_index2{$cluster_id}."\n";
	print $cluster_id."\tspan_coverage1\t".$span_coverage{$cluster_id}{"0"}."\n";
	print $cluster_id."\tspan_coverage2\t".$span_coverage{$cluster_id}{"1"}."\n";
	print $cluster_id."\texpression1\t".get_expression($gene1)."\n";
	print $cluster_id."\texpression2\t".get_expression($gene2)."\n";
	
	print $cluster_id."\tgene_location1\t".$gene_location1{$cluster_id}."\n";
	print $cluster_id."\tgene_location2\t".$gene_location2{$cluster_id}."\n";
	print $cluster_id."\torf\t".$orf."\n";
	print $cluster_id."\texonboundaries\t".$exonboundaries."\n";
	print $cluster_id."\tadjacent\t".$adjacent."\n";
	print $cluster_id."\tinterchromosomal\t".$interchromosomal."\n";
	print $cluster_id."\tinversion\t".$inversion."\n";
	print $cluster_id."\teversion\t".$eversion."\n";
	print $cluster_id."\tdeletion\t".$deletion."\n";
	print $cluster_id."\tread_through\t".$read_through."\n";
	print $cluster_id."\taltsplice\t".$altsplice{$cluster_id}."\n";
	
	print $cluster_id."\tspan_count\t".$fusion_span_count{$cluster_id}."\n";
	print $cluster_id."\tgenome_breakseqs_percident\t".$breakseqs_percident{genome}{$cluster_id}."\n";
	print $cluster_id."\tcdna_breakseqs_percident\t".$breakseqs_percident{cdna}{$cluster_id}."\n";
	print $cluster_id."\test_breakseqs_percident\t".$breakseqs_percident{est}{$cluster_id}."\n";
	print $cluster_id."\tbreakseqs_estislands_percident\t".$breakseqs_percident{estisland}{$cluster_id}."\n";

	print $cluster_id."\tbreak_adj_entropy1\t".$break_adj_entropy1{$cluster_id}."\n";
	print $cluster_id."\tbreak_adj_entropy2\t".$break_adj_entropy2{$cluster_id}."\n";
	print $cluster_id."\tbreakpoint_homology\t".$breakpoint_homology{$cluster_id}."\n";
	
	print $cluster_id."\tbreak_adj_entropy_min\t".min($break_adj_entropy1{$cluster_id},$break_adj_entropy2{$cluster_id})."\n";
	print $cluster_id."\tspan_coverage_min\t".min($span_coverage{$cluster_id}{"0"},$span_coverage{$cluster_id}{"1"})."\n";
	print $cluster_id."\tspan_coverage_max\t".max($span_coverage{$cluster_id}{"0"},$span_coverage{$cluster_id}{"1"})."\n";

	print $cluster_id."\trepeat_proportion1\t".$fusion_repeat_proportion1{$cluster_id}."\n";
	print $cluster_id."\trepeat_proportion2\t".$fusion_repeat_proportion2{$cluster_id}."\n";
	print $cluster_id."\tmax_repeat_proportion\t".max($fusion_repeat_proportion1{$cluster_id},$fusion_repeat_proportion2{$cluster_id})."\n";
	print $cluster_id."\tsplice_score\t".$fusion_splice_score{$cluster_id}."\n";
	print $cluster_id."\tnum_splice_variants\t".$num_splice_variants."\n";
	
	print $cluster_id."\tmin_map_count\t".$mapping_stats{$cluster_id}{min_map_count}."\n";
	print $cluster_id."\tmax_map_count\t".$mapping_stats{$cluster_id}{max_map_count}."\n";
	print $cluster_id."\tmean_map_count\t".$mapping_stats{$cluster_id}{mean_map_count}."\n";
	print $cluster_id."\tnum_multi_map\t".$mapping_stats{$cluster_id}{num_multi_map}."\n";
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

sub read_break_predict
{
	my $break_predict_filename = shift;
	my $break_predict_hash_ref = shift;
	
	open BP, $break_predict_filename or die "Error: Unable to find $break_predict_filename: $!\n";
	while (<BP>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		
		$break_predict_hash_ref->{$cluster_id}{break_predict} = $fields[1];
	}
	close BP;
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
	my $mapping_stats_ref = shift;
	
	my $mapping_stats_filename = $output_directory."/mapping.stats";
	$runner->run("$calc_map_stats_script -c $config_filename -d $dataset_directory -o $output_directory > #>1", [$clusters_sc], [$mapping_stats_filename]);
	
	open MPS, $mapping_stats_filename or die "Error: Unable to find $mapping_stats_filename: $!\n";
	while (<MPS>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $cluster_id = $fields[0];
		my $anno_type = $fields[1];
		my $anno_value = $fields[2];
		
		$mapping_stats_ref->{$cluster_id}{$anno_type} = $anno_value;
	}
	close MPS;
}

sub calculate_break_concordant
{
	my $breaks_filename = shift;
	my $break_concordant_ref = shift;

	my $break_concordant_filename = $breaks_filename.".concordant.counts";

	$runner->run("$break_concordant_script -c $config_filename -d $dataset_directory -o $output_directory -b #<1 > #>1", [$breaks_filename], [$break_concordant_filename]);

	open BRC, $break_concordant_filename or die "Error: Unable to open $break_concordant_filename: $!\n";
	while (<BRC>)
	{
		chomp;
		my ($cluster_id, $cluster_end, $break_concordant) = split /\t/;
		$break_concordant_ref->{$cluster_id}{$cluster_end}{break_concordant} = $break_concordant;
	}
	close BRC;
}

sub calculate_interrupted
{
	my $breaks_filename = shift;
	my $interrupted_ref = shift;

	my $interrupted_filename = $breaks_filename.".interrupted.counts";

	$runner->run("$interrupted_script -c $config_filename -d $dataset_directory -o $output_directory -b #<1 > #>1", [$breaks_filename], [$interrupted_filename]);

	open BIN, $interrupted_filename or die "Error: Unable to open $interrupted_filename: $!\n";
	while (<BIN>)
	{
		chomp;
		my ($cluster_id, $cluster_end, $gene, $size_before, $size_after, $count_before, $count_after) = split /\t/;
		$interrupted_ref->{$cluster_id}{$cluster_end}{gene} = $gene;
		$interrupted_ref->{$cluster_id}{$cluster_end}{size_before} = $size_before;
		$interrupted_ref->{$cluster_id}{$cluster_end}{size_after} = $size_after;
		$interrupted_ref->{$cluster_id}{$cluster_end}{count_before} = $count_before;
		$interrupted_ref->{$cluster_id}{$cluster_end}{count_after} = $count_after;
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


