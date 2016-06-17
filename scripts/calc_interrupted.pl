#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];
use Cwd qw[abs_path];

use FindBin;
use lib "$FindBin::RealBin";
use configdata;
use gene_models;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -b, --breaks    Breaks filename\n";

my $help;
my $config_filename;
my $dataset_directory;
my $output_directory;
my $breaks_filename;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'output=s'    => \$output_directory,
	'breaks=s'    => \$breaks_filename,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $dataset_directory or die @usage;
defined $output_directory or die @usage;
defined $breaks_filename or die @usage;

my $source_directory = abs_path("$FindBin::RealBin/../");

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

# Config values
my $gene_models_filename	= $config->get_value("gene_models");
my $splice_bias             = $config->get_value("splice_bias");
my $samtools_bin			= $config->get_value("samtools_bin");

# Require concordant alignments bam
my $cdna_bam_filename = $output_directory."/cdna.pair.bam";

# Read in the read stats
my $read_stats = $output_directory."/concordant.read.stats";
my %read_stat_values;
get_stats($read_stats, \%read_stat_values);

# Approximate max fragment length
my $max_fragment_length = int($read_stat_values{"fraglength_mean"} + 3 * $read_stat_values{"fraglength_stddev"});

# Read in the gene models
my $gene_models = gene_models->new($gene_models_filename);

# Read the breakpoints file
my %breaks;
read_breaks($breaks_filename, \%breaks);

my %transcript_fusion_position;
my %fusion_gene;
my %fusion_strand;
my %exons_before_size;
my %exons_after_size;
foreach my $cluster_id (keys %breaks)
{
	foreach my $cluster_end ("0","1")
	{
		my $reference = $breaks{$cluster_id}{$cluster_end}{reference};
		my $strand = $breaks{$cluster_id}{$cluster_end}{strand};
		my $breakpos = $breaks{$cluster_id}{$cluster_end}{breakpos};
		
		my $gene_id = $gene_models->calc_gene($reference, $breakpos);
		my $gene_location = $gene_models->calc_gene_location($gene_id, $breakpos);
		
		next if $gene_location eq "upstream" or $gene_location eq "downstream";
		
		$fusion_gene{$cluster_id}{$cluster_end} = $gene_id;
		$fusion_strand{$cluster_id}{$cluster_end} = $gene_models->{genes}{$gene_id}{strand};
		
		# Find the genomic position of the break with bias
		my $breakpos_genomic;
		if ($strand eq "+")
		{
			$breakpos_genomic = $gene_models->calc_genomic_position($reference, $breakpos - $splice_bias) + $splice_bias;
		}
		elsif ($strand eq "-")
		{
			$breakpos_genomic = $gene_models->calc_genomic_position($reference, $breakpos + $splice_bias) - $splice_bias;
		}
			
		# Find exons that occur before and after the break
		# Find position of break in each cdna
		my @exonsbefore;
		my @exonsafter;
		foreach my $transcript_id (keys %{$gene_models->{genes}{$gene_id}{transcripts}})
		{
			my @exons = @{$gene_models->{transcripts}{$transcript_id}{exons}};
			
			# Find position of break in cdna space
			my $breakpos_transcript = $gene_models->calc_transcript_position($transcript_id, $breakpos_genomic);
			die if not defined $breakpos_transcript or $breakpos_transcript < 0;
				
			# Add the cdna pos for this transcript and fusion
			$transcript_fusion_position{$transcript_id}{$cluster_id}{$cluster_end} = $breakpos_transcript;
	
			# Find exons before and after
			foreach my $exon (@exons)
			{
				if ($exon->[1] < $breakpos_genomic)
				{
					push @exonsbefore, $exon;
				}
				elsif ($exon->[0] > $breakpos_genomic)
				{
					push @exonsafter, $exon;
				}
				else
				{
					push @exonsbefore, [$exon->[0],$breakpos_genomic];
					push @exonsafter, [$breakpos_genomic,$exon->[1]];
				}
			}
		}
		
		my $exons_before_size = 0;
		my $exons_after_size = 0;
		$exons_before_size = regions_length(merge_regions(@exonsbefore)) if scalar @exonsbefore > 0;
		$exons_after_size = regions_length(merge_regions(@exonsafter)) if scalar @exonsafter > 0;
		
		if ($gene_models->{genes}{$gene_id}{strand} eq "-")
		{
			($exons_before_size,$exons_after_size) = ($exons_after_size,$exons_before_size);
		}
		
		$exons_before_size{$cluster_id}{$cluster_end} = $exons_before_size;
		$exons_after_size{$cluster_id}{$cluster_end} = $exons_after_size;
	}
}

my %counts_before;
my %counts_after;
foreach my $transcript_id (keys %transcript_fusion_position)
{
	my $gene_id = $gene_models->{transcripts}{$transcript_id}{gene};
	
	open TA, "$samtools_bin view $cdna_bam_filename '$transcript_id' |" or die "Error: Unable to run samtools on $cdna_bam_filename: $!\n";
	while (<TA>)
	{
		my @sam_fields = split /\t/;
			
		my $rname = $sam_fields[2];
		my $pos = $sam_fields[3];
		my $seq = $sam_fields[9];
			
		if ($transcript_id ne $rname)
		{
			warn "Error: samtools retrieived alignments to $rname when alignments to $transcript_id were requested\n";
			next;
		}
		
		my $start = $pos;
		my $end = $pos + length($seq) - 1;

		foreach my $cluster_id (keys %{$transcript_fusion_position{$transcript_id}})
		{
			foreach my $cluster_end (keys %{$transcript_fusion_position{$transcript_id}{$cluster_id}})
			{
				my $fusion_cdnapos = $transcript_fusion_position{$transcript_id}{$cluster_id}{$cluster_end};
	
				if ($end < $fusion_cdnapos)
				{
					$counts_before{$cluster_id}{$cluster_end} += length($seq);
				}
				elsif ($start > $fusion_cdnapos)
				{
					$counts_after{$cluster_id}{$cluster_end} += length($seq);
				}
				else
				{
					$counts_before{$cluster_id}{$cluster_end} += $fusion_cdnapos - $start;
					$counts_after{$cluster_id}{$cluster_end} += $end - $fusion_cdnapos;				
				}
			}
		}
	}
	close TA;
	($? >> 8) == 0 or die "Error: Unable to run samtools on $cdna_bam_filename\n";
}

foreach my $cluster_id (keys %breaks)
{
	foreach my $cluster_end (keys %{$breaks{$cluster_id}})
	{
		my $gene_id = $fusion_gene{$cluster_id}{$cluster_end};
		my $strand = $fusion_strand{$cluster_id}{$cluster_end};
		
		next if not defined $gene_id;

		my $count_before = $counts_before{$cluster_id}{$cluster_end};
		my $count_after = $counts_after{$cluster_id}{$cluster_end};
		
		my $size_before = $exons_before_size{$cluster_id}{$cluster_end};
		my $size_after = $exons_after_size{$cluster_id}{$cluster_end};
		
		if ($strand eq "-")
		{
			($count_before,$count_after) = ($count_after,$count_before);
			($size_before,$size_after) = ($size_after,$size_before);
		}

		$count_before = 0 if not defined $count_before;
		$count_after = 0 if not defined $count_after;
		
		print $cluster_id."\t".$cluster_end."\t".$gene_id."\t".$size_before."\t".$size_after."\t".$count_before."\t".$count_after."\n";
	}
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

