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

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -b, --breaks    Breaks filename\n";
push @usage, "  -g, --genspce   Breaks in Genomic Space\n";

my $help;
my $config_filename;
my $output_directory;
my $breaks_filename;
my $genomic_space;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
	'breaks=s'    => \$breaks_filename,
	'genspce'     => \$genomic_space,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $output_directory or die @usage;
defined $breaks_filename or die @usage;

my $config = configdata->new();
$config->read($config_filename);

# Config values
my $cdna_gene_regions		= $config->get_value("cdna_gene_regions");
my $cdna_regions			= $config->get_value("cdna_regions");
my $gene_tran_list			= $config->get_value("gene_tran_list");
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

# Read in the gene transcripts
my %gene_transcripts;
read_gene_transcript($gene_tran_list, \%gene_transcripts);

# Read in the cdna regions
my %cdna;
read_regions($cdna_regions, \%cdna);

# Read in the cdna gene regions
my %cdna_gene;
read_regions($cdna_gene_regions, \%cdna_gene);

# Read the breakpoints file
my %breaks;
read_breaks($breaks_filename, \%breaks);

my %transcript_fusion_position;
my %fusion_gene_exons_before_size;
my %fusion_gene_exons_after_size;
foreach my $cluster_id (keys %breaks)
{
	foreach my $cluster_end ("0","1")
	{
		my $reference = $breaks{$cluster_id}{$cluster_end}{reference};
		my $strand = $breaks{$cluster_id}{$cluster_end}{strand};
		my $position = $breaks{$cluster_id}{$cluster_end}{breakpos};
		
		$reference =~ /(ENSG\d+)/;
		my $gene = $1;
		
		# Breakpoint positions can be defined in either genomic or transcriptome space
		my $breakpos_genomic;
		if (defined $genomic_space)
		{
			# The genomic position of the break is given
			$breakpos_genomic = $position;
			
			# Translate breakpoint strand to transcriptome space
			if ($cdna_gene{$reference}{strand} eq $strand)
			{
				$strand = "+";
			}
			else
			{
				$strand = "-";
			}
		}
		else
		{
			# Find the genomic position of the break with bias
			if ($strand eq "+")
			{
				$breakpos_genomic = calc_genomic_position($position - $splice_bias, $cdna_gene{$reference}) + $splice_bias;
			}
			elsif ($strand eq "-")
			{
				$breakpos_genomic = calc_genomic_position($position + $splice_bias, $cdna_gene{$reference}) - $splice_bias;
			}
		}
			
		# Find exons that occur before and after the break
		# Find position of break in each cdna
		my @exonsbefore;
		my @exonsafter;
		foreach my $transcript (@{$gene_transcripts{$gene}})
		{
			my @exons = @{$cdna{$transcript}{exons}};
			
			# Find position of break in cdna space
			my $breakpos_transcript = calc_transcript_position($breakpos_genomic, $cdna{$transcript});
			die if not defined $breakpos_transcript or $breakpos_transcript < 0;
				
			# Add the cdna pos for this transcript and fusion
			$transcript_fusion_position{$transcript}{$cluster_id} = $breakpos_transcript;
	
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
		
		if ($cdna_gene{$reference}{strand} eq "-")
		{
			($exons_before_size,$exons_after_size) = ($exons_after_size,$exons_before_size);
		}
		elsif ($cdna_gene{$reference}{strand} ne "+")
		{
			die;
		}
		
		$fusion_gene_exons_before_size{$cluster_id}{$gene} = $exons_before_size;
		$fusion_gene_exons_after_size{$cluster_id}{$gene} = $exons_after_size;
	}
}

my %fusion_gene_counts_before;
my %fusion_gene_counts_after;
foreach my $transcript (keys %transcript_fusion_position)
{
	$transcript =~ /(ENSG\d+)/;
	my $gene = $1;

	open TA, "$samtools_bin view $cdna_bam_filename '$transcript' |" or die "Error: Unable to run samtools on $cdna_bam_filename: $!\n";
	while (<TA>)
	{
		my @sam_fields = split /\t/;
			
		my $rname = $sam_fields[2];
		my $pos = $sam_fields[3];
		my $seq = $sam_fields[9];
			
		$transcript eq $rname or die "Error: samtools retrieived alignments to $rname when alignments to $transcript were requested\n";
			
		my $start = $pos;
		my $end = $pos + length($seq) - 1;

		foreach my $cluster_id (keys %{$transcript_fusion_position{$transcript}})
		{
			my $fusion_cdnapos = $transcript_fusion_position{$transcript}{$cluster_id};

			if ($end < $fusion_cdnapos)
			{
				$fusion_gene_counts_before{$cluster_id}{$gene} += length($seq);
			}
			elsif ($start > $fusion_cdnapos)
			{
				$fusion_gene_counts_after{$cluster_id}{$gene} += length($seq);
			}
			else
			{
				$fusion_gene_counts_before{$cluster_id}{$gene} += $fusion_cdnapos - $start;
				$fusion_gene_counts_after{$cluster_id}{$gene} += $end - $fusion_cdnapos;				
			}
		}
	}
	close TA;
}

foreach my $cluster_id (keys %breaks)
{
	foreach my $cluster_end (keys %{$breaks{$cluster_id}})
	{
		my $reference = $breaks{$cluster_id}{$cluster_end}{reference};
		my $strand = $breaks{$cluster_id}{$cluster_end}{strand};
		my $position = $breaks{$cluster_id}{$cluster_end}{breakpos};
		
		$reference =~ /(ENSG\d+)/;
		my $gene = $1;

		my $count_before = $fusion_gene_counts_before{$cluster_id}{$gene};
		my $count_after = $fusion_gene_counts_after{$cluster_id}{$gene};
		
		my $size_before = $fusion_gene_exons_before_size{$cluster_id}{$gene};
		my $size_after = $fusion_gene_exons_after_size{$cluster_id}{$gene};
		
		if ($strand eq "-")
		{
			($count_before,$count_after) = ($count_after,$count_before);
			($size_before,$size_after) = ($size_after,$size_before);
		}

		$count_before = 0 if not defined $count_before;
		$count_after = 0 if not defined $count_after;
		
		print $cluster_id."\t".$cluster_end."\t".$gene."\t".$size_before."\t".$size_after."\t".$count_before."\t".$count_after."\n";
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


