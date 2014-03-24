#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use List::Util qw[min max];

$| = 1;

use lib dirname($0);
use configdata;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -o, --output    Output Directory\n";

my $help;
my $config_filename;
my $output_directory;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $output_directory or die @usage;

my $config = configdata->new();
$config->read($config_filename);

# Config values
my $samtools_bin = $config->get_value("samtools_bin");
my $cdna_gene_regions = $config->get_value("cdna_gene_regions");

my $bin_spacing = 200000;

my $clusters_sc_filename = $output_directory."/clusters.sc.sam";
my $discordant_bam = $output_directory."/discordant.aligned.bam";

my %solution_clusters;
read_solution_clusters($clusters_sc_filename, \%solution_clusters);
count_alignments($discordant_bam, \%solution_clusters);

my %alignment_counts;
foreach my $fragment_index (keys %solution_clusters)
{
	my $cluster_id = $solution_clusters{$fragment_index}{cluster_id};
	
	my $alignment_count_end_1 = $solution_clusters{$fragment_index}{alignment_counts}{"1"};
	my $alignment_count_end_2 = $solution_clusters{$fragment_index}{alignment_counts}{"2"};
	
	my $alignment_count = $alignment_count_end_1 * $alignment_count_end_2;
	
	push @{$alignment_counts{$cluster_id}}, $alignment_count;
}

foreach my $cluster_id (keys %alignment_counts)
{
	my $min_count = min(@{$alignment_counts{$cluster_id}});
	my $max_count = max(@{$alignment_counts{$cluster_id}});
	my $mean_count = calc_mean(@{$alignment_counts{$cluster_id}});
	my $num_multi = count_multi(@{$alignment_counts{$cluster_id}});
	
	print $cluster_id."\tmin_map_count\t".$min_count."\n";
	print $cluster_id."\tmax_map_count\t".$max_count."\n";
	print $cluster_id."\tmean_map_count\t".$mean_count."\n";
	print $cluster_id."\tnum_multi_map\t".$num_multi."\n";
}

sub read_solution_clusters
{
	my $solution_clusters_filename = shift;
	my $solution_clusters_ref = shift;
	
	open SC, $solution_clusters_filename or die "Error: Unable to open $solution_clusters_filename: $!\n";
	while (<SC>)
	{
		chomp;
		my ($cluster_id,$read_id) = split /\t/;		
		
		$read_id =~ /(.*)\/([12])/;
		my $fragment_index = $1;
		my $read_end = $2;
		
		$solution_clusters_ref->{$fragment_index}{cluster_id} = $cluster_id;
	}
	close SC;
}

sub count_alignments
{
	my $bam_filename = shift;
	my $cluster_reads_hash_ref = shift;
	
	# Read in the cdna regions
	my %cdna_gene;
	read_regions($cdna_gene_regions, \%cdna_gene);
	
	my %genome_positions;
	open SAM, "$samtools_bin view $bam_filename |" or die "Error: Unable to view $bam_filename with samtools: $!\n";
	while (<SAM>)
	{
		chomp;
		my @fields = split /\t/;
		
		my $read_id = $fields[0];
		my $rname = $fields[2];
		my $pos = $fields[3];
		my $seq = $fields[9];

		$read_id =~ /(.*)\/([12])/;
		my $fragment_index = $1;
		my $read_end = $2;
		
		next unless defined $cluster_reads_hash_ref->{$fragment_index};
		
		my $start = $pos;
		my $end = $pos + length($seq) - 1;
		
		my $genome_start = calc_genomic_position($start, $cdna_gene{$rname});
		my $genome_end = calc_genomic_position($end, $cdna_gene{$rname});
		
		if ($genome_start > $genome_end)
		{
			($genome_start,$genome_end) = ($genome_end,$genome_start);
		}
		
		push @{$genome_positions{$read_id}}, [$cdna_gene{$rname}{chromosome},$genome_start,$genome_end];
	}
	close SAM;
	
	foreach my $read_id (keys %genome_positions)
	{
		my %read_bins;
		foreach my $align_index (0..$#{$genome_positions{$read_id}})
		{
			my $chromosome = $genome_positions{$read_id}->[$align_index]->[0];
			my $start = $genome_positions{$read_id}->[$align_index]->[1];
			my $end = $genome_positions{$read_id}->[$align_index]->[2];
			
			foreach my $bin (get_bins($start,$end))
			{
				push @{$read_bins{$chromosome}{$bin}}, [$start,$end,$align_index];
			}
		}
		
		my %overlapping;
		foreach my $chromosome (keys %read_bins)
		{
			foreach my $bin (keys %{$read_bins{$chromosome}})
			{
				foreach my $align1 (@{$read_bins{$chromosome}{$bin}})
				{
					foreach my $align2 (@{$read_bins{$chromosome}{$bin}})
					{
						if (overlap($align1,$align2))
						{
							$overlapping{$align1->[2]}{$align2->[2]} = 1;
						}
					}
				}
			}
		}
		
		my $alignment_count = 0;
		while (scalar keys %overlapping > 0)
		{
			my $align_index = (keys %overlapping)[0];
			my @other_align_indices = keys %{$overlapping{$align_index}};
			
			foreach my $other_align_index (@other_align_indices)
			{
				delete $overlapping{$other_align_index};
			}
			
			$alignment_count++;
		}
		
		$read_id =~ /(.*)\/([12])/;
		my $fragment_index = $1;
		my $read_end = $2;
		
		$cluster_reads_hash_ref->{$fragment_index}{alignment_counts}{$read_end} = $alignment_count;
	}
}

sub calc_mean
{
	my $mean = 0;
	foreach my $value (@_)
	{
		$mean += $value;
	}
	
	$mean /= scalar @_;
	
	return $mean;
}

sub count_multi
{
	my $num_multi = 0;
	foreach my $value (@_)
	{
		$num_multi++ if $value > 1;
	}
	
	return $num_multi;
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

# Calculate bins overlapped by a region
sub get_bins
{
	my $start = shift;
	my $end = shift;
	
	my $start_bin = int($start / $bin_spacing);
	my $end_bin = int($end / $bin_spacing);
	
	return ($start_bin .. $end_bin);
}

# Check for overlap between regions
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

