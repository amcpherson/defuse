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
use gene_models;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Annotate fusions\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -o, --output    Output Directory\n";

my $help;
my $config_filename;
my $dataset_directory;
my $output_directory;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'output=s'    => \$output_directory,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $dataset_directory or die @usage;
defined $output_directory or die @usage;

my $source_directory = abs_path("$FindBin::RealBin/../");

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

# Config values
my $gene_models_filename = $config->get_value("gene_models");

my $bin_spacing = 200000;

# Read in the gene models
my $gene_models = gene_models->new($gene_models_filename);

my %solution_clusters;
my $clusters_sc_filename = $output_directory."/clusters.sc";
read_solution_clusters($clusters_sc_filename, \%solution_clusters);

my $reads_split_catalog = $output_directory."/reads.split.catalog";
my @split_fastq_prefixes = readsplitcatalog($reads_split_catalog);
count_alignments(\@split_fastq_prefixes, \%solution_clusters);

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
		my ($cluster_id,$cluster_end,$fragment_index,$read_end) = split /\t/;
		
		$solution_clusters_ref->{$fragment_index}{cluster_id} = $cluster_id;
	}
	close SC;
}

sub count_alignments
{
	my $split_fastq_prefixes_ref = shift;
	my $cluster_reads_hash_ref = shift;
	
	my %genome_positions;
	foreach my $split_fastq_prefix (@{$split_fastq_prefixes_ref})
	{
		my $spanning_filelist = $split_fastq_prefix.".spanning.filelist";
		open SFL, $spanning_filelist or die;
		while (<SFL>)
		{
			chomp;
			my ($chr1,$chr2,$filename) = split /\t/;
			
			open SAL, $filename or die "Error: Unable to open $filename: $!\n";
			while (<SAL>)
			{
				chomp;
				my @fields = split /\t/;
				
				my $fragment_index = $fields[0];
				my $read_end = $fields[1] + 1;
				my $rname = $fields[2];
				my $strand = $fields[3];
				my $start = $fields[4];
				my $end = $fields[5];
				
				next unless defined $cluster_reads_hash_ref->{$fragment_index};
				
				my $read_id = $fragment_index."/".$read_end;
				my $chromosome = $gene_models->calc_genomic_chromosome($rname);
				my $genome_start = $gene_models->calc_genomic_position($rname, $start);
				my $genome_end = $gene_models->calc_genomic_position($rname, $end);
				
				if ($genome_start > $genome_end)
				{
					($genome_start,$genome_end) = ($genome_end,$genome_start);
				}
				
				push @{$genome_positions{$read_id}}, [$chromosome,$genome_start,$genome_end];
			}
			close SAL;
		}
		close SFL;
	}
	
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

