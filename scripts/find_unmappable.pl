#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin";
use gene_models;

my $gene_models_filename = shift;
my $max_alignments = shift;

defined $gene_models_filename or die "Usage: $0 gene_models_filename max_alignments < sam_alignments > unmappable_readids\n";

my $gene_models = gene_models->new($gene_models_filename);

my $bin_spacing = 200000;

my %genome_positions;
while (<>)
{
	next if /^\@/;

	chomp;
	my @fields = split /\t/;

	my $read_id = $fields[0];
	my $flag = $fields[1];
	my $rname = $fields[2];
	my $pos = $fields[3];
	my $seq = $fields[9];
	
	foreach my $info_idx (11..$#fields)
	{
		if ($fields[$info_idx] =~ /^XM:i:(\d+)$/)
		{
			if ($1 > 0)
			{
				$read_id =~ /(.*)\/([12])/;
				my $fragment_index = $1;
				print $fragment_index."\n";
				last;
			}
		}
	}
	
	next if ($flag & hex('0x0004'));

	my $start = $pos;
	my $end = $pos + length($seq) - 1;

	my $chromosome = $gene_models->calc_genomic_chromosome($rname);
	my $genome_start = $gene_models->calc_genomic_position($rname, $start);
	my $genome_end = $gene_models->calc_genomic_position($rname, $end);
	
	if ($genome_start > $genome_end)
	{
		($genome_start,$genome_end) = ($genome_end,$genome_start);
	}
	
	push @{$genome_positions{$read_id}}, [$chromosome,$genome_start,$genome_end];
}

my %read_alignment_counts;
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
	
	$read_alignment_counts{$fragment_index}{$read_end} = $alignment_count;
}

foreach my $fragment_index (keys %read_alignment_counts)
{
	my $count_1 = $read_alignment_counts{$fragment_index}{'1'};
	my $count_2 = $read_alignment_counts{$fragment_index}{'2'};

	$count_1 = 1 if not defined $count_1;
	$count_2 = 1 if not defined $count_2;

	if ($count_1 * $count_2 > $max_alignments)
	{
		print $fragment_index."\n";
	}
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


