package gene_models;

use strict;
use warnings FATAL => 'all';
use List::Util qw[min max];

# Gene models member functions

# Gene models constructor
sub new
{
	my %accepted_feature_types;
	$accepted_feature_types{CDS} = 1;
	$accepted_feature_types{exon} = 1;
	$accepted_feature_types{start_codon} = 1;
	$accepted_feature_types{stop_codon} = 1;

	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $gene_models_filename = shift;

	my $self = {};
	
	my $line_number = 1;
	open GFF, $gene_models_filename or die "Error: Unable to open $gene_models_filename\n";
	while (<GFF>)
	{
		chomp;
		next if /^#/;
		my @gff_fields = split /\t/;
	
		my $chromosome = $gff_fields[0];
		my $source = $gff_fields[1];
		my $feature_type = $gff_fields[2];
		my $start = $gff_fields[3];
		my $end = $gff_fields[4];
		my $strand = $gff_fields[6];
		my @features = split /;/, $gff_fields[8];

		next unless $accepted_feature_types{$feature_type};
		
		my $gene_id;
		my $transcript_id;
		my $gene_name;
		foreach my $feature (@features)
		{
			if ($feature =~ /(\S+)\s+(.*)/)
			{
				my $key = $1;
				my $value = $2;
	
				$value =~ s/"//g;
	
				$gene_id = $value if $key eq "gene_id";
				$transcript_id = $value if $key eq "transcript_id";
				$gene_name = $value if $key eq "gene_name";
			}
		}
		
		defined $gene_id or die "Error: gtf line $line_number has no gene_id\n";
		defined $transcript_id or die "Error: gtf line $line_number has no transcript_id\n";
		defined $gene_name or die "Error: gtf line $line_number has no gene_name\n";
		
		$transcript_id = $gene_id."|".$transcript_id;		
		
		$self->{transcripts}{$transcript_id}{gene} = $gene_id;
		$self->{transcripts}{$transcript_id}{chromosome} = $chromosome;
		$self->{transcripts}{$transcript_id}{strand} = $strand;
		$self->{transcripts}{$transcript_id}{source} = $source;
		push @{$self->{transcripts}{$transcript_id}{exons}}, [$start,$end] if $feature_type eq "exon";
		push @{$self->{transcripts}{$transcript_id}{cds}}, [$start,$end] if $feature_type eq "CDS";
		
		$self->{genes}{$gene_id}{name} = $gene_name;
		$self->{genes}{$gene_id}{chromosome} = $chromosome;
		$self->{genes}{$gene_id}{strand} = $strand;
		$self->{genes}{$gene_id}{source} = $source;
		$self->{genes}{$gene_id}{transcripts}{$transcript_id} = 1;
		
		$self->{chromosomes}{$chromosome}{genes}{$gene_id} = 1;
		
		$self->{sources}{$source}{$gene_id} = 1;
		
		$line_number++;
	}
	close GFF;
	
	# Create null gene entry
	$self->{transcripts}{""}{gene} = "";
	$self->{transcripts}{""}{chromosome} = "";
	$self->{transcripts}{""}{strand} = "";
	$self->{transcripts}{""}{source} = "";
	push @{$self->{transcripts}{""}{exons}}, [0,0];
	push @{$self->{transcripts}{""}{cds}}, [0,0];
	$self->{genes}{""}{name} = "";
	$self->{genes}{""}{chromosome} = "";
	$self->{genes}{""}{strand} = "";
	$self->{genes}{""}{source} = "";
	$self->{genes}{""}{transcripts}{""} = 1;
	
	# Sort the exons and cds for each transcript
	foreach my $transcript_id (keys %{$self->{transcripts}})
	{
		@{$self->{transcripts}{$transcript_id}{exons}} = sort { $a->[0] <=> $b->[0] } (@{$self->{transcripts}{$transcript_id}{exons}});
		
		if (defined $self->{transcripts}{$transcript_id}{cds})
		{
			@{$self->{transcripts}{$transcript_id}{cds}} = sort { $a->[0] <=> $b->[0] } (@{$self->{transcripts}{$transcript_id}{cds}});
		}
	}
	
	# Calculate gene regions
	foreach my $gene_id (keys %{$self->{genes}})
	{
		my @starts;
		my @ends;
		foreach my $transcript_id (keys %{$self->{genes}{$gene_id}{transcripts}})
		{
			my @exons = @{$self->{transcripts}{$transcript_id}{exons}};
			
			push @starts, $exons[0]->[0];
			push @ends, $exons[$#exons]->[1];
		}
		
		$self->{genes}{$gene_id}{region} = [min(@starts), max(@ends)];
	}
	
	# Create gene region binning
	# Create nearest gene binning
	# Create adjacencies
	foreach my $chromosome (keys %{$self->{chromosomes}})
	{
		my @gene_ids = keys %{$self->{chromosomes}{$chromosome}{genes}};
		
		create_binning(\@gene_ids, $self->{genes}, "region", 10000, \%{$self->{chromosomes}{$chromosome}{overlap_genes}});
		create_nearest_binning(\@gene_ids, $self->{genes}, "region", 10000, \%{$self->{chromosomes}{$chromosome}{nearest_gene}});
		create_adjacency(\@gene_ids, $self->{genes}, "region", \%{$self->{adjacent_gene}});
	}
	
	# Calculate 3 prime and 5 prime utrs
	foreach my $transcript_id (keys %{$self->{transcripts}})
	{
		next if not defined $self->{transcripts}{$transcript_id}{cds};
	
		my @cds = @{$self->{transcripts}{$transcript_id}{cds}};
		my @exons = @{$self->{transcripts}{$transcript_id}{exons}};
		my $strand = $self->{transcripts}{$transcript_id}{strand};
	
		my $coding_start = $cds[0]->[0];
		my $coding_end = $cds[$#cds]->[1];
	
		foreach my $exon (@exons)
		{
			if ($exon->[0] < $coding_start)
			{
				my $utr_exon_start = $exon->[0];
				my $utr_exon_end = min($exon->[1], $coding_start - 1);
	
				if ($strand eq "+")
				{
					push @{$self->{transcripts}{$transcript_id}{utr5p}}, [$utr_exon_start,$utr_exon_end];
				}
				else
				{
					push @{$self->{transcripts}{$transcript_id}{utr3p}}, [$utr_exon_start,$utr_exon_end];
				}
			}
			elsif ($exon->[1] > $coding_end)
			{
				my $utr_exon_start = max($exon->[0], $coding_end + 1);
				my $utr_exon_end = $exon->[1];
	
				if ($strand eq "+")
				{
					push @{$self->{transcripts}{$transcript_id}{utr3p}}, [$utr_exon_start,$utr_exon_end];
				}
				else
				{
					push @{$self->{transcripts}{$transcript_id}{utr5p}}, [$utr_exon_start,$utr_exon_end];
				}
			}
		}
	}
	
	bless($self,$class);
	return $self;
}


sub calc_nearest_gene
{
	my $self = shift;
	my $chromosome = shift;
	my $break_pos = shift;
	
	if (not defined $self->{chromosomes}{$chromosome})
	{
		return "";
	}
	
	my @gene_ids = retrieve_nearest($self->{chromosomes}{$chromosome}{nearest_gene}, [$break_pos,$break_pos]);
	
	my $nearest_gene_id;
	my $nearest_gene_dist;
	foreach my $gene_id (@gene_ids)
	{
		my $gene_region = $self->{genes}{$gene_id}{region};

		my $dist = 0;
		if ($break_pos < $gene_region->[0])
		{
			$dist = $gene_region->[0] - $break_pos;
		}
		elsif ($break_pos > $gene_region->[1])
		{
			$dist = $break_pos - $gene_region->[1];
		}
		
		if (not defined $nearest_gene_dist or $dist < $nearest_gene_dist)
		{
			$nearest_gene_dist = $dist;
			$nearest_gene_id = $gene_id;
		}
	}
	
	return $nearest_gene_id;
}

sub calc_gene
{
	my $self = shift;
	my $ref_name = shift;
	my $break_pos = shift;
	
	if (not defined $self->{chromosomes}{$ref_name} and not defined $self->{transcripts}{$ref_name})
	{
		return "";
	}
	elsif (defined $self->{transcripts}{$ref_name})
	{
		return $self->{transcripts}{$ref_name}{gene};
	}
	else
	{
		return $self->calc_nearest_gene($ref_name, $break_pos);
	}
}

sub calc_overlapping_genes
{
	my $self = shift;
	my $ref_name = shift;
	my $region = shift;
	
	if (not defined $self->{chromosomes}{$ref_name} and not defined $self->{transcripts}{$ref_name})
	{
		return ();
	}
	
	my $chromosome = $self->calc_genomic_chromosome($ref_name);
	my @genomic_regions = $self->calc_genomic_regions($ref_name, $region);
	
	my %overlapping_gene_ids;
	foreach my $genomic_region (@genomic_regions)
	{
		my @gene_ids = retrieve_binning($self->{chromosomes}{$chromosome}{overlap_genes}, $genomic_region);
		foreach my $gene_id (@gene_ids)
		{
			if (overlap($genomic_region, $self->{genes}{$gene_id}{region}))
			{
				$overlapping_gene_ids{$gene_id} = 1;
			}
		}
	}
	
	return keys %overlapping_gene_ids;
}

sub calc_gene_location
{
	my $self = shift;
	my $gene_id = shift;
	my $break_pos = shift;
	
	if ($gene_id eq "")
	{
		return "";
	}
	
	my $gene_region = $self->{genes}{$gene_id}{region};
	my $strand = $self->{genes}{$gene_id}{strand};
	
	if (($break_pos < $gene_region->[0] && $strand eq "+") or ($break_pos > $gene_region->[1] && $strand eq "-"))
	{
		return "upstream";
	}
	
	if (($break_pos > $gene_region->[1] && $strand eq "+") or ($break_pos < $gene_region->[0] && $strand eq "-"))
	{
		return "downstream";
	}
	
	my %gene_location;
	foreach my $transcript_id (keys %{$self->{genes}{$gene_id}{transcripts}})
	{
		foreach my $cds (@{$self->{transcripts}{$transcript_id}{cds}})
		{
			if ($break_pos >= $cds->[0] and $break_pos <= $cds->[1])
			{
				$gene_location{cds} = 1;
			}
		}
		
		foreach my $utr5p (@{$self->{transcripts}{$transcript_id}{utr5p}})
		{
			if ($break_pos >= $utr5p->[0] and $break_pos <= $utr5p->[1])
			{
				$gene_location{utr5p} = 1;
			}
		}
		
		foreach my $utr3p (@{$self->{transcripts}{$transcript_id}{utr3p}})
		{
			if ($break_pos >= $utr3p->[0] and $break_pos <= $utr3p->[1])
			{
				$gene_location{utr3p} = 1;
			}
		}
	}
	
	if (defined $gene_location{cds})
	{
		return "coding";
	}
	elsif (defined $gene_location{utr5p})
	{
		return "utr5p";
	}
	elsif (defined $gene_location{utr3p})
	{
		return "utr3p";
	}
	else
	{
		return "intron";
	}
}

# Check if a reference sequence is a transcript
sub is_transcript
{
	my $self = shift;
	my $reference_id = shift;
	
	return defined $self->{transcripts}{$reference_id};
}

# Find chromosome in genome given a reference name from the transcriptome or genome
sub calc_genomic_chromosome
{
	my $self = shift;
	my $transcript_id = shift;
	
	if (not $self->is_transcript($transcript_id))
	{
		return $transcript_id;
	}
	
	return $self->{transcripts}{$transcript_id}{chromosome};
}

# Find position in genome given a position in the transcriptome or genome
sub calc_genomic_position
{
	my $self = shift;
	my $transcript_id = shift;
	my $position = shift;
	
	if (not $self->is_transcript($transcript_id))
	{
		return $position;
	}

	my $transcript_ref = $self->{transcripts}{$transcript_id};
	
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

# Find genomic regions given a region and the strand and exons of the transcript
sub calc_genomic_regions
{
	my $self = shift;
	my $transcript_id = shift;
	my $region = shift;
	
	if (not $self->is_transcript($transcript_id))
	{
		return $region;
	}
	
	my $transcript_ref = $self->{transcripts}{$transcript_id};
	
	my $strand = $transcript_ref->{strand};
	my $exons = $transcript_ref->{exons};
	
	my $transcript_length = regions_length(@{$exons});
	
	if ($strand eq "-")
	{
		$region = [$transcript_length - $region->[1] + 1, $transcript_length - $region->[0] + 1];
	}
	
	if ($region->[0] < 1)
	{
		$region = [1, $region->[1]];
	}

	if ($region->[1] > $transcript_length)
	{
		$region = [$region->[0], $transcript_length];
	}
	
	my @genomic;
	my $local_offset = 0;
	foreach my $exon (@{$exons})
	{
		my $exonsize = $exon->[1] - $exon->[0] + 1;
		
		my $local_start = $region->[0] - $local_offset;
		my $local_end = $region->[1] - $local_offset;
		
		my $overlap_start = max(1,$local_start) + $exon->[0] - 1;
		my $overlap_end = min($exonsize,$local_end) + $exon->[0] - 1;
		
		if ($overlap_start <= $overlap_end)
		{
			push @genomic, [$overlap_start,$overlap_end];
		}
		
		$local_offset += $exonsize;
	}
	
	return @genomic;
}

# Find position in genome given a position in an exon
sub exon_to_genome
{
	my $self = shift;
	my $exon_id = shift;
	my $position = shift;
	
	my @exon_id_fields = split /\|/, $exon_id;
	
	scalar @exon_id_fields == 3 or die "Error: $exon_id is not an exon id\n";
	
	my $transcript_id = $exon_id_fields[0]."|".$exon_id_fields[1];
	my $exon_number = $exon_id_fields[2];
	
	$self->is_transcript($transcript_id) or die "Error: unable to recognize transcript for $exon_id\n";
	
	my $transcript_ref = $self->{transcripts}{$transcript_id};
	
	my $strand = $transcript_ref->{strand};
	my $exons = $transcript_ref->{exons};
	
	$exon_number < scalar @{$exons} or die "Error: invalid exon number for $exon_id\n";
	
	my $exon = $exons->[$exon_number];
	my $exon_size = $exon->[1] - $exon->[0] + 1;
	
	if ($strand eq "-")
	{
		$position = $exon_size - $position + 1;
	}
	
	return $exon->[0] + $position - 1;
}

# Find strand in genome given a strand in transcriptome or genome
sub calc_genomic_strand
{
	my $self = shift;
	my $transcript_id = shift;
	my $strand = shift;
	
	if (not $self->is_transcript($transcript_id))
	{
		return $strand;
	}
	elsif ($self->{transcripts}{$transcript_id}{strand} eq $strand)
	{
		return "+";
	}
	else
	{
		return "-";
	}
}

# Find position in a transcript given a genomic position and strand and exons of the transcript
# This version returns the position of the beginning of the next exon if the genomic position is intronic
sub calc_transcript_position
{
	my $self = shift;
	my $transcript_id = shift;
	my $position = shift;
	
	my $transcript_ref = $self->{transcripts}{$transcript_id};
	
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

# Find strand in transcript given a strand in the genome
sub calc_transcript_strand
{
	my $self = shift;
	my $transcript_id = shift;
	my $strand = shift;
	
	if ($self->{transcripts}{$transcript_id}{strand} eq $strand)
	{
		return "+";
	}
	else
	{
		return "-";
	}
}



# Gene models helper functions

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

# Expand to the largest region containing both regions
sub expand
{
	my $region1 = $_[0];
	my $region2 = $_[1];
	
	my $expanded_start = min($region1->[0], $region2->[0]);
	my $expanded_end = max($region1->[1], $region2->[1]);

	return [$expanded_start, $expanded_end];
}

# Calculate bins overlapped by a region
sub get_bins
{
	my $region_ref = shift;
	my $bin_spacing = shift;

	my $start = $region_ref->[0];
	my $end = $region_ref->[1];
	
	my $start_bin = int($start / $bin_spacing);
	my $end_bin = int($end / $bin_spacing);
	
	return ($start_bin .. $end_bin);
}

# Create a binning structure
sub create_binning
{
	my $ids_ref = shift;
	my $regions_ref = shift;
	my $regions_name = shift;
	my $spacing = shift;
	my $bins_ref = shift;
	
	foreach my $region_id (@{$ids_ref})
	{
		foreach my $bin (get_bins($regions_ref->{$region_id}{$regions_name}, $spacing))
		{
			push @{$bins_ref->{bins}{$bin}}, $region_id;
		}
	}
	
	$bins_ref->{spacing} = $spacing;
	$bins_ref->{maxbin} = max(keys %{$bins_ref->{bins}});
}

# Create a binning structure, including nearest regions
sub create_nearest_binning
{
	my $ids_ref = shift;
	my $regions_ref = shift;
	my $regions_name = shift;
	my $spacing = shift;
	my $bins_ref = shift;
	
	create_binning($ids_ref, $regions_ref, $regions_name, $spacing, \%{$bins_ref});
	
	my $max_bin = $bins_ref->{maxbin};
	
	my %nearest;		

	my $current_gene;
	my $current_gene_end;
	foreach my $bin (0..$max_bin)
	{
		if (defined $current_gene)
		{
			push @{$nearest{$bin}}, $current_gene;
		}
		
		foreach my $gene_id (@{$bins_ref->{bins}{$bin}})
		{
			if (not defined $current_gene_end or $bins_ref->{genes}{$gene_id}{region}->[1] > $current_gene_end)
			{
				$current_gene = $gene_id;
				$current_gene_end = $bins_ref->{genes}{$gene_id}{region}->[1];
			}
		}
	}

	undef $current_gene;
	my $current_gene_start;
	foreach my $bin (reverse(0..$max_bin))
	{
		if (defined $current_gene)
		{
			push @{$nearest{$bin}}, $current_gene;
		}
		
		foreach my $gene_id (@{$bins_ref->{bins}{$bin}})
		{
			if (not defined $current_gene_start or $bins_ref->{genes}{$gene_id}{region}->[0] < $current_gene_start)
			{
				$current_gene = $gene_id;
				$current_gene_start = $bins_ref->{genes}{$gene_id}{region}->[0];
			}
		}
	}
	
	foreach my $bin (keys %nearest)
	{
		push @{$bins_ref->{bins}{$bin}}, @{$nearest{$bin}};
	}
}

# Retrieve ids from a binning structure
sub retrieve_binning
{
	my $bins_ref = shift;
	my $region_ref = shift;
	
	my $spacing = $bins_ref->{spacing};
	
	my @region_ids;
	foreach my $bin (get_bins($region_ref, $spacing))
	{
		push @region_ids, @{$bins_ref->{bins}{$bin}} if defined $bins_ref->{bins}{$bin};
	}
	
	return @region_ids;
}

# Retrieve nearest id from a binning structure
sub retrieve_nearest
{
	my $bins_ref = shift;
	my $region_ref = shift;
	
	my $spacing = $bins_ref->{spacing};
	
	my @region_ids;
	foreach my $bin (get_bins($region_ref, $spacing))
	{
		next if not defined $bins_ref->{bins}{$bin};
		push @region_ids, @{$bins_ref->{bins}{$bin}};
	}
	
	if (scalar @region_ids == 0)
	{
		push @region_ids, @{$bins_ref->{bins}{$bins_ref->{maxbin}}};
	}
	
	return @region_ids;
}

# Create adjacency structure for list of regions
sub create_adjacency
{
	my $ids_ref = shift;
	my $regions_ref = shift;
	my $regions_name = shift;
	my $adjacency_ref = shift;
	
	my @sorted_ids = sort { $regions_ref->{$a}{$regions_name}->[0] <=> $regions_ref->{$b}{$regions_name}->[0] } (@{$ids_ref});

	foreach my $index1 (0..$#sorted_ids)
	{
		my $id1 = $sorted_ids[$index1];
		my $region1 = $regions_ref->{$id1}{$regions_name};

		my @neighbours;
		push @neighbours, $id1;

		my $neighbourhood = $region1;
		my $bridged_gaps = 0;

		foreach my $index2 ($index1 + 1..$#sorted_ids)
		{
			my $id2 = $sorted_ids[$index2];
			my $region2 = $regions_ref->{$id2}{$regions_name};

			if (not overlap($neighbourhood,$region2))
			{
				$bridged_gaps++;
			}

			last if $bridged_gaps == 2;

			$neighbourhood = expand($neighbourhood,$region2);

			push @neighbours, $id2;
		}

		foreach my $neighbour1 (@neighbours)
		{
			foreach my $neighbour2 (@neighbours)
			{
				next if $neighbour1 eq $neighbour2;
				$adjacency_ref->{$neighbour1}{$neighbour2} = 1;
			}
		}
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


1;
