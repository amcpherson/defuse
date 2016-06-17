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

use lib "$FindBin::RealBin/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Generate expression plot for a specific gene.\n";
push @usage, "  -h, --help        Displays this information\n";
push @usage, "  -c, --config      Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -o, --output      Library Output directory\n";
push @usage, "  -r, --res         Results filename for fusion breakpoint position (default: results.tsv in Output Directory)\n";
push @usage, "  -f, --fusid       Fusion ID (optional)\n";
push @usage, "  -g, --gene        Ensembl Gene ID\n";
push @usage, "  -p, --pdf         Interrupted Expression PDF\n";

my $help;
my $config_filename;
my $dataset_directory;
my $output_directory;
my $results_filename;
my $fusion_id;
my $gene_id;
my $expression_pdf;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'output=s'    => \$output_directory,
	'res=s'       => \$results_filename,
	'fusid=s'     => \$fusion_id,
	'gene=s'      => \$gene_id,
	'pdf=s'       => \$expression_pdf,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
defined $dataset_directory or die @usage;
defined $output_directory or die @usage;
defined $gene_id or die @usage;
defined $expression_pdf or die @usage;

# Results filename for fusion breakpoint position
if (not defined $results_filename)
{
	$results_filename = $output_directory."/results.tsv";
}

my $source_directory = abs_path("$FindBin::RealBin/../");

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

my $cdna_regions			= $config->get_value("cdna_regions");
my $gene_tran_list			= $config->get_value("gene_tran_list");
my $splice_bias             = $config->get_value("splice_bias");
my $samtools_bin			= $config->get_value("samtools_bin");
my $scripts_directory		= $config->get_value("scripts_directory");

my $expression_rscript = $scripts_directory."/expression_plot.R";

# Find fusion breakpoint
my $fusion_breakpos;
my $fusion_breakstrand;
if (defined $fusion_id)
{
	my $first_line = 1;
	my %fi;
	open RES, $results_filename or die "Error: Unable to open $results_filename: $!\n";
	while (<RES>)
	{
		chomp;
		my @fields = split /\t/;
	
		if ($first_line)
		{
			foreach my $field_index (0..$#fields)
			{
				$fi{$fields[$field_index]} = $field_index;
			}
	
			$first_line = 0;
	
			defined $fi{"cluster_id"} or die "Error: No valid header for results file $results_filename\n";
	
			next;
		}
	
		next if $fields[$fi{"cluster_id"}] eq "cluster_id";
		
		my $cluster_id = $fields[$fi{"cluster_id"}];
		
		next unless $cluster_id == $fusion_id;
		
		my $gene1 = $fields[$fi{"gene1"}];
		my $gene2 = $fields[$fi{"gene2"}];
		my $genomic_break_pos1 = $fields[$fi{"genomic_break_pos1"}];
		my $genomic_break_pos2 = $fields[$fi{"genomic_break_pos2"}];
		my $genomic_strand1 = $fields[$fi{"genomic_strand1"}];
		my $genomic_strand2 = $fields[$fi{"genomic_strand2"}];

		if ($gene1 eq $gene_id)
		{
			$fusion_breakpos = $genomic_break_pos1;
			$fusion_breakstrand = $genomic_strand1;
		}
		elsif ($gene2 eq $gene_id)
		{
			$fusion_breakpos = $genomic_break_pos2;
			$fusion_breakstrand = $genomic_strand2;
		}
		else
		{
			die "Error: fusion $fusion_id is between $gene1 and $gene2\n";
		}
	}
	close RES;
}

if (defined $fusion_id and not defined $fusion_breakpos)
{
	die "Error: Unable to find fusion $fusion_id\n";
}

# Require concordant alignments bam
my $cdna_bam_filename = $output_directory."/cdna.pair.bam";

# Read in the gene transcripts
my %gene_transcripts;
read_gene_transcript($gene_tran_list, \%gene_transcripts);

# Read in the cdna regions
my %cdna;
read_regions($cdna_regions, \%cdna);

# Calculate exonic coverage
my %coverage;
my %covered_positions;
my @exonic_regions;
my $gene_strand;
foreach my $transcript (@{$gene_transcripts{$gene_id}})
{
	# Get the coverage of this transcript from the bam file
	my %transcript_coverage;
	get_coverage($cdna_bam_filename, $transcript, \%transcript_coverage);
	
	# Store the coverage at genomic positions
	foreach my $transcript_position (keys %transcript_coverage)
	{
		my $genomic_position = calc_genomic_position($transcript_position, $cdna{$transcript});

		$coverage{$genomic_position} += $transcript_coverage{$transcript_position};
		$covered_positions{$genomic_position} = 1;
	}
	
	# Create union of all exonic regions
	@exonic_regions = merge_regions(@exonic_regions, @{$cdna{$transcript}{exons}});
	
	# Store strand
	not defined $gene_strand or $gene_strand eq $cdna{$transcript}{strand} or die;
	$gene_strand = $cdna{$transcript}{strand};
}

# Create a fake transcript of exonic regions
my %transcript;
$transcript{strand} = $gene_strand;
$transcript{exons} = [@exonic_regions];

# Remap to merged exonic positions
my %exonic_coverage;
foreach my $genomic_position (keys %coverage)
{
	my $exonic_position = calc_transcript_position($genomic_position, \%transcript);
	$exonic_coverage{$exonic_position} = $coverage{$genomic_position};
}

# Output exonic position coverage to tmp file for reading by R
my $expression_tmp = $expression_pdf.".expr.tmp";
open TMP, ">".$expression_tmp or die "Error: Unable to open $expression_tmp: $!\n";
print TMP "position\texpression\n";
my $exonic_length = regions_length(@exonic_regions);
foreach my $exonic_position (1..$exonic_length)
{
	my $expression = $exonic_coverage{$exonic_position};
	$expression = 0 if not defined $expression;
		
	print TMP $exonic_position."\t";
	print TMP $expression."\n";
}
close TMP;

# Run the r script to create the pdf
if (defined $fusion_breakpos)
{
	# Translate fusion breakpos
	my $exonic_fusion_breakpos = calc_transcript_position($fusion_breakpos, \%transcript);
	my $exonic_fusion_breakstrand = ($fusion_breakstrand eq $gene_strand) ? "+" : "-";
	
	# Translate strand to an integer
	if ($exonic_fusion_breakstrand eq "+")
	{
		$exonic_fusion_breakstrand = "1";
	}
	else
	{
		$exonic_fusion_breakstrand = "-1";
	}

	my $rscript_result = system "Rscript $expression_rscript $expression_tmp $expression_pdf $exonic_fusion_breakpos $exonic_fusion_breakstrand";
	die "Error: Rscript failed\n" if $rscript_result != 0;
}
else
{
	my $rscript_result = system "Rscript $expression_rscript $expression_tmp $expression_pdf";
	die "Error: Rscript failed\n" if $rscript_result != 0;
}

# Remove files
unlink $expression_tmp;

# Get coverage of a specific transcript in a bam file
sub get_coverage
{
	my $bam_filename = shift;
	my $transcript = shift;
	my $coverage = shift;
	
	open PU, "$samtools_bin view -b $bam_filename '$transcript' | samtools pileup -c - |" or die "Error: Unable to get pileup for $transcript in $bam_filename: $!\n";
	while (<PU>)
	{
		my @pileup_fields = split /\t/;
		
		my $pos = $pileup_fields[1];
		my $numreads = $pileup_fields[7];
		
		$coverage->{$pos} = $numreads;
	}
	close PU;
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



