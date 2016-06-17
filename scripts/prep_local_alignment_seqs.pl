#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin";
use gene_models;

use lib "$FindBin::RealBin/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Create sequences for local realignment\n";
push @usage, "  -h, --help          Displays this information\n";
push @usage, "  -r, --reference     Reference Fasta\n";
push @usage, "  -g, --genemodels    Gene Models GTF\n";
push @usage, "  -c, --clusters      Clusters Filename\n";
push @usage, "  -s, --seqrng        Sequence Range\n";

my $help;
my $reference_fasta;
my $gene_models_filename;
my $clusters_filename;
my $sequence_range;

GetOptions
(
	'help'          => \$help,
	'reference=s'   => \$reference_fasta,
	'genemodels=s'  => \$gene_models_filename,
	'clusters=s'    => \$clusters_filename,
	'seqrng=i'      => \$sequence_range,
);

not defined $help or die @usage;

defined $reference_fasta or die @usage;
defined $gene_models_filename or die @usage;
defined $clusters_filename or die @usage;
defined $sequence_range or die @usage;

my $gene_models = gene_models->new($gene_models_filename);
my $reference_db = Bio::DB::Fasta->new($reference_fasta);

my %clusters;
my %alignments;
read_clusters($clusters_filename, \%clusters, \%alignments);

sub revcomp
{
	my $sequence = shift;
	my $revcomp = reverse($sequence);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

sub getseq
{
	my $cluster_ref = shift;

	my $ref_name = $cluster_ref->{ref_name};
	my $start = $cluster_ref->{start};
	my $end = $cluster_ref->{end};

	my $seq = $reference_db->seq($ref_name, $start, $end);

	return $seq;
}

sub getcount
{
	my $cluster_ref = shift;

	return scalar keys %{$cluster_ref};
}

sub other_end
{
	my $cluster_end = shift;
	
	return "0" if $cluster_end eq "1";
	return "1" if $cluster_end eq "0";
}

sub output
{
	my $cluster_id = shift;
	my $ref_name = shift;
	my $midpoint = shift;
	my $strand = shift;
	my $other_seq = shift;
	my $other_strand = shift;
	
	my $start;
	my $end;
	if ($strand eq "+")
	{
		$start = $midpoint;
		$end = $midpoint + $sequence_range;
	}
	else
	{
		$start = $midpoint - $sequence_range;
		$end = $midpoint;
	}
	
	my $sequence = $reference_db->seq($ref_name, $start, $end);
	
	return if not defined $sequence;
	
	$sequence = revcomp($sequence) if $strand eq $other_strand;
				
	print $cluster_id."\t".$sequence."\t".$other_seq."\n";
}

foreach my $cluster_id (keys %clusters)
{
	foreach my $cluster_end ("0","1")
	{
		my $other_seq = getseq($clusters{$cluster_id}{other_end($cluster_end)});
		my $other_strand = $clusters{$cluster_id}{other_end($cluster_end)}{strand};
		
		my $ref_name = $clusters{$cluster_id}{$cluster_end}{ref_name};
		my $midpoint = ($clusters{$cluster_id}{$cluster_end}{start} + $clusters{$cluster_id}{$cluster_end}{end}) / 2;
		my $strand = $clusters{$cluster_id}{$cluster_end}{strand};
		
		my $chromosome = $gene_models->calc_genomic_chromosome($ref_name);
		my $genomic_midpoint = $gene_models->calc_genomic_position($ref_name, $midpoint);
		my $genomic_strand = $gene_models->calc_genomic_strand($ref_name, $strand);
		
		output($cluster_id, $chromosome, $genomic_midpoint, $genomic_strand, $other_seq, $other_strand);
		
		my @overlapping_genes = $gene_models->calc_overlapping_genes($clusters{$cluster_id}{$cluster_end}{ref_name}, [$clusters{$cluster_id}{$cluster_end}{start}, $clusters{$cluster_id}{$cluster_end}{end}]);
		foreach my $gene_id (@overlapping_genes)
		{
			my $gene_location = $gene_models->calc_gene_location($gene_id, $genomic_midpoint);
			
			next unless $gene_location eq "coding" or $gene_location eq "utr5p" or $gene_location eq "utr3p";
			
			foreach my $transcript_id (keys %{$gene_models->{genes}{$gene_id}{transcripts}})
			{
				my $transcript_midpoint = $gene_models->calc_transcript_position($transcript_id, $genomic_midpoint);
				my $transcript_strand = $gene_models->calc_transcript_strand($transcript_id, $genomic_strand);
				
				output($cluster_id, $transcript_id, $transcript_midpoint, $transcript_strand, $other_seq, $other_strand);
			}
		}
	}
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
		
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{ref_name} = $ref_name;
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{strand} = $strand;
		
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{start} = $start if not defined $clusters_hash_ref->{$cluster_id}{$cluster_end}{start};
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{end} = $end if not defined $clusters_hash_ref->{$cluster_id}{$cluster_end}{end};

		$clusters_hash_ref->{$cluster_id}{$cluster_end}{start} = min($start, $clusters_hash_ref->{$cluster_id}{$cluster_end}{start});
		$clusters_hash_ref->{$cluster_id}{$cluster_end}{end} = max($end, $clusters_hash_ref->{$cluster_id}{$cluster_end}{end});
	}
	close CLU;
}

