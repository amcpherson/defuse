#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use List::Util qw[min max];

use lib dirname($0)."/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my $genome = shift;

die "Usage: $0 genome\n" if not defined $genome;

# Create an indexable genome
my $genome_db = Bio::DB::Fasta->new($genome);

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
	
	if ($strand eq "+")
	{
		return $genome_db->seq($chromosome, $position + 1, $position + 2);
	}
	elsif ($strand eq "-")
	{
		return revcomp($genome_db->seq($chromosome, $position - 2, $position - 1));
	}
}

sub calc_edit_dist
{
	my $seq1 = shift;
	my $seq2 = shift;
	
	die "Error: expected $seq1 and $seq2 to be same size\n" if length($seq1) != length$seq2;
	
	my @nt1 = split //, $seq1;
	my @nt2 = split //, $seq2;
	
	my $dist = 0;
	foreach my $ntindex (0..$#nt1)
	{
		$dist++ if $nt1[$ntindex] ne $nt2[$ntindex];
	}
	
	return $dist;
}

my $found_header = 0;
my %ind;
while (<>)
{
	chomp;
	my @fields = split /\t/;

	if (not $found_header)
	{
		foreach my $field_index (0..$#fields)
		{
			$ind{$fields[$field_index]} = $field_index;
		}

		$found_header = 1;

		print $_."\n";
		next;
	}

	if ($found_header and $fields[0] eq "cluster_id")
	{
		next;
	}

	my $cluster_id = $fields[$ind{"cluster_id"}];
	my $chromosome1 = $fields[$ind{"gene_chromosome1"}];
	my $chromosome2 = $fields[$ind{"gene_chromosome2"}];
	my $genomic_break_pos1 = $fields[$ind{"genomic_break_pos1"}];
	my $genomic_break_pos2 = $fields[$ind{"genomic_break_pos2"}];
	my $genomic_strand1 = $fields[$ind{"genomic_strand1"}];
	my $genomic_strand2 = $fields[$ind{"genomic_strand2"}];
	
	my $splice_seq1 = get_splice_seq($chromosome1, $genomic_break_pos1, $genomic_strand1);
	my $splice_seq2 = get_splice_seq($chromosome2, $genomic_break_pos2, $genomic_strand2);

	my $splice_seqf = $splice_seq1.revcomp($splice_seq2);
	my $splice_seqr = $splice_seq2.revcomp($splice_seq1);

	my @scores;
	push @scores, calc_edit_dist("GTAG", $splice_seqf);
	push @scores, calc_edit_dist("GTAG", $splice_seqr);
	push @scores, calc_edit_dist("ATAC", $splice_seqf);
	push @scores, calc_edit_dist("ATAC", $splice_seqr);
	
	my $score = 4 - min(@scores);
	
	print $cluster_id."\t".$splice_seqf."\t".$splice_seqr."\t".$score."\n";
}
