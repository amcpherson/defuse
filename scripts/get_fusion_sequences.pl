#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

use FindBin;
use lib "$FindBin::RealBin/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Create approximate sequences for paired end clusters\n";
push @usage, "  -h, --help        Displays this information\n";
push @usage, "  -r, --reference   Reference Fasta\n";
push @usage, "  -c, --clusters    Clusters Filename\n";

my $help;
my $reference_fasta;
my $clusters_filename;

GetOptions
(
	'help'         => \$help,
	'reference=s'  => \$reference_fasta,
	'clusters=s'   => \$clusters_filename,
);

not defined $help or die @usage;

defined $reference_fasta or die @usage;
defined $clusters_filename or die @usage;

my $ref_db = Bio::DB::Fasta->new($reference_fasta);

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

	my $seq = $ref_db->seq($ref_name, $start, $end);

	return $seq;
}

sub getcount
{
	my $cluster_ref = shift;

	return scalar keys %{$cluster_ref};
}

foreach my $cluster_id (keys %clusters)
{
	my $strand1 = $clusters{$cluster_id}{0}{strand};
	my $strand2 = $clusters{$cluster_id}{1}{strand};

	die if not defined $strand1 or not defined $strand2;

	my $seq1 = getseq($clusters{$cluster_id}{0});
	my $seq2 = getseq($clusters{$cluster_id}{1});

	my $seq;
	if ($strand1 eq "+" and $strand2 eq "-")
	{
		$seq = $seq1."N".$seq2;
	}
	elsif ($strand1 eq "-" and $strand2 eq "+")
	{
		$seq = $seq2."N".$seq1;
	}
	elsif ($strand1 eq "-" and $strand2 eq "-")
	{
		$seq = revcomp($seq1)."N".$seq2;
	}
	elsif ($strand1 eq "+" and $strand2 eq "+")
	{
		$seq = $seq1."N".revcomp($seq2);
	}
	else
	{
		die;
	}

	print ">".$cluster_id."\n".$seq."\n";
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

