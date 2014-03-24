#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;

my @usage;
push @usage, "Usage: $0 [options]\n";
push @usage, "Divide alignments by chromosome pair.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -t, --trans     Transcript chromosome map\n";
push @usage, "  -p, --prefix    Prefix for divided files\n";

my $help;
my $trans_chr_filename;
my $output_prefix;

GetOptions
(
	'help'        => \$help,
	'trans=s'     => \$trans_chr_filename,
	'prefix=s'    => \$output_prefix,
);

not defined $help or die @usage;

defined $trans_chr_filename or die @usage;
defined $output_prefix or die @usage;

my %trans_chr;

open TRANS, $trans_chr_filename or die "Error: Unable to open $trans_chr_filename: $!\n";
while (<TRANS>)
{
	chomp;
	my @fields = split /\t/;
	$trans_chr{$fields[0]."|".$fields[1]} = $fields[2];
}
close TRANS;

my %alignment_buffer;
my %output_filename;

my $max_cached = 10000;

sub flush
{
	my $chr1 = shift;
	my $chr2 = shift;
	
	my $filename = $output_prefix.$chr1."-".$chr2;
	
	if (not defined $output_filename{$chr1}{$chr2})
	{
		unlink $filename;
	}
	
	$output_filename{$chr1}{$chr2} = $filename;
	
	open OUT, ">>".$filename or die "Error: Unable to write to $filename\n $!\n";
	foreach my $alignment (@{$alignment_buffer{$chr1}{$chr2}})
	{
		my $line = join "\t",@{$alignment};
		print OUT $line."\n";
	}
	close OUT;
	
	@{$alignment_buffer{$chr1}{$chr2}} = ();
}

sub output
{
	my $chr1 = shift;
	my $chr2 = shift;
	my $alignment = shift;
	
	push @{$alignment_buffer{$chr1}{$chr2}}, $alignment;
	
	if (scalar @{$alignment_buffer{$chr1}{$chr2}} > $max_cached)
	{
		flush($chr1,$chr2);
	}
}

sub process_alignments
{
	my $alignments_ref = shift;
	
	if (scalar keys %{$alignments_ref} < 2)
	{
		return;
	}
	
	foreach my $chr1 (keys %{$alignments_ref->{"1"}})
	{
		foreach my $chr2 (keys %{$alignments_ref->{"2"}})
		{
			my ($chr1_sorted,$chr2_sorted) = sort ($chr1,$chr2);
			
			foreach my $alignment (@{$alignments_ref->{"1"}{$chr1}})
			{
				output($chr1_sorted, $chr2_sorted, $alignment);
			}
			foreach my $alignment (@{$alignments_ref->{"2"}{$chr2}})
			{
				output($chr1_sorted, $chr2_sorted, $alignment);
			}
		}
	}
}

my $current_fragment_id;
my %current_alignments;
while (<>)
{
	my $line = $_;

	chomp;
	next if /^\@/;

	my @sam_info = split /\t/;
	
	my $qname = $sam_info[0];
	my $flag = $sam_info[1];
	my $rname = $sam_info[2];
	my $pos = $sam_info[3];
	my $mapq = $sam_info[4];
	my $cigar = $sam_info[5];
	my $mrnm = $sam_info[6];
	my $mpos = $sam_info[7];
	my $isize = $sam_info[8];
	my $seq = $sam_info[9];
	my $qual = $sam_info[10];
	my $opt = $sam_info[11];
	
	my $chr = $rname;
	if (defined $trans_chr{$rname})
	{
		$chr = $trans_chr{$rname};
	}
	
	$qname =~ /(.*)\/([12])/;
	my $fragment_id = $1;
	my $read_end = $2;
	
	my $strand = "+";
	if ($flag & hex('0x0010'))
	{
		$strand = "-";
	}
	
	if (defined $current_fragment_id and $current_fragment_id ne $fragment_id)
	{
		process_alignments(\%current_alignments);
		%current_alignments = ();
	}
	
	$current_fragment_id = $fragment_id;
	push @{$current_alignments{$read_end}{$chr}}, [$fragment_id,$read_end-1,$rname,$strand,$pos,$pos+length($seq)-1];
}

if (defined $current_fragment_id)
{
	process_alignments(\%current_alignments);
}

foreach my $chr1 (sort keys %alignment_buffer)
{
	foreach my $chr2 (sort keys %{$alignment_buffer{$chr1}})
	{
		flush($chr1,$chr2);
		print $chr1."\t".$chr2."\t".$output_filename{$chr1}{$chr2}."\n";
	}
}


