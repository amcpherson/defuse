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
use cmdrunner;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Generate interrupted expression statistics.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -r, --results   Results filename\n";
push @usage, "  -l, --libs      Library list\n";
push @usage, "  -o, --output	Output directory\n";
push @usage, "  -s, --submit    Submitter Type (default: direct)\n";
push @usage, "  -p, --parallel  Maximum Number of Parallel Jobs (default: 1)\n";

my $help;
my $config_filename;
my $results_filename;
my $libraries_filename;
my $output_directory;
my $submitter_type;
my $max_parallel;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'results=s'   => \$results_filename,
	'libs=s'      => \$libraries_filename,
	'output=s'    => \$output_directory,
	'submit=s'    => \$submitter_type,
	'parallel=s'  => \$max_parallel,
);

not defined $help or die @usage;

defined $config_filename or die @usage;
defined $results_filename or die @usage;
defined $libraries_filename or die @usage;
defined $output_directory or die @usage;

# Submitter type defaults to direct
if (not defined $submitter_type)
{
	$submitter_type = "direct";
}

# Max parallel defaults to 1 for direct submitter
if (not defined $max_parallel)
{
	if ($submitter_type eq "direct")
	{
		$max_parallel = 1;
	}
	else
	{
		$max_parallel = 200;
	}
}

-e $config_filename or die "Error: Unable to find config file $config_filename\n";
-d $output_directory or die "Error: Unable to find output directory $output_directory\n";

$config_filename = abs_path($config_filename);
$output_directory = abs_path($output_directory);

my $config = configdata->new();
$config->read($config_filename);

my $scripts_directory = $config->get_value("scripts_directory");

my $interrupted_script = $scripts_directory."/calc_interrupted.pl";

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/interrupted";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("interrupted");
$runner->prefix($log_prefix);
$runner->maxparallel($max_parallel);
$runner->submitter($submitter_type);
$runner->jobmem(3000000000);

# Read in libraries of interest
open LIB, $libraries_filename or die "Error: Unable to open $libraries_filename: $!\n";
my @libraries = <LIB>;
chomp(@libraries);
close LIB;

# Create breaks file
my $breaks_filename = $output_directory."/interrupted.breaks.txt";
my $breaks_temp_filename = $output_directory."/interrupted.breaks.txt.tmp";

open SEQ, ">".$breaks_temp_filename or die "Error: Unable to open file $breaks_temp_filename\n";

my %library_cluster_genes;

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
	
	my $library_name = $fields[$fi{"library_name"}];
	my $cluster_id = $fields[$fi{"cluster_id"}];
	my $gene1 = $fields[$fi{"gene1"}];
	my $gene2 = $fields[$fi{"gene2"}];
	my $genomic_break_pos1 = $fields[$fi{"genomic_break_pos1"}];
	my $genomic_break_pos2 = $fields[$fi{"genomic_break_pos2"}];
	my $genomic_strand1 = $fields[$fi{"genomic_strand1"}];
	my $genomic_strand2 = $fields[$fi{"genomic_strand2"}];
	
	my $library_cluster_id = $library_name."-".$cluster_id;
	
	print SEQ $library_cluster_id."\t".$gene1."\t".$genomic_strand1."\t".$genomic_break_pos1."\n";
	print SEQ $library_cluster_id."\t".$gene2."\t".$genomic_strand2."\t".$genomic_break_pos2."\n";
	
	push @{$library_cluster_genes{$library_cluster_id}}, $gene1;
	push @{$library_cluster_genes{$library_cluster_id}}, $gene2;
}
close RES;
close SEQ;

# Replace old file if we generated a different one
cmdrunner::replaceifdifferent($breaks_temp_filename, $breaks_filename);

# Iterate through all libraries and calculate interrupted for each
foreach my $library_name (@libraries)
{
	my $library_output_directory = $output_directory."/".$library_name;
	my $library_interrupted_filename = $output_directory."/".$library_name.".interrupted.breaks";
	
	$runner->padd("$interrupted_script -g -c $config_filename -o $library_output_directory -b #<1 > #>1", [$breaks_filename], [$library_interrupted_filename]);
}
$runner->prun();

# Iterate through all libraries and collect results
my %interrupted;
foreach my $library_name (@libraries)
{
	my $library_interrupted_filename = $output_directory."/".$library_name.".interrupted.breaks";
	
	read_interrupted($library_interrupted_filename, \%{$interrupted{$library_name}});
}

print STDERR "Outputting\n";

print "\t\t\t";
foreach my $library_name (@libraries)
{
	print $library_name."\t";
}
print "\n";

foreach my $library_cluster_id (keys %library_cluster_genes)
{
	$library_cluster_id =~ /^(.+)-(\d+)$/;
	
	my $fusion_library_name = $1;
	my $cluster_id = $2;
	
	foreach my $gene (@{$library_cluster_genes{$library_cluster_id}})
	{
		print $fusion_library_name."\t".$cluster_id."\t".$gene."\t";

		foreach my $library_name (@libraries)
		{
			my $size_before = $interrupted{$library_name}{$library_cluster_id}{$gene}{size_before};
			my $size_after = $interrupted{$library_name}{$library_cluster_id}{$gene}{size_after};
			my $count_before = $interrupted{$library_name}{$library_cluster_id}{$gene}{count_before};
			my $count_after = $interrupted{$library_name}{$library_cluster_id}{$gene}{count_after};

			print $count_before."/".$size_before."/".$count_after."/".$size_after."\t";
		}
		
		print "\n";
	}
}

sub read_interrupted
{
	my $interrupted_filename = shift;
	my $interrupted_ref = shift;
	
	open BIN, $interrupted_filename or die "Error: Unable to open $interrupted_filename: $!\n";
	while (<BIN>)
	{
		chomp;
		my ($cluster_id, $gene, $size_before, $size_after, $count_before, $count_after) = split /\t/;
		$interrupted_ref->{$cluster_id}{$gene}{size_before} = $size_before;
		$interrupted_ref->{$cluster_id}{$gene}{size_after} = $size_after;
		$interrupted_ref->{$cluster_id}{$gene}{count_before} = $count_before;
		$interrupted_ref->{$cluster_id}{$gene}{count_after} = $count_after;
	}
	close BIN;
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
sub get_genomic_position
{
	my $position = shift;
	my $strand = shift;
	my $exons = shift;
	
	if ($strand eq "-")
	{
		$position = regions_length(@{$exons}) - $position + 1;
	}
	elsif ($strand ne "+")
	{
		die;
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
	
	return $position - $local_offset + ${$exons}[$#{$exons}]->[1];
}

# Find position in a transcript given a genomic position and strand and exons of the transcript
# This version returns the position of the beginning of the next exon if the genomic position is intronic
sub get_transcript_position
{
	my $position = shift;
	my $strand = shift;
	my $exons = shift;
	
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
	elsif ($strand ne "+")
	{
		die;
	}
	
	return $transcript_position;
}


