#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];

use lib dirname($0);
use configdata;
use cmdrunner;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Create integrated analysis of multiple libraries\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -n, --names     Library Names Filename\n";

my $help;
my $config_filename;
my $output_directory;
my $library_names;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'output=s'    => \$output_directory,
	'names=s'     => \$library_names,
);

not defined $help or usage() and exit;

defined $config_filename or die @usage;
defined $output_directory or die @usage;
defined $library_names or die @usage;

# Config values
my $config = configdata->new();
$config->read($config_filename);
my $cdna_gene_fasta = $config->get_value("cdna_gene_fasta");
my $cdna_fasta = $config->get_value("cdna_fasta");
my $samtools_bin = $config->get_value("samtools_bin");
my $scripts_directory = $config->get_value("scripts_directory");

# Necessary scripts
my $merge_read_stats_script = "$scripts_directory/merge_read_stats.pl";
my $merge_expression_script = "$scripts_directory/merge_expression.pl";

# Samtools fasta indices required
my $cdna_fasta_index = $cdna_fasta.".fai";
my $cdna_gene_fasta_index = $cdna_gene_fasta.".fai";

# Ensure required files exist
verify_file_exists($cdna_fasta);
verify_file_exists($cdna_fasta_index);
verify_file_exists($cdna_gene_fasta);
verify_file_exists($cdna_gene_fasta_index);

# Create cmdrunner
my $log_directory = dirname($0)."/log";
my $log_prefix = $log_directory."/multilib_integrate";
mkdir $log_directory if not -e $log_directory;
my $runner = cmdrunner->new();
$runner->name("multilib_integrate");
$runner->prefix($log_prefix);

# Read in library info
my %libraries;
open LNS, $library_names or die "Error: Unable to open $library_names\n$!\n";
while (<LNS>)
{
	chomp;
	my ($library_name,$library_directory) = split /\t/;
	$libraries{$library_name} = $library_directory;
}
close LNS;

my $reads_end_1_fastq = $output_directory."/reads.1.fastq";
my $reads_end_2_fastq = $output_directory."/reads.2.fastq";
my $read_stats = $output_directory."/concordant.read.stats";
my $expression = $output_directory."/expression.txt";
my $cdna_pair_bam = $output_directory."/cdna.pair.bam";
my $discordant_aligned_bam = $output_directory."/discordant.aligned.bam";
my $discordant_unaligned_bam = $output_directory."/discordant.unaligned.bam";
my $fragment_translate = $output_directory."/fragment.translate";

my @all_reads_end_1_fastq;
my @all_reads_end_2_fastq;
my @all_read_stats;
my @all_expression;
my @all_cdna_pair_sam;
my @all_discordant_aligned_sam;
my @all_discordant_unaligned_sam;

# Translate read ids
my $current_fragment_offset = 0;
open FTR, ">".$fragment_translate or die "Error: Unable to write to $fragment_translate\n$!\n";
foreach my $library_name (keys %libraries)
{
	my $library_directory = $libraries{$library_name};

	my $library_reads_end_1_fastq = $library_directory."/reads.1.fastq";
	my $library_reads_end_2_fastq = $library_directory."/reads.2.fastq";

	my $library_read_stats = $library_directory."/concordant.read.stats";
	my $library_expression = $library_directory."/expression.txt";
	
	my $library_cdna_pair_bam = $library_directory."/cdna.pair.bam";
	my $library_discordant_aligned_bam = $library_directory."/discordant.aligned.bam";
	my $library_discordant_unaligned_bam = $library_directory."/discordant.unaligned.bam";
	
	my $translated_reads_end_1_fastq = $output_directory."/".$library_name.".reads.1.fastq";
	my $translated_reads_end_2_fastq = $output_directory."/".$library_name.".reads.2.fastq";
	my $translated_cdna_pair_sam = $output_directory."/".$library_name.".cdna.pair.sam";
	my $translated_discordant_aligned_sam = $output_directory."/".$library_name.".discordant.aligned.sam";
	my $translated_discordant_unaligned_sam = $output_directory."/".$library_name.".discordant.unaligned.sam";
	
	my $max_fragment_index1 = translate_fastq($library_reads_end_1_fastq, $translated_reads_end_1_fastq, $current_fragment_offset);
	my $max_fragment_index2 = translate_fastq($library_reads_end_2_fastq, $translated_reads_end_2_fastq, $current_fragment_offset);
	
	translate_readids($library_cdna_pair_bam, $translated_cdna_pair_sam, $current_fragment_offset);
	translate_readids($library_discordant_aligned_bam, $translated_discordant_aligned_sam, $current_fragment_offset);
	translate_readids($library_discordant_unaligned_bam, $translated_discordant_unaligned_sam, $current_fragment_offset);
	
	die "Error: mismatched reads for $library_name\n" if $max_fragment_index1 != $max_fragment_index2;
	my $max_fragment_index = $max_fragment_index1;
	
	print FTR $library_name."\t".$current_fragment_offset."\t".$max_fragment_index."\n";
	
	$current_fragment_offset = $max_fragment_index + 1;
	
	push @all_reads_end_1_fastq, $translated_reads_end_1_fastq;
	push @all_reads_end_2_fastq, $translated_reads_end_2_fastq;
	push @all_read_stats, $library_read_stats;
	push @all_expression, $library_expression;
	push @all_cdna_pair_sam, $translated_cdna_pair_sam;
	push @all_discordant_aligned_sam, $translated_discordant_aligned_sam;
	push @all_discordant_unaligned_sam, $translated_discordant_unaligned_sam;
}
close FTR;

# Merge reads
$runner->run("cat #<A > #>1", [@all_reads_end_1_fastq], [$reads_end_1_fastq]);
$runner->run("cat #<A > #>1", [@all_reads_end_2_fastq], [$reads_end_2_fastq]);

# Merge read statistics
$runner->run("$merge_read_stats_script #<A > #>1", [@all_read_stats], [$read_stats]);

# Merge expression statistics
$runner->run("$merge_expression_script #<A > #>1", [@all_expression], [$expression]);

# Create bam files
my $cdna_pair_bam_prefix = $cdna_pair_bam.".sort";
my $discordant_aligned_bam_prefix = $discordant_aligned_bam.".sort";
my $discordant_unaligned_bam_prefix = $discordant_unaligned_bam.".sort";
$runner->run("cat #<A | $samtools_bin view -bt $cdna_fasta_index - | $samtools_bin sort -o - $cdna_pair_bam_prefix > #>1", [@all_cdna_pair_sam], [$cdna_pair_bam]);
$runner->run("cat #<A | $samtools_bin view -bt $cdna_gene_fasta_index - | $samtools_bin sort -o - $discordant_aligned_bam_prefix > #>1", [@all_discordant_aligned_sam], [$discordant_aligned_bam]);
$runner->run("cat #<A | $samtools_bin view -bt $cdna_gene_fasta_index - | $samtools_bin sort -o - $discordant_unaligned_bam_prefix > #>1", [@all_discordant_unaligned_sam], [$discordant_unaligned_bam]);

sub translate_readids
{
	my $bam_input_filename = shift;
	my $sam_output_filename = shift;
	my $current_fragment_offset = shift;
	
	open OUT, ">".$sam_output_filename or die "Error: Unable to write to $sam_output_filename\n$!\n";
	open IN, "$samtools_bin view $bam_input_filename|" or die "Error: Unable to run samtools on $bam_input_filename\n$!\n";
	while (<IN>)
	{
		chomp;
		my @sam_fields = split /\t/;
		
		next if $sam_fields[0] =~ /^\@/g;
		
		my $fragment_index;
		my $read_end;
		if ($sam_fields[0] =~ /^(\d*)(\/[12])$/)
		{
			$fragment_index = $1;
			$read_end = $2;
		}
		elsif ($sam_fields[0] =~ /^(\d*)$/)
		{
			$fragment_index = $1;
			$read_end = "";
		}
		
		$fragment_index += $current_fragment_offset;
		
		$sam_fields[0] = $fragment_index.$read_end;
		
		my $sam_line = join "\t", @sam_fields;
		
		print OUT $sam_line."\n";
	}
	close IN;
	close OUT;
}

sub translate_fastq
{
	my $fastq_input_filename = shift;
	my $fastq_output_filename = shift;
	my $current_fragment_offset = shift;

	my $current_max_fragment_index = -1;
	
	open FQ, $fastq_input_filename or die "Error: Unable to read $fastq_input_filename\n$!\n";
	open OUT, ">".$fastq_output_filename or die "Error: Unable to write $fastq_output_filename\n$!\n";
	while (1)
	{
		my $readid = <FQ>;
		my $sequence = <FQ>;
		my $comment = <FQ>;
		my $quality = <FQ>;
		
		last if not defined $quality;
		
		chomp($readid);
		chomp($sequence);
		chomp($comment);
		chomp($quality);
		
		$readid =~ /^@(.+)(\/[12])/;
		my $fragment_index = $1;
		my $read_end = $2;
		
		$fragment_index += $current_fragment_offset;
		
		$current_max_fragment_index = max($current_max_fragment_index, $fragment_index);
		
		print OUT "\@$fragment_index$read_end\n$sequence\n$comment\n$quality\n";
	}
	close FQ;
	close OUT;
	
	return $current_max_fragment_index;
}

sub verify_file_exists
{
	my $filename = shift;
	
	if (not -e $filename)
	{
		die "Error: Required file $filename does not exist\n";
	}
}

