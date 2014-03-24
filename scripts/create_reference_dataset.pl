#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use List::Util qw[min max];
use File::Basename;

use lib dirname($0);
use cmdrunner;
use configdata;

use lib dirname($0)."/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

sub usage
{
	print("Usage: ".basename($0)." [options]\n");
	print("Create the reference dataset.\n");
	print("  -h, --help      Displays this information\n");
	print("  -c, --config    Configuration Filename\n");
}

my $help;
my $config_filename;

GetOptions
(
 	'help'        => \$help,
	'config=s'    => \$config_filename,
);

not defined $help or usage() and exit;
defined $config_filename or usage() and exit;

my $config = configdata->new();
$config->read($config_filename);

# Filenames for reference files
my $gene_models			        = $config->get_value("gene_models");
my $genome_fasta                = $config->get_value("genome_fasta");
my $chromosome_prefix           = $config->get_value("chromosome_prefix");
my $exons_fasta					= $config->get_value("exons_fasta");
my $cds_fasta			        = $config->get_value("cds_fasta");
my $rrna_fasta			        = $config->get_value("rrna_fasta");
my $cdna_regions		        = $config->get_value("cdna_regions");
my $cdna_fasta			        = $config->get_value("cdna_fasta");
my $est_fasta			        = $config->get_value("est_fasta");
my $est_split_catalog           = $config->get_value("est_split_catalog");
my $ig_gene_list		        = $config->get_value("ig_gene_list");
my $reference_fasta             = $config->get_value("reference_fasta");
my $repeats_filename            = $config->get_value("repeats_filename");
my $repeats_regions             = $config->get_value("repeats_regions");
my %gene_sources		        = $config->get_hash("gene_sources");
my %chromosomes			        = $config->get_hash("chromosomes");
my %ig_gene_sources             = $config->get_hash("ig_gene_sources");
my %rrna_gene_sources	        = $config->get_hash("rrna_gene_sources");
my @prefilter_fastas	        = $config->get_list("prefilter");
my $bowtie_build_bin	        = $config->get_value("bowtie_build_bin");
my $samtools_bin		        = $config->get_value("samtools_bin");
my $fatotwobit_bin				= $config->get_value("fatotwobit_bin");
my $scripts_directory			= $config->get_value("scripts_directory");

my $split_fasta_script = $scripts_directory."/split_fasta.pl";

sub verify_file_exists
{
	my $filename = shift;
	
	if (not -e $filename)
	{
		die "Error: Required file $filename does not exist\n";
	}
}

# Check for the required files
verify_file_exists($gene_models);
verify_file_exists($genome_fasta);

# Create cmdrunner for running bowtie-build and samtools faidx
my $log_directory = dirname($0)."/log";
my $log_prefix = $log_directory."/create_reference_dataset";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("defuse");
$runner->prefix($log_prefix);

# Create an indexable genome
my $dna_db = Bio::DB::Fasta->new($genome_fasta);

# Merge overlapping regions
sub merge_regions
{
	my @regions = @_;
	my @merged;

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

# Remove duplicate regions
sub remove_duplicates
{
	my @regions = @_;
	my @unique;
	my %uniquecheck;

	foreach my $region (@regions)
	{
		next if defined $uniquecheck{$region->[0]}{$region->[1]};

		$uniquecheck{$region->[0]}{$region->[1]} = 1;
		push @unique, $region;
	}

	return @unique;
}

# Find gaps between regions
sub region_gaps
{
	my @regions = @_;
	my @gaps;

	my $gap_start;
	foreach my $region (@regions)
	{
		if (defined $gap_start)
		{
			if ($gap_start <= $region->[0] - 1)
			{
				push @gaps, [$gap_start, $region->[0] - 1];
			}
		}

		$gap_start = $region->[1] + 1;
	}

	return @gaps;
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

# Check for overlap between regions
sub overlap
{
	my $region1 = shift;
	my $region2 = shift;
	
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
	my $region1 = shift;
	my $region2 = shift;
	
	my $expanded_start = min($region1->[0], $region2->[0]);
	my $expanded_end = max($region1->[1], $region2->[1]);

	return [$expanded_start, $expanded_end];
}

my $line_number = 1;

my %candidate_genes;
my %candidate_transcripts;
my %rrna_genes;
my %rrna_transcripts;

my %gene_name;
my %gene_transcripts;
my %gene_chromosome;
my %gene_strand;

my %transcript_gene;
my %transcript_chromosome;
my %transcript_strand;
my %transcript_cds;
my %transcript_exon;
my %transcript_tss;
my %transcript_tts;

my %chromosome_genes;

my %ig_genes;

print "Reading gene models\n";
open GFF, $gene_models or die "Error: Unable to open $gene_models\n";
while (<GFF>)
{
	chomp;
	my @gff_fields = split /\t/;

	my $chromosome = $gff_fields[0];
	my $source = $gff_fields[1];
	my $feature_type = $gff_fields[2];
	my $start = $gff_fields[3];
	my $end = $gff_fields[4];
	my $strand = $gff_fields[6];
	my @features = split /;/, $gff_fields[8];

	my $gene_id;
	my $transcript_id;
	my $exon_number;
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
			$exon_number = $value if $key eq "exon_number";
			$gene_name = $value if $key eq "gene_name";
		}
	}

	defined $gene_id or die "Error: line $line_number has no gene_id\n";
	defined $transcript_id or die "Error: line $line_number has no transcript_id\n";
	defined $exon_number or die "Error: line $line_number has no exon_number\n";
	defined $gene_name or die "Error: line $line_number has no gene_name\n";

	# Keep a list of ig genes
	$ig_genes{$gene_id} = 1 if $ig_gene_sources{$source};
	
	# Only keep requested genes from requested sources
	next unless defined $gene_sources{$source} or $rrna_gene_sources{$source};

	# Only keep requested genes from requested chromosomes
	next unless defined $chromosomes{$chromosome};
	
	# Keep a list of candidate fusion partner genes and transcripts
	$candidate_genes{$gene_id} = 1 if $gene_sources{$source};
	$candidate_transcripts{$transcript_id} = 1 if $gene_sources{$source};

	# Keep a list of rrna genes and transcripts
	$rrna_genes{$gene_id} = 1 if $rrna_gene_sources{$source};
	$rrna_transcripts{$transcript_id} = 1 if $rrna_gene_sources{$source};

	# Store gene info
	$gene_name{$gene_id} = $gene_name;
	$gene_transcripts{$gene_id}{$transcript_id} = 1;
	$gene_chromosome{$gene_id} = $chromosome;
	$gene_strand{$gene_id} = $strand;

	$chromosome_genes{$chromosome}{$gene_id} = 1;

	# Store transcript info
	$transcript_gene{$transcript_id} = $gene_id;
	$transcript_chromosome{$transcript_id} = $chromosome;
	$transcript_strand{$transcript_id} = $strand;

	push @{$transcript_cds{$transcript_id}}, [$start,$end] if $feature_type eq "CDS";
	push @{$transcript_exon{$transcript_id}}, [$start,$end] if $feature_type eq "exon";

	$transcript_tss{$transcript_id} = [$start,$end] if $feature_type eq "start_codon";
	$transcript_tts{$transcript_id} = [$start,$end] if $feature_type eq "stop_codon";

	$line_number++;
}
close GFF;

# Sort exons by start position
foreach my $transcript_id (keys %transcript_exon)
{
	$transcript_exon{$transcript_id} = [sort { $a->[0] <=> $b->[0] } (@{$transcript_exon{$transcript_id}})];
}

# Sort cds by start position
foreach my $transcript_id (keys %transcript_cds)
{
	$transcript_cds{$transcript_id} = [sort { $a->[0] <=> $b->[0] } (@{$transcript_cds{$transcript_id}})];
}

sub output_spliced_sequences
{
	my $filename = shift;
	my $transcripts = shift;
	my $genes = shift;
	my $chromosomes = shift;
	my $strands = shift;
	my $regions = shift;

	my $seqio = Bio::SeqIO->new('-file' => ">$filename", '-format' => 'fasta');

	foreach my $transcript (keys %{$transcripts})
	{
		next if not defined $regions->{$transcript};
		
		my $gene = $genes->{$transcript};
		my $chromosome = $chromosomes->{$transcript};
		my $strand = $strands->{$transcript};
		
		my $seq_name = $gene."|".$transcript;
		
		my $seq = "";
		foreach my $region (@{$regions->{$transcript}})
		{
			$seq = $seq.$dna_db->seq($chromosome, $region->[0], $region->[1]);
		}
	
		my $seqobj = Bio::Seq->new(-seq => $seq, -id => $seq_name);
		$seqobj = $seqobj->revcom() if defined $strand and $strand eq "-";
		$seqio->write_seq($seqobj);
	}
}

sub output_spliced_regions
{
	my $filename = shift;
	my $transcripts = shift;
	my $genes = shift;
	my $chromosomes = shift;
	my $strands = shift;
	my $regions = shift;

	open REG, ">".$filename or die "Error: Unable to write to $filename\n";

	foreach my $transcript (keys %{$transcripts})
	{
		next if not defined $regions->{$transcript};
		
		my $gene = $genes->{$transcript};
		my $chromosome = $chromosomes->{$transcript};
		my $strand = $strands->{$transcript};

		$strand = "+" if not defined $strand;

		print REG $gene."\t";
		print REG $transcript."\t";
		print REG $chromosome."\t";
		print REG $strand."\t";

		foreach my $region (@{$regions->{$transcript}})
		{
			print REG $region->[0]."\t";
			print REG $region->[1]."\t";
		}

		print REG "\n";
	}

	close REG;
}

sub output_unspliced_sequences
{
	my $filename = shift;
	my $transcripts = shift;
	my $genes = shift;
	my $chromosomes = shift;
	my $strands = shift;
	my $regions = shift;

	my $seqio = Bio::SeqIO->new('-file' => ">$filename", '-format' => 'fasta');

	foreach my $transcript (keys %{$transcripts})
	{
		next if not defined $regions->{$transcript};

		my $region_num = 0;
		foreach my $region (@{$regions->{$transcript}})
		{
			my $gene = $genes->{$transcript};
			my $chromosome = $chromosomes->{$transcript};
			my $strand = $strands->{$transcript};
			
			my $seq_name = $gene."|".$transcript."|".$region_num;
			
			my $seq = $dna_db->seq($chromosome, $region->[0], $region->[1]);
			my $seqobj = Bio::Seq->new(-seq => $seq, -id => $seq_name);
			$seqobj = $seqobj->revcom() if defined $strand and $strand eq "-";
			$seqio->write_seq($seqobj);

			$region_num++;
		}
	}
}

if (not cmdrunner::uptodate([$genome_fasta,$gene_models], [$exons_fasta]))
{
	print "Writing exon sequences\n";
	output_unspliced_sequences($exons_fasta, \%candidate_transcripts, \%transcript_gene, \%transcript_chromosome, \%transcript_strand, \%transcript_exon);
}

if (not cmdrunner::uptodate([$genome_fasta,$gene_models], [$cds_fasta]))
{
	print "Writing cds sequences\n";
	output_spliced_sequences($cds_fasta, \%candidate_transcripts, \%transcript_gene, \%transcript_chromosome, \%transcript_strand, \%transcript_cds);
}
	
if (not cmdrunner::uptodate([$genome_fasta,$gene_models], [$cdna_fasta,$cdna_regions]))
{
	print "Writing cdna sequences and regions\n";
	output_spliced_sequences($cdna_fasta, \%candidate_transcripts, \%transcript_gene, \%transcript_chromosome, \%transcript_strand, \%transcript_exon);
	output_spliced_regions($cdna_regions, \%candidate_transcripts, \%transcript_gene, \%transcript_chromosome, \%transcript_strand, \%transcript_exon);
}
	
if (not cmdrunner::uptodate([$genome_fasta,$gene_models], [$ig_gene_list]))
{
	print "Writing ig gene list\n";
	open IGL, ">".$ig_gene_list or die "Error: Unable to write to $ig_gene_list\n";
	foreach my $gene_id (keys %ig_genes)
	{
		print IGL $gene_id."\n";
	}
	close IGL;
}

if (not cmdrunner::uptodate([$genome_fasta,$gene_models], [$reference_fasta]))
{
	print "Writing dna and cdna sequences\n";
	$runner->run("cat $genome_fasta $cdna_fasta > $reference_fasta", [], []);
}

if (not cmdrunner::uptodate([$genome_fasta,$gene_models], [$rrna_fasta]))
{
	print "Writing rRNA gene sequences\n";
	output_spliced_sequences($rrna_fasta, \%rrna_transcripts, \%transcript_gene, \%transcript_chromosome, \%transcript_strand, \%transcript_exon);
}

print "Creating bowtie indices\n";
$runner->run("$bowtie_build_bin $cdna_fasta $cdna_fasta", [$cdna_fasta], [$cdna_fasta.".1.ebwt"]);
$runner->run("$bowtie_build_bin $genome_fasta $genome_fasta", [$genome_fasta], [$genome_fasta.".1.ebwt"]);
$runner->run("$bowtie_build_bin $rrna_fasta $rrna_fasta", [$rrna_fasta], [$rrna_fasta.".1.ebwt"]);
foreach my $prefilter_fasta (@prefilter_fastas)
{
	$runner->run("$bowtie_build_bin $prefilter_fasta $prefilter_fasta", [$prefilter_fasta], [$prefilter_fasta.".1.ebwt"]);
}

print "Creating samtools fasta indices\n";
$runner->run("$samtools_bin faidx $cdna_fasta", [$cdna_fasta], [$cdna_fasta.".fai"]);
$runner->run("$samtools_bin faidx $reference_fasta", [$reference_fasta], [$reference_fasta.".fai"]);

print "Splitting est file\n";
if (not cmdrunner::uptodate([$est_fasta], [$est_split_catalog]))
{
	unlink(<$est_fasta.split.*.fa>);
	unlink(<$est_fasta.split.*.fa.2bit>);
	$runner->run("$split_fasta_script #<1 1000000 $est_fasta > #>1", [$est_fasta], [$est_split_catalog]);
}
my @est_split_filenames = `cat $est_split_catalog`;
chomp(@est_split_filenames);

print "Splitting genome file\n";
my $genome_in = Bio::SeqIO->new(-format => "fasta", -file => $genome_fasta);
while (my $seq = $genome_in->next_seq())
{
	my $chromosome_fasta = $chromosome_prefix.".".$seq->id().".fa";
	if (not cmdrunner::uptodate([$genome_fasta], [$chromosome_fasta]))
	{
		my $chromosome_out = Bio::SeqIO->new(-format => "fasta", -file => ">".$chromosome_fasta);
		$chromosome_out->write_seq($seq);
		$chromosome_out->close();
	}
}
$genome_in->close();

print "Creating 2bit files\n";
$runner->run("$fatotwobit_bin #<1 #>1", [$genome_fasta], [$genome_fasta.".2bit"]);
$runner->run("$fatotwobit_bin #<1 #>1", [$cdna_fasta], [$cdna_fasta.".2bit"]);
$runner->run("$fatotwobit_bin #<1 #>1", [$rrna_fasta], [$rrna_fasta.".2bit"]);
$runner->run("$fatotwobit_bin #<1 #>1", [$cds_fasta], [$cds_fasta.".2bit"]);
$runner->run("$fatotwobit_bin #<1 #>1", [$exons_fasta], [$exons_fasta.".2bit"]);
foreach my $est_split_filename (@est_split_filenames)
{
	$runner->run("$fatotwobit_bin #<1 #>1", [$est_split_filename], [$est_split_filename.".2bit"]);
}
foreach my $chromosome (keys %chromosomes)
{
	my $chromosome_fasta = $chromosome_prefix.".".$chromosome.".fa";
	$runner->run("$fatotwobit_bin #<1 #>1", [$chromosome_fasta], [$chromosome_fasta.".2bit"]);
}

if (not cmdrunner::uptodate([$repeats_filename], [$repeats_regions]))
{
	print "Converting repeats\n";
	open REP, $repeats_filename or die "Error: Unable to open $repeats_filename: $!\n";
	open RER, ">".$repeats_regions or die "Error: Unable to open $repeats_regions: $!\n";
	while (<REP>)
	{
		chomp;
		next if /^#/;
		
		my @fields = split /\t/;
		
		my $chr = $fields[5];
		my $start = $fields[6];
		my $end = $fields[7];
		my $type = $fields[11];
		
		$chr =~ s/chr//;
		$start++;
		
		print RER $chr."\t".$start."\t".$end."\t".$type."\n";
	}
	close REP;
	close RER;
}

