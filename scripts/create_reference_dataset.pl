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
my $gene_models			= $config->get_value("gene_models");
my $genome_fasta		= $config->get_value("genome_fasta");
my $intron_regions		= $config->get_value("intron_regions");
my $intronic_fasta		= $config->get_value("intronic_fasta");
my $exons_fasta			= $config->get_value("exons_fasta");
my $cds_fasta			= $config->get_value("cds_fasta");
my $rrna_fasta			= $config->get_value("rrna_fasta");
my $utr5p_fasta			= $config->get_value("utr5p_fasta");
my $utr3p_fasta			= $config->get_value("utr3p_fasta");
my $cdna_regions		= $config->get_value("cdna_regions");
my $cdna_fasta			= $config->get_value("cdna_fasta");
my $cdna_ext_fasta		= $config->get_value("cdna_ext_fasta");
my $cdna_ext_regions	= $config->get_value("cdna_ext_regions");
my $gene_regions		= $config->get_value("gene_regions");
my $gene_fasta			= $config->get_value("gene_fasta");
my $ig_gene_list		= $config->get_value("ig_gene_list");
my $utr_extend			= $config->get_value("utr_extend");
my $cdna_utr_extend     = $config->get_value("cdna_utr_extend");
my $rrna_extend         = $config->get_value("rrna_extend");
my $cdna_gene_fasta		= $config->get_value("cdna_gene_fasta");
my $cdna_gene_regions   = $config->get_value("cdna_gene_regions");
my $tss_regions			= $config->get_value("tss_regions");
my $tts_regions			= $config->get_value("tts_regions");
my $upstream_fasta		= $config->get_value("upstream_fasta");
my $downstream_fasta	= $config->get_value("downstream_fasta");
my %gene_sources		= $config->get_hash("gene_sources");
my %chromosomes			= $config->get_hash("chromosomes");
my %ig_gene_sources		= $config->get_hash("ig_gene_sources");
my %rrna_gene_sources	= $config->get_hash("rrna_gene_sources");
my $gene_info_list		= $config->get_value("gene_info_list");
my $gene_adj_list		= $config->get_value("gene_adj_list");
my $gene_tran_list      = $config->get_value("gene_tran_list");
my @prefilter_fastas	= $config->get_list("prefilter");
my $bowtie_build_bin	= $config->get_value("bowtie_build_bin");
my $samtools_bin		= $config->get_value("samtools_bin");

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
my $genome_db = Bio::DB::Fasta->new($genome_fasta);

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
			push @gaps, [$gap_start, $region->[0] - 1];
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

my $line_number = 1;

my %candidate_genes;
my %candidate_transcripts;
my %rrna_genes;
my %rrna_transcripts;

my %gene_name;
my %gene_displayid;
my %gene_transcripts;
my %gene_chromosome;
my %gene_strand;
my %gene_extend;

my %transcript_displayid;
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

	# Only keep requested genes from requested sources
	next unless defined $gene_sources{$source} or $rrna_gene_sources{$source};

	# Only keep requested genes from requested sources
	next unless defined $chromosomes{$chromosome};
	
	# Keep a list of candidate fusion partner genes and transcripts
	$candidate_genes{$gene_id} = 1 if $gene_sources{$source};
	$candidate_transcripts{$transcript_id} = 1 if $gene_sources{$source};

	# Keep a list of rrna genes and transcripts
	$rrna_genes{$gene_id} = 1 if $rrna_gene_sources{$source};
	$rrna_transcripts{$transcript_id} = 1 if $rrna_gene_sources{$source};

	# Keep a list of ig genes
	$ig_genes{$gene_id} = 1 if $ig_gene_sources{$source};
	
	# Set gene extension dependent on source
	$gene_extend{$gene_id} = $utr_extend if $gene_sources{$source};
	$gene_extend{$gene_id} = $rrna_extend if $rrna_gene_sources{$source};

	# Store gene info
	$gene_name{$gene_id} = $gene_name;
	$gene_displayid{$gene_id} = $gene_id;
	$gene_transcripts{$gene_id}{$transcript_id} = 1;
	$gene_chromosome{$gene_id} = $chromosome;
	$gene_strand{$gene_id} = $strand;

	$chromosome_genes{$chromosome}{$gene_id} = 1;

	# Store transcript info
	$transcript_displayid{$transcript_id} = $gene_id."|".$transcript_id;
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

# Create extended cdna
my %transcript_extend_exons;
foreach my $transcript_id (keys %transcript_exon)
{
	my @extend_exons;
	foreach my $exon (@{$transcript_exon{$transcript_id}})
	{
		push @extend_exons, [$exon->[0],$exon->[1]];
	}

	$extend_exons[0]->[0] -= $cdna_utr_extend;
	$extend_exons[$#extend_exons]->[1] += $cdna_utr_extend;

	$transcript_extend_exons{$transcript_id} = [@extend_exons];
}

# Sort cds by start position
foreach my $transcript_id (keys %transcript_cds)
{
	$transcript_cds{$transcript_id} = [sort { $a->[0] <=> $b->[0] } (@{$transcript_cds{$transcript_id}})];
}

# Calculate extended gene regions
my %gene_region;
my %upstream_region;
my %downstream_region;
foreach my $gene_id (keys %gene_transcripts)
{
	my $strand = $gene_strand{$gene_id};
	
	my $start_min;
	my $end_max;
	foreach my $transcript_id (keys %{$gene_transcripts{$gene_id}})
	{
		my @exons = @{$transcript_exon{$transcript_id}};

		$start_min = $exons[0]->[0] if not defined $start_min;
		$end_max = $exons[$#exons]->[1] if not defined $end_max;

		$start_min = min($start_min, $exons[0]->[0]);
		$end_max = max($end_max, $exons[$#exons]->[1]);
	}

	my $gene_start = $start_min - $gene_extend{$gene_id};
	my $gene_end = $end_max + $gene_extend{$gene_id};

	push @{$gene_region{$gene_id}}, [$gene_start, $gene_end];
	
	if ($strand eq "+")
	{
		push @{$upstream_region{$gene_id}}, [$gene_start, $start_min - 1];
		push @{$downstream_region{$gene_id}}, [$end_max + 1, $gene_end];
	}
	elsif ($strand eq "-")
	{
		push @{$upstream_region{$gene_id}}, [$end_max + 1, $gene_end];
		push @{$downstream_region{$gene_id}}, [$gene_start, $start_min - 1];
	}
}

# Calculate adjacent candidate genes
my %gene_adjacency;
foreach my $chromosome (keys %chromosomes)
{
	my @chromosome_candidate_genes;
	foreach my $gene_id (keys %{$chromosome_genes{$chromosome}})
	{
		push @chromosome_candidate_genes, $gene_id if defined $candidate_genes{$gene_id};
	}

	my @genes_sorted = sort { $gene_region{$a}->[0]->[0] <=> $gene_region{$b}->[0]->[0] } (@chromosome_candidate_genes);

	foreach my $gene_index1 (0..$#genes_sorted)
	{
		my $gene_id1 = $genes_sorted[$gene_index1];
		my $gene_region1 = $gene_region{$gene_id1}->[0];

		my @neighbours;
		push @neighbours, $gene_id1;

		my $neighbourhood = $gene_region1;
		my $bridged_gaps = 0;

		foreach my $gene_index2 ($gene_index1 + 1..$#genes_sorted)
		{
			my $gene_id2 = $genes_sorted[$gene_index2];
			my $gene_region2 = $gene_region{$gene_id2}->[0];

			if (not overlap($neighbourhood,$gene_region2))
			{
				$bridged_gaps++;
			}

			last if $bridged_gaps == 2;

			$neighbourhood = expand($neighbourhood,$gene_region2);

			push @neighbours, $gene_id2;
		}

		foreach my $neighbour1 (@neighbours)
		{
			foreach my $neighbour2 (@neighbours)
			{
				next if $neighbour1 eq $neighbour2;
				$gene_adjacency{$neighbour1}{$neighbour2} = 1;
			}
		}
	}
}

# Calculating gene region context transcription start site list
my %transcript_gene_tss;
foreach my $transcript_id (keys %transcript_tss)
{
	my $gene_id = $transcript_gene{$transcript_id};
	my $gene_seq_start = $gene_region{$gene_id}->[0]->[0];
	my $gene_seq_end = $gene_region{$gene_id}->[0]->[1];

	my $chromosome = $transcript_chromosome{$transcript_id};
	my $strand = $transcript_strand{$transcript_id};
	
	my $tss_start = $transcript_tss{$transcript_id}->[0] - $gene_seq_start + 1;
	my $tss_end = $transcript_tss{$transcript_id}->[1] - $gene_seq_start + 1;
	if ($strand eq "-")
	{
		$tss_start = $gene_seq_end - $transcript_tss{$transcript_id}->[1] + 1;
		$tss_end = $gene_seq_end - $transcript_tss{$transcript_id}->[0] + 1;
	}
	
	push @{$transcript_gene_tss{$transcript_id}}, [$tss_start, $tss_end];
}

# Calculating gene region context transcription termination site list
my %transcript_gene_tts;
foreach my $transcript_id (keys %transcript_tts)
{
	my $gene_id = $transcript_gene{$transcript_id};
	my $gene_seq_start = $gene_region{$gene_id}->[0]->[0];
	my $gene_seq_end = $gene_region{$gene_id}->[0]->[1];
	
	my $chromosome = $transcript_chromosome{$transcript_id};
	my $strand = $transcript_strand{$transcript_id};
	
	my $tts_start = $transcript_tts{$transcript_id}->[0] - $gene_seq_start + 1;
	my $tts_end = $transcript_tts{$transcript_id}->[1] - $gene_seq_start + 1;
	if ($strand eq "-")
	{
		$tts_start = $gene_seq_end - $transcript_tts{$transcript_id}->[1] + 1;
		$tts_end = $gene_seq_end - $transcript_tts{$transcript_id}->[0] + 1;
	}
	
	push @{$transcript_gene_tts{$transcript_id}}, [$tts_start, $tts_end];
}

# Calculate 3 prime and 5 prime utrs
my %transcript_utr5p;
my %transcript_utr3p;
foreach my $transcript_id (keys %transcript_gene)
{
	next if not defined $transcript_cds{$transcript_id};

	my @cdss = @{$transcript_cds{$transcript_id}};
	my @exons = @{$transcript_exon{$transcript_id}};
	my $strand = $transcript_strand{$transcript_id};

	my $coding_start = $cdss[0]->[0];
	my $coding_end = $cdss[$#cdss]->[1];

	foreach my $exon (@exons)
	{
		if ($exon->[0] < $coding_start)
		{
			my $utr_exon_start = $exon->[0];
			my $utr_exon_end = min($exon->[1], $coding_start - 1);

			if ($strand eq "+")
			{
				push @{$transcript_utr5p{$transcript_id}}, [$utr_exon_start,$utr_exon_end];
			}
			else
			{
				push @{$transcript_utr3p{$transcript_id}}, [$utr_exon_start,$utr_exon_end];
			}
		}
		elsif ($exon->[1] > $coding_end)
		{
			my $utr_exon_start = max($exon->[0], $coding_end + 1);
			my $utr_exon_end = $exon->[1];

			if ($strand eq "+")
			{
				push @{$transcript_utr3p{$transcript_id}}, [$utr_exon_start,$utr_exon_end];
			}
			else
			{
				push @{$transcript_utr5p{$transcript_id}}, [$utr_exon_start,$utr_exon_end];
			}
		}
	}
}

# Calculate introns
my %transcript_introns;
foreach my $transcript_id (keys %transcript_gene)
{
	my @exons = @{$transcript_exon{$transcript_id}};

	# Create introns
	my @introns = region_gaps(@exons);

	$transcript_introns{$transcript_id} = [@introns];
}

# Calculate exons and intronic regions
my %gene_exons;
my %gene_intronic;
foreach my $gene_id (keys %gene_transcripts)
{
	my @exons;

	foreach my $transcript_id (keys %{$gene_transcripts{$gene_id}})
	{
		push @exons, @{$transcript_exon{$transcript_id}};
	}

	# Sort exons by start
	@exons = sort { $a->[0] <=> $b->[0] } (@exons);

	# Remove duplicate exons
	@exons = remove_duplicates(@exons);

	# Create exonic regions
	my @exonic = merge_regions(@exons);

	# Create intronic regions
	my @intronic = region_gaps(@exonic);

	$gene_exons{$gene_id} = [@exons];
	$gene_intronic{$gene_id} = [@intronic];
}

sub output_spliced_sequences
{
	my $filename = $_[0];
	my $ids = $_[1];
	my $names = $_[2];
	my $chromosomes = $_[3];
	my $strands = $_[4];
	my $regions = $_[5];

	my $seqio = Bio::SeqIO->new('-file' => ">$filename", '-format' => 'fasta');

	foreach my $id (keys %{$ids})
	{
		next if not defined $regions->{$id};

		my $seq_name = $names->{$id};
		my $chromosome = $chromosomes->{$id};
		my $strand = $strands->{$id};

		my $seq = "";
		foreach my $region (@{$regions->{$id}})
		{
			$seq = $seq.$genome_db->seq($chromosome, $region->[0], $region->[1]);
		}
	
		my $seqobj = Bio::Seq->new(-seq => $seq, -id => $seq_name);
		$seqobj = $seqobj->revcom() if defined $strand and $strand eq "-";
		$seqio->write_seq($seqobj);
	}
}

sub output_spliced_regions
{
	my $filename = $_[0];
	my $ids = $_[1];
	my $names = $_[2];
	my $chromosomes = $_[3];
	my $strands = $_[4];
	my $regions = $_[5];

	open REG, ">".$filename or die "Error: Unable to write to $filename\n";

	foreach my $id (keys %{$ids})
	{
		next if not defined $regions->{$id};

		my $seq_name = $names->{$id};
		my $chromosome = $chromosomes->{$id};
		my $strand = $strands->{$id};

		$strand = "+" if not defined $strand;

		print REG $seq_name."\t";
		print REG $chromosome."\t";
		print REG $strand."\t";

		foreach my $region (@{$regions->{$id}})
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
	my $filename = $_[0];
	my $ids = $_[1];
	my $names = $_[2];
	my $chromosomes = $_[3];
	my $strands = $_[4];
	my $regions = $_[5];

	my $seqio = Bio::SeqIO->new('-file' => ">$filename", '-format' => 'fasta');

	foreach my $id (keys %{$ids})
	{
		next if not defined $regions->{$id};

		my $region_num = 0;
		foreach my $region (@{$regions->{$id}})
		{
			my $seq_name = $names->{$id}."|".$region_num;
			my $chromosome = $chromosomes->{$id};
			my $strand = $strands->{$id};

			my $seq = $genome_db->seq($chromosome, $region->[0], $region->[1]);
			my $seqobj = Bio::Seq->new(-seq => $seq, -id => $seq_name);
			$seqobj = $seqobj->revcom() if defined $strand and $strand eq "-";
			$seqio->write_seq($seqobj);

			$region_num++;
		}
	}
}

sub output_unspliced_regions
{
	my $filename = $_[0];
	my $ids = $_[1];
	my $names = $_[2];
	my $chromosomes = $_[3];
	my $strands = $_[4];
	my $regions = $_[5];

	open REG, ">".$filename or die "Error: Unable to write to $filename\n";

	foreach my $id (keys %{$ids})
	{
		next if not defined $regions->{$id};

		my $region_num = 0;
		foreach my $region (@{$regions->{$id}})
		{
			my $name = $names->{$id}."|".$region_num;
			my $chromosome = $chromosomes->{$id};
			my $strand = $strands->{$id};

			$strand = "+" if not defined $strand;

			print REG $name."\t";
			print REG $chromosome."\t";
			print REG $strand."\t";
			print REG $region->[0]."\t";
			print REG $region->[1]."\n";

			$region_num++;
		}
	}

	close REG;
}

print "Writing exon sequences\n";
output_unspliced_sequences($exons_fasta, \%candidate_genes, \%gene_displayid, \%gene_chromosome, \%gene_strand, \%gene_exons);

# Output cds sequences
print "Writing cds sequences\n";
output_spliced_sequences($cds_fasta, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_cds);

# Output utr5p sequences
print "Writing 5 prime utr sequences\n";
output_spliced_sequences($utr5p_fasta, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_utr5p);

# Output utr3p sequences
print "Writing 3 prime utr sequences\n";
output_spliced_sequences($utr3p_fasta, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_utr3p);

# Output cdna sequences and regions
print "Writing cdna sequences and regions\n";
output_spliced_sequences($cdna_fasta, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_exon);
output_spliced_regions($cdna_regions, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_exon);

# Output extended cdna sequences and regions
print "Writing extended cdna sequences and regions\n";
output_spliced_sequences($cdna_ext_fasta, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_extend_exons);
output_spliced_regions($cdna_ext_regions, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_extend_exons);

# Output gene sequences and regions
print "Writing gene sequences and regions\n";
output_spliced_sequences($gene_fasta, \%candidate_genes, \%gene_displayid, \%gene_chromosome, \%gene_strand, \%gene_region);
output_spliced_regions($gene_regions, \%candidate_genes, \%gene_displayid, \%gene_chromosome, \%gene_strand, \%gene_region);

# Output intronic sequences
print "Writing intronic sequences\n";
output_unspliced_sequences($intronic_fasta, \%candidate_genes, \%gene_displayid, \%gene_chromosome, \%gene_strand, \%gene_intronic);

# Output intron regions
print "Writing intron regions\n";
output_unspliced_regions($intron_regions, \%candidate_transcripts, \%transcript_displayid, \%transcript_chromosome, \%transcript_strand, \%transcript_introns);

# Output upstream and downstream sequences
print "Writing upstream and downstream sequences\n";
output_spliced_sequences($upstream_fasta, \%candidate_genes, \%gene_displayid, \%gene_chromosome, \%gene_strand, \%upstream_region);
output_spliced_sequences($downstream_fasta, \%candidate_genes, \%gene_displayid, \%gene_chromosome, \%gene_strand, \%downstream_region);

# Output transcription start site regions
print "Writing transcription start site list\n";
output_unspliced_regions($tss_regions, \%candidate_transcripts, \%transcript_displayid, \%transcript_gene, {}, \%transcript_gene_tss);

# Output transcription termination site regions
print "Writing transcription termination site list\n";
output_unspliced_regions($tts_regions, \%candidate_transcripts, \%transcript_displayid, \%transcript_gene, {}, \%transcript_gene_tts);

# Output gene transcript list
print "Writing gene transcript list\n";
open GT, ">".$gene_tran_list or die "Error: Unable to write to $gene_tran_list\n";
foreach my $transcript_id (keys %candidate_transcripts)
{
	my $transcript_length = regions_length(@{$transcript_exon{$transcript_id}});
	print GT $transcript_gene{$transcript_id}."\t".$transcript_displayid{$transcript_id}."\t".$transcript_length."\n";
}
close GT;

# Output ig gene list
print "Writing ig gene list\n";
open IGL, ">".$ig_gene_list or die "Error: Unable to write to $ig_gene_list\n";
foreach my $gene_id (keys %ig_genes)
{
	print IGL $gene_id."\n";
}
close IGL;

# Create concatenated cdna genes reference
print "Writing cdna and gene sequences\n";
system "cat $cdna_ext_fasta $gene_fasta > $cdna_gene_fasta";
system "cat $cdna_ext_regions $gene_regions > $cdna_gene_regions";

# Formatting genome sequence
my $temp_genome = $genome_fasta.".tmp";
my $genome_in = Bio::SeqIO->new(-format => "fasta", -file => $genome_fasta);
my $genome_out = Bio::SeqIO->new(-format => "fasta", -file => ">".$temp_genome);
while (my $seq = $genome_in->next_seq())
{
	$seq->description("");
	$genome_out->write_seq($seq);
}
rename $temp_genome, $genome_fasta;

# Output gene info list
print "Writing gene info\n";
open GIN, ">".$gene_info_list or die "Error: Unable to write to $gene_info_list\n";
foreach my $gene_id (keys %gene_name)
{
	print GIN $gene_id."\t";
	print GIN $gene_name{$gene_id}."\t";
	print GIN $gene_chromosome{$gene_id}."\t";
	print GIN $gene_strand{$gene_id}."\t";
	print GIN $gene_region{$gene_id}->[0]->[0]."\t";
	print GIN $gene_region{$gene_id}->[0]->[1]."\n";
}
close GIN;

# Output gene adjacency list
open GAJ, ">".$gene_adj_list or die "Error: Unable to write to $gene_adj_list\n";
foreach my $gene_id1 (keys %gene_adjacency)
{
	foreach my $gene_id2 (keys %{$gene_adjacency{$gene_id1}})
	{
		print GAJ $gene_id1."\t";
		print GAJ $gene_id2."\n";
	}
}
close GAJ;

# Output rrna gene sequences
print "Writing rRNA gene sequences\n";
output_spliced_sequences($rrna_fasta, \%rrna_genes, \%gene_displayid, \%gene_chromosome, \%gene_strand, \%gene_region);

# Make bowtie indices
print "Creating bowtie indices\n";
$runner->run("$bowtie_build_bin $cdna_fasta $cdna_fasta", [], []);
$runner->run("$bowtie_build_bin $cdna_gene_fasta $cdna_gene_fasta", [], []);
$runner->run("$bowtie_build_bin $rrna_fasta $rrna_fasta", [], []);
foreach my $prefilter_fasta (@prefilter_fastas)
{
	$runner->run("$bowtie_build_bin $prefilter_fasta $prefilter_fasta", [], []);
}

# Make samtools faidx files
print "Creating samtools fasta indices\n";
$runner->run("$samtools_bin faidx $cdna_fasta", [], []);
$runner->run("$samtools_bin faidx $cdna_gene_fasta", [], []);




