#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use List::Util qw[min max];
use File::Basename;
use Cwd qw[abs_path];

use FindBin;
use lib "$FindBin::RealBin";
use cmdrunner;
use configdata;

use lib "$FindBin::RealBin/../external/BioPerl-1.6.1";
use Bio::DB::Fasta;
use Bio::SeqIO;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Create the reference dataset.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";

my $help;
my $config_filename;
my $dataset_directory;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
);

not defined $help or die @usage;

defined $dataset_directory or die @usage;

my $source_directory = abs_path("$FindBin::RealBin/../");

if (not defined $config_filename)
{
	$config_filename = $source_directory."/scripts/config.txt";
}

-e $config_filename or die "Error: Unable to find config file $config_filename\n";

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

# Filenames for reference files
my $ensembl_version             = $config->get_value("ensembl_version");
my $ensembl_genome_version      = $config->get_value("ensembl_genome_version");
my $ucsc_genome_version         = $config->get_value("ucsc_genome_version");
my $gene_models                 = $config->get_value("gene_models");
my $genome_fasta                = $config->get_value("genome_fasta");
my $est_alignments              = $config->get_value("est_alignments");
my $unigene_fasta               = $config->get_value("unigene_fasta");
my %chromosomes                 = $config->get_hash("chromosomes");
my $chromosome_prefix           = $config->get_value("chromosome_prefix");
my $mt_chromosome               = $config->get_value("mt_chromosome");
my $exons_fasta                 = $config->get_value("exons_fasta");
my $cds_fasta                   = $config->get_value("cds_fasta");
my $rrna_fasta                  = $config->get_value("rrna_fasta");
my $cdna_regions                = $config->get_value("cdna_regions");
my $cdna_fasta                  = $config->get_value("cdna_fasta");
my $est_fasta                   = $config->get_value("est_fasta");
my @est_split_fastas            = $config->get_list("est_split_fasta");
my $ig_gene_list                = $config->get_value("ig_gene_list");
my $reference_fasta             = $config->get_value("reference_fasta");
my $repeats_filename            = $config->get_value("repeats_filename");
my $repeats_regions             = $config->get_value("repeats_regions");
my %gene_biotypes               = $config->get_hash("gene_biotypes");
my %ig_gene_biotypes            = $config->get_hash("ig_gene_biotypes");
my %rrna_gene_biotypes          = $config->get_hash("rrna_gene_biotypes");
my @prefilter_fastas            = $config->get_list("prefilter");
my $bowtie_build_bin            = $config->get_value("bowtie_build_bin");
my $samtools_bin                = $config->get_value("samtools_bin");
my $fatotwobit_bin              = $config->get_value("fatotwobit_bin");
my $gmap_build_bin              = $config->get_value("gmap_build_bin");
my $gmap_index_directory        = $config->get_value("gmap_index_directory");
my $scripts_directory           = $config->get_value("scripts_directory");

my $divide_fasta_script = $scripts_directory."/divide_fasta.pl";
my $remove_fasta_description_script = $scripts_directory."/remove_fasta_description.pl";


# Create cmdrunner for running bowtie-build and samtools faidx
my $log_directory = dirname($0)."/log";
my $log_prefix = $log_directory."/create_reference_dataset";

mkdir $log_directory if not -d $log_directory;

my $runner = cmdrunner->new();
$runner->name("defuse");
$runner->prefix($log_prefix);

mkdir $dataset_directory if not -d $dataset_directory;

sub system_with_check
{
	my $command = shift;
	
	my $result = system $command;
	die "Error: Unable to execute $command" unless $result == 0;
}

sub wget_gunzip
{
	my $url = shift;
	my $filename = shift;
	
	my $filename_gz = $filename.".gz";
	
	system_with_check("wget $url -O $filename_gz");
	system_with_check("gunzip $filename_gz");
}

sub cat_files
{
	my $input_filenames_ref = shift;
	my $output_filename = shift;
	
	my $input_filenames_list = join " ", @{$input_filenames_ref};
	my $output_filename_tmp = $output_filename.".tmp";
	
	system_with_check("cat $input_filenames_list > $output_filename_tmp");
	rename $output_filename_tmp, $output_filename;
}

# Retrieve the ensembl chromosome fastas
my %chromosome_fastas;
foreach my $chromosome (keys %chromosomes)
{
	my $chromosome_fasta = $chromosome_prefix.".".$chromosome.".fa";
	
	$chromosome_fastas{$chromosome} = $chromosome_fasta;
	
	next if -e $chromosome_fasta;
	
	my $chromosome_tmp = $chromosome_prefix.".".$chromosome.".tmp.fa";

	if ($ensembl_version < 76)
	{
		wget_gunzip("ftp://ftp.ensembl.org/pub/release-$ensembl_version/fasta/homo_sapiens/dna/Homo_sapiens.$ensembl_genome_version.$ensembl_version.dna.chromosome.$chromosome.fa.gz", $chromosome_tmp);
	}
	else
	{
		wget_gunzip("ftp://ftp.ensembl.org/pub/release-$ensembl_version/fasta/homo_sapiens/dna/Homo_sapiens.$ensembl_genome_version.dna.chromosome.$chromosome.fa.gz", $chromosome_tmp);
	}

	my $chromosome_tmp2 = $chromosome_prefix.".".$chromosome.".tmp2.fa";
	
	system_with_check("$remove_fasta_description_script < $chromosome_tmp > $chromosome_tmp2");
	
	rename $chromosome_tmp2, $chromosome_fasta;
	unlink 	$chromosome_tmp;
}

# Create genome fasta
if (not -e $genome_fasta)
{
	cat_files([values(%chromosome_fastas)], $genome_fasta);
}

# Retrieve the gene models
if (not -e $gene_models)
{
	wget_gunzip("ftp://ftp.ensembl.org/pub/release-$ensembl_version/gtf/homo_sapiens/Homo_sapiens.$ensembl_genome_version.$ensembl_version.gtf.gz", $gene_models);
}

# Retrieve the repeats
if (not -e $repeats_filename)
{
	if ($ucsc_genome_version eq "hg18")
	{
		my @chromosome_repeats_all;
		foreach my $chromosome (keys %chromosomes)
		{
			$chromosome = "M" if $chromosome eq $mt_chromosome;
			
			my $chromosome_repeats = $repeats_filename.".".$chromosome.".tmp";
			
			wget_gunzip("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chr$chromosome\_rmsk.txt.gz", $chromosome_repeats);
			
			push @chromosome_repeats_all, $chromosome_repeats;
		}
		
		cat_files(\@chromosome_repeats_all, $repeats_filename);
		unlink @chromosome_repeats_all;
	}
	else
	{
		wget_gunzip("ftp://hgdownload.cse.ucsc.edu/goldenPath/$ucsc_genome_version/database/rmsk.txt.gz", $repeats_filename);
	}
}

# Retrieve the est fasta
if (not -e $est_fasta)
{
	wget_gunzip("ftp://hgdownload.cse.ucsc.edu/goldenPath/$ucsc_genome_version/bigZips/est.fa.gz", $est_fasta);
}

# Retrieve the est alignments
if (not -e $est_alignments)
{
	if ($ucsc_genome_version eq "hg18")
	{
		my @chromosome_est_alignments_all;
		foreach my $chromosome (keys %chromosomes)
		{
			$chromosome = "M" if $chromosome eq $mt_chromosome;
			
			my $chromosome_est_alignments = $est_alignments.".".$chromosome.".tmp";
			
			wget_gunzip("ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chr$chromosome\_intronEst.txt.gz", $chromosome_est_alignments);
			
			push @chromosome_est_alignments_all, $chromosome_est_alignments;
		}
		
		cat_files(\@chromosome_est_alignments_all, $est_alignments);
		unlink @chromosome_est_alignments_all;
	}
	else
	{
		wget_gunzip("ftp://hgdownload.cse.ucsc.edu/goldenPath/$ucsc_genome_version/database/intronEst.txt.gz", $est_alignments);
	}
}

# Retrieve the unigene fasta
if (not -e $unigene_fasta)
{
	wget_gunzip("ftp://ftp.ncbi.nlm.nih.gov/repository/UniGene/Homo_sapiens/Hs.seq.uniq.gz", $unigene_fasta);
}


# Create an indexable genome
my $dna_db = Bio::DB::Fasta->new($genome_fasta);


my $line_number = 0;

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

my %accepted_feature_types;
$accepted_feature_types{CDS} = 1;
$accepted_feature_types{exon} = 1;
$accepted_feature_types{start_codon} = 1;
$accepted_feature_types{stop_codon} = 1;

print "Reading gene models\n";
open GFF, $gene_models or die "Error: Unable to open $gene_models\n";
while (<GFF>)
{
	$line_number++;

	chomp;
	next if /^#/;
	my @gff_fields = split /\t/;

	my $chromosome = $gff_fields[0];
	my $feature_type = $gff_fields[2];
	my $start = $gff_fields[3];
	my $end = $gff_fields[4];
	my $strand = $gff_fields[6];
	my @features = split /;/, $gff_fields[8];

	next unless $accepted_feature_types{$feature_type};

	my $gene_id;
	my $transcript_id;
	my $exon_number;
	my $gene_name;
	my $biotype;
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
			$biotype = $value if $key eq "gene_biotype";
		}
	}

	defined $gene_id or die "Error: line $line_number has no gene_id\n";
	defined $transcript_id or die "Error: line $line_number has no transcript_id\n";
	defined $exon_number or die "Error: line $line_number has no exon_number\n";
	defined $gene_name or die "Error: line $line_number has no gene_name\n";
	defined $biotype or die "Error: line $line_number has no biotype\n";

	# Keep a list of ig genes
	$ig_genes{$gene_id} = 1 if $ig_gene_biotypes{$biotype};
	
	# Only keep requested genes from requested biotypes
	next unless defined $gene_biotypes{$biotype} or $rrna_gene_biotypes{$biotype};

	# Only keep requested genes from requested chromosomes
	next unless defined $chromosomes{$chromosome};
	
	# Keep a list of candidate fusion partner genes and transcripts
	$candidate_genes{$gene_id} = 1 if $gene_biotypes{$biotype};
	$candidate_transcripts{$transcript_id} = 1 if $gene_biotypes{$biotype};

	# Keep a list of rrna genes and transcripts
	$rrna_genes{$gene_id} = 1 if $rrna_gene_biotypes{$biotype};
	$rrna_transcripts{$transcript_id} = 1 if $rrna_gene_biotypes{$biotype};

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

sub revcomp
{
        my $sequence = shift;
        my $revcomp = reverse($sequence);
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

sub output_spliced_sequences_polya
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
		
		$seq = revcomp($seq) if defined $strand and $strand eq "-";
		$seq .= "A" x 50;

		my $seqobj = Bio::Seq->new(-seq => $seq, -id => $seq_name);
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
	output_spliced_sequences_polya($cdna_fasta, \%candidate_transcripts, \%transcript_gene, \%transcript_chromosome, \%transcript_strand, \%transcript_exon);
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
$runner->run("$divide_fasta_script #<1 #>A", [$est_fasta], [@est_split_fastas]);

print "Creating 2bit files\n";
$runner->run("$fatotwobit_bin #<1 #>1", [$cds_fasta], [$cds_fasta.".2bit"]);
$runner->run("$fatotwobit_bin #<1 #>1", [$exons_fasta], [$exons_fasta.".2bit"]);

mkdir $gmap_index_directory if not -d $gmap_index_directory;

sub create_gmap_indices
{
	my $fasta = shift;
	my $name = shift;
	
	print "Creating $name gmap index\n";
	
	my $gmap_build_log = $gmap_index_directory."/defuse.gmap_build.$name.log";
	$runner->run("$gmap_build_bin -D $gmap_index_directory -d $name #<1 > #>1", [$fasta], [$gmap_build_log]);
}

create_gmap_indices($cdna_fasta, "cdna");

foreach my $est_split_index (0..$#est_split_fastas)
{
	create_gmap_indices($est_split_fastas[$est_split_index], "est".$est_split_index);
}

foreach my $chromosome (keys %chromosome_fastas)
{
	create_gmap_indices($chromosome_fastas{$chromosome}, "chr".$chromosome);
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

