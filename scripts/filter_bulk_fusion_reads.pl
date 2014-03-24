#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];
use File::Temp qw[tempfile tempdir];
use DBI;

sub usage_exit
{
	print STDERR "Usage: $0 [options]\n";
	print STDERR "Filter for potential fusion reads.\n";
	print STDERR "  -h, --help      Displays this information\n";
	print STDERR "  -w, --working   Working Directory\n";
	print STDERR "  -n, --name      Library Name\n";
	exit;
}

my $help;
my $working_directory;
my $library_name;
my $reference_pairs_filename;

GetOptions
(
	'help'        => \$help,
	'working=s'   => \$working_directory,
	'name=s'      => \$library_name,
	'refs=s'      => \$reference_pairs_filename,
);

not defined $help or usage_exit();

defined $working_directory or usage_exit();
defined $library_name or usage_exit();

my $library_working_directory = $working_directory."/".$library_name;
my $fusion_fastq_directory = $library_working_directory."/fusion_seqs";

-e $working_directory or die "Error: Working directory $working_directory does not exist.\n";
-e $library_working_directory or die "Error: Library working directory $library_working_directory does not exist.\n";

my $paired_end_1_fastq = $library_working_directory."/paired.1.fastq";
my $paired_end_2_fastq = $library_working_directory."/paired.2.fastq";

my $cdna_end_1_bwtout = $library_working_directory."/cdna.end.1.bwtout";
my $cdna_end_2_bwtout = $library_working_directory."/cdna.end.2.bwtout";

# Connect to fusions database
my $fusions_db_name = "ovarian_fusions";
my $fusions_db_host = "shannon";
my $fusions_db_user = "ovmut_rw";
my $fusions_db_password = "o9v8m7";

my $fusions_dsn = "DBI:mysql:database=$fusions_db_name;host=$fusions_db_host";
my $fusions_dbh = DBI->connect($fusions_dsn,$fusions_db_user);#,$fusions_db_password);
if (not $fusions_dbh)
{
	print STDERR "Failed to connect to $fusions_db_name database\n";
}

# Select library id
my $sth_select_library_id = $fusions_dbh->prepare
("
	select library_id
	from library
	where library_name = ?
");

# Select genes
my $sth_select_genes = $fusions_dbh->prepare
("
	select distinct gene
	from read_pair natural join read_sequence natural join alignment
	where library_id = ?
");

$sth_select_library_id->execute($library_name);
my $library_id = $sth_select_library_id->fetchrow_array();

# Find all genes involved in potential fusions
print "Finding all genes involved in potential fusions\n";
my %fusion_gene;
$sth_select_genes->execute($library_id);
while (my $gene = $sth_select_genes->fetchrow_array())
{
	$fusion_gene{$gene} = 1;
}

# Removing previous fusion fastq directory
if (-d $fusion_fastq_directory)
{
	print "Removing previous $fusion_fastq_directory directory\n";
	system "rm -rf $fusion_fastq_directory";
}
mkdir $fusion_fastq_directory;

# Get read ids for genes of interest
print "Retrieving read ids of reads with alignments to fusion genes\n";
my %gene_reads;
get_gene_reads($cdna_end_1_bwtout);
get_gene_reads($cdna_end_2_bwtout);

# Filter fastq files for genes of interest
print "Filtering read sequences with alignments to fusion genes\n";
filter_read_sequences($paired_end_1_fastq,"1");
filter_read_sequences($paired_end_2_fastq,"2");

sub get_gene_reads
{
	my $bwtout_filename = $_[0];
	
	open BWTOUT, $bwtout_filename or die "Error: Unable to open bowtie alignment file $bwtout_filename: $!\n";
	while (<BWTOUT>)
	{
		chomp;
		my @align_info = split /\t/;
		
		my $readid = $align_info[0];
		my $ref = $align_info[2];
		my $seq = $align_info[4];
		my $qual = $align_info[5];
	
		# Split read name into id and end
		$readid =~ /^@?(.+)\/([12])$/ or print STDERR "Unable to interpret read id $readid\n" and next;
		my $read_pair_id = $1;
		my $read_pair_end = $2;

		# Interpret reference name as ensembl gene|transcript
		$ref =~ /^>?(ENSG\d+)\|(ENST\d+)/ or print STDERR "Unable to interpret reference name $ref\n" and next;
		my $gene = $1;
		my $transcript = $2;
		
		next unless $fusion_gene{$gene};

		push @{$gene_reads{$read_pair_id}}, $gene;
	}
	close BWTOUT;
}

sub filter_read_sequences
{
	my $fastq_filename = $_[0];
	my $read_end = $_[1];
	
	open FASTQ, $fastq_filename or die "Error: Unable to open fastq file $fastq_filename: $!\n";
	while (1)
	{	
		my $readid = <FASTQ>;
		my $sequence = <FASTQ>;
		my $comment = <FASTQ>;
		my $quality = <FASTQ>;
		
		last if not defined $quality;
		
		chomp($readid);
		chomp($sequence);
		chomp($comment);
		chomp($quality);
		
		# Split read name into id and end
		$readid =~ /^@?(.+)\/([12])$/ or print STDERR "Unable to interpret read id $readid\n" and next;
		my $read_pair_id = $1;
		my $read_pair_end = $2;
		
		if ($gene_reads{$read_pair_id})
		{
			foreach my $gene (@{$gene_reads{$read_pair_id}})
			{
				my $gene_fastq = $fusion_fastq_directory."/paired.$read_end.$gene.fastq";
				
				open GFQ, ">>".$gene_fastq or die "Error: Unable to open $gene_fastq\n";
				print GFQ $readid."\n";
				print GFQ $sequence."\n";
				print GFQ $comment."\n";
				print GFQ $quality."\n";
				close GFQ;
			}
		}
	}
	close FASTQ;
}


