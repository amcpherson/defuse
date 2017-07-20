#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use File::Spec;
use List::Util qw[min max];
use Cwd qw[abs_path];

use FindBin;
use lib "$FindBin::RealBin";
use configdata;
use parsers;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Retrieve reads supporting a fusion.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --dataset   Dataset Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -i, --id        Cluster ID\n";

my $help;
my $config_filename;
my $dataset_directory;
my $output_directory;
my $query_cluster_id;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'dataset=s'   => \$dataset_directory,
	'output=s'    => \$output_directory,
	'id=s'        => \$query_cluster_id,
);

not defined $help or die @usage;

defined $dataset_directory or die @usage;
defined $output_directory or die @usage;

my $source_directory = abs_path("$FindBin::RealBin/../");

if (not defined $config_filename)
{
	$config_filename = $source_directory."/scripts/config.txt";
}

-e $config_filename or die "Error: Unable to find config file $config_filename\n";

my $config = configdata->new();
$config->read($config_filename, $dataset_directory, $source_directory);

# Config values
my $reference_fasta = $config->get_value("reference_fasta");
my $cdna_regions = $config->get_value("cdna_regions");
my $scripts_directory = $config->get_value("scripts_directory");
my $tools_directory = $config->get_value("tools_directory");
my $samtools_bin = $config->get_value("samtools_bin");

my $splitseq_bin = "$tools_directory/splitseq";

my $reads_prefix = $output_directory."/reads";
my $read_stats = $output_directory."/concordant.read.stats";
my $clusters_sc = $output_directory."/clusters.sc";
my $reads_split_catalog = $output_directory."/reads.split.catalog";

my @split_fastq_prefixes = readsplitcatalog($reads_split_catalog);

my %read_stat_values;
parsers::get_stats($read_stats, \%read_stat_values);

my $read_length_min = $read_stat_values{"readlength_min"};
my $read_length_max = $read_stat_values{"readlength_max"};
my $fragment_mean = $read_stat_values{"fraglength_mean"};
my $fragment_stddev = $read_stat_values{"fraglength_stddev"};

my %cluster_info;
read_clusters($clusters_sc, $query_cluster_id, \%cluster_info);

if (not defined $cluster_info{"0"}{ref_name})
{
	die "Unable to find cluster $query_cluster_id\n";
}

print "Split Reads:\n";
my $clusters_sc_regions = $output_directory."/clusters.sc.regions";
my $splitreads_predalign = $output_directory."/splitreads.predalign";
system("$splitseq_bin -r $clusters_sc_regions -n $read_length_min -x $read_length_max -u $fragment_mean -s $fragment_stddev -e $cdna_regions -f $reference_fasta -p $reads_prefix -a $splitreads_predalign -i $query_cluster_id") == 0 or die;
print "\n";

print "Spanning Reads:\n";
my %cluster_lines;
my %found_fragments;
foreach my $split_fastq_prefix (@split_fastq_prefixes)
{
	my $spanning_filelist = $split_fastq_prefix.".spanning.filelist";
	open SFL, $spanning_filelist or die;
	while (<SFL>)
	{
		chomp;
		my ($chr1,$chr2,$filename) = split /\t/;
		
		open SAL, $filename or die "Error: Unable to open $filename: $!\n";
		while (<SAL>)
		{
			my $line = $_;
			
			chomp;
			my @fields = split /\t/;
			
			my $fragment_index = $fields[0];
			my $read_end = $fields[1];
			my $rname = $fields[2];
			my $strand = $fields[3];
			my $start = $fields[4];
			my $end = $fields[5];
			
			foreach my $cluster_end ("0","1")
			{
				if (defined $cluster_info{$cluster_end}{fragments}{$fragment_index})
				{
					if ($strand eq $cluster_info{$cluster_end}{strand} and $start <= $cluster_info{$cluster_end}{end} and $end >= $cluster_info{$cluster_end}{start})
					{
						$cluster_lines{$fragment_index}{$cluster_end} = $line;
						$found_fragments{$fragment_index} = 1;
					}
				}
			}
		}
		close SAL;
	}
	close SFL;
}

foreach my $fragment_index (sort {$a <=> $b} keys %found_fragments)
{
	foreach my $cluster_end ("0","1")
	{
		die "Error: Could not match reads\n" if not defined $cluster_lines{$fragment_index}{$cluster_end};
		print $cluster_lines{$fragment_index}{$cluster_end};
	}
}

sub readsplitcatalog
{
        my $split_catalog = shift;
        my @split_prefixes;

        open SC, $split_catalog or die "Error: Unable to open $split_catalog: $!\n";
        while (<SC>)
        {
                chomp;
                my @fields = split /\t/;

                push @split_prefixes, $fields[0];
        }
        close SC;

        return @split_prefixes;
}

sub read_clusters
{
	my $clusters_filename = shift;
	my $query_cluster_id = shift;
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
		
		next unless $cluster_id == $query_cluster_id;

		$clusters_hash_ref->{$cluster_end}{ref_name} = $ref_name;
		$clusters_hash_ref->{$cluster_end}{strand} = $strand;

		$clusters_hash_ref->{$cluster_end}{start} = $start if not defined $clusters_hash_ref->{$cluster_end}{start};
		$clusters_hash_ref->{$cluster_end}{end} = $end if not defined $clusters_hash_ref->{$cluster_end}{end};

		$clusters_hash_ref->{$cluster_end}{start} = min($start, $clusters_hash_ref->{$cluster_end}{start});
		$clusters_hash_ref->{$cluster_end}{end} = max($end, $clusters_hash_ref->{$cluster_end}{end});
		
		$clusters_hash_ref->{$cluster_end}{fragments}{$fragment_id} = 1;
	}
	close CLU;
}

