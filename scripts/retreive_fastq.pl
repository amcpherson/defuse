#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];

use lib dirname($0);
use configdata;
use cmdrunner;

my @usage;
push @usage, "Usage: $0 [options]\n";
push @usage, "Retrieve the fastq files.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --data      Source Data Directory\n";
push @usage, "  -o, --output    Output Directory\n";
push @usage, "  -s, --submit    Submitter Type\n";

my $help;
my $config_filename;
my $source_directory;
my $output_directory;
my $submitter_type;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'data=s'      => \$source_directory,
	'output=s'    => \$output_directory,
	'submit=s'    => \$submitter_type,
);

not defined $help or die @usage;
defined $config_filename or die @usage;
defined $source_directory or die @usage;
defined $output_directory or die @usage;
defined $submitter_type or die @usage;

my $config = configdata->new();
$config->read($config_filename);

# Config
my $scripts_directory = $config->get_value("scripts_directory");

# Script paths
my $filter_fastq_garbage_script = "$scripts_directory/filter_paired_fastq_garbage.pl";
my $index_fastq_script = "$scripts_directory/index_paired_fastq.pl";

$source_directory = abs_path($source_directory);
$output_directory = abs_path($output_directory);

-e $source_directory or die "Error: Source directory $source_directory does not exist.\n";
-e $output_directory or die "Error: Output directory $output_directory does not exist.\n";

my $log_directory = $output_directory."/log";
my $log_prefix = $log_directory."/retrieve_fastq";

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("retrieve_fastq");
$runner->prefix($log_prefix);
$runner->submitter($submitter_type);

# Read in all data regex's and conversion scripts
my @data_lane_regexs;
my @data_end_regexs;
my @data_compress_regexs;
my @data_end1_converters;
my @data_end2_converters;
my $converter_num = 1;
while (1)
{
	my $data_lane_regex = "data_lane_regex_$converter_num";
	my $data_end_regex = "data_end_regex_$converter_num";
	my $data_compress_regex = "data_compress_regex_$converter_num";
	my $data_converter = "data_converter_$converter_num";
	my $data_end1_converter = "data_end1_converter_$converter_num";
	my $data_end2_converter = "data_end2_converter_$converter_num";
	
	if (not $config->has_value($data_lane_regex))
	{
		last;
	}
	
	if ($config->has_value($data_end_regex))
	{
		push @data_lane_regexs, $config->get_value("$data_lane_regex");
		push @data_end_regexs, $config->get_value("$data_end_regex");
		push @data_compress_regexs, $config->get_value("$data_compress_regex");
		push @data_end1_converters, $config->get_value("$data_converter");
		push @data_end2_converters, $config->get_value("$data_converter");
	}
	else
	{
		push @data_lane_regexs, $config->get_value("$data_lane_regex");
		push @data_end_regexs, "";
		push @data_compress_regexs, $config->get_value("$data_compress_regex");
		push @data_end1_converters, $config->get_value("$data_end1_converter");
		push @data_end2_converters, $config->get_value("$data_end2_converter");
	}
	
	$converter_num++;
}

my $reads_end_1_fastq = $output_directory."/reads.1.fastq";
my $reads_end_2_fastq = $output_directory."/reads.2.fastq";
my $reads_sources = $output_directory."/reads.sources";
my $reads_prefix = $output_directory."/reads";

# Look for a previous source list
my @previous_source_list;
if (-e $reads_sources)
{
	open RSRC, $reads_sources or die "Error: Unable to read $reads_prefix: $!\n";
	@previous_source_list = <RSRC>;
	close RSRC;
	
	print "Found previous source list including the following files\n";
	print @previous_source_list;
	chomp(@previous_source_list);
}
else
{
	print "No previous source list found\n";
}

print "Seaching for source data files\n";

# Find all data files
my @raw_listings = <$source_directory/*>;

# Find all filenames and conversion commands
my %convert_filename;
my %convert_command;
foreach my $raw_listing (@raw_listings)
{
	my $raw_filename = fileparse($raw_listing);
	
	my $convertable = 0;
	foreach my $converter_index (0..$#data_lane_regexs)
	{
		my $data_lane_regex = $data_lane_regexs[$converter_index];
		my $data_end_regex = $data_end_regexs[$converter_index];
		my $data_compress_regex = $data_compress_regexs[$converter_index];

		my %data_converter;
		$data_converter{'1'} = $data_end1_converters[$converter_index];
		$data_converter{'2'} = $data_end2_converters[$converter_index];
		
		if ($raw_filename =~ /$data_lane_regex/)
		{
			my $lane = $1;
			
			$raw_filename =~ /$data_compress_regex/;
			my $ext = $1;
			
			defined $ext or die "Error: Unable to determine extension for file $raw_filename\n";
			
			my $decompress_command;
			if ($ext eq ".bz2")
			{
				$decompress_command = "bzcat $raw_listing ";
			}
			elsif ($ext eq ".gz")
			{
				$decompress_command = "gunzip -c $raw_listing ";
			}
			else
			{
				$decompress_command = "cat $raw_listing ";
			}
			
			if ($data_end_regex ne "")
			{
				$raw_filename =~ /$data_end_regex/;
				my $end = $1;
				
				defined $end or die "Error: Unable to determine end for file $raw_filename\n";
				$end eq "1" or $end eq "2" or die "Error: Invalid end for file $raw_filename\n";
				
				$convert_command{$lane}{$end} = $decompress_command."| ".$data_converter{$end};
				$convert_filename{$lane}{$end} = $raw_filename;
			}
			else
			{
				$convert_command{$lane}{'1'} = $decompress_command."| ".$data_converter{'1'};
				$convert_filename{$lane}{'1'} = $raw_filename;
				
				$convert_command{$lane}{'2'} = $decompress_command."| ".$data_converter{'2'};
				$convert_filename{$lane}{'2'} = $raw_filename;
			}
			
			$convertable = 1;
			last;
		}
	}
	
	if (not $convertable)
	{
		print "Found unconvertable file $raw_filename\n";
	}
}

# Remove unpaired filenames and commands
foreach my $lane (keys %convert_filename)
{
	foreach my $end ('1','2')
	{
		my $other_end = ($end eq '1') ? '2' : '1';

		my $filename = $convert_filename{$lane}{$end};
		
		if (not defined $convert_filename{$lane}{$other_end})
		{
			print "Found unpaired file $filename\n";
			delete $convert_filename{$lane}{$end};
		}
	}
}

# Create a sorted list of filenames
my @source_list = sort {$a cmp $b} keys %convert_filename;

# Compare with existing source file list
my $same_source_files = 1;
if (scalar @source_list != scalar @previous_source_list)
{
	$same_source_files = 0;
}
else
{
	foreach my $filename_index (0..$#source_list)
	{
		if ($source_list[$filename_index] ne $previous_source_list[$filename_index])
		{
			$same_source_files = 0;
			last;
		}
	}
}

# Exit if these are the same source files and the fastq files exist
if ($same_source_files and -e $reads_end_1_fastq and -e $reads_end_2_fastq)
{
	print "No new data files found\n";
	exit;
}

# Remove previous fastq files and sources file
unlink $reads_end_1_fastq;
unlink $reads_end_2_fastq;
unlink $reads_sources;

# Fastq files hashed by end
my %fastq_end_filenames;
$fastq_end_filenames{'1'} = $reads_end_1_fastq;
$fastq_end_filenames{'2'} = $reads_end_2_fastq;	

# Run conversion command for each source file
foreach my $lane (keys %convert_command)
{
	die if not defined $convert_command{$lane}{'1'} or not defined $convert_command{$lane}{'2'};
	
	foreach my $end ('1','2')
	{
		print "Converting raw file $convert_filename{$lane}{$end} containing $end end for $lane\n";

		$runner->run("$convert_command{$lane}{$end} >> $fastq_end_filenames{$end}", [], []);
	}
}

# Heuristic filter of reads that are likely garbage
print "Filtering garbage\n";
$runner->run("$filter_fastq_garbage_script $reads_prefix", [], []);

# Index the fastq files
print "Indexing\n";
$runner->run("$index_fastq_script $reads_prefix", [], []);

# Write out source list
my $reads_sources_temp = $reads_sources.".tmp";
open RSRC, ">".$reads_sources_temp or die "Error: Unable to write to source list $reads_sources_temp: $!\n";
foreach my $source_filename (@source_list)
{
	print RSRC $source_filename."\n";
}
close RSRC;
rename $reads_sources_temp, $reads_sources;

print "Finished Retrieval\n";

