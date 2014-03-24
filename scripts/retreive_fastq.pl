#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;

$| = 1;

use lib dirname($0);
use configdata;

my @usage;
push @usage, "Usage: $0 [options]\n";
push @usage, "Retrieve the fastq files.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -d, --data      Source Data Directory\n";
push @usage, "  -1,             Fastq file for end 1\n";
push @usage, "  -2,             Fastq file for end 2\n";
push @usage, "  -i, --index     Fastq index file\n";
push @usage, "  -n, --names     Fastq names file\n";
push @usage, "  -s, --sources   Fastq sources file\n";
push @usage, "  -f, --filter    Filter Garbage\n";

my $help;
my $config_filename;
my $source_directory;
my $reads_end_1_fastq;
my $reads_end_2_fastq;
my $reads_index_filename;
my $reads_names_filename;
my $reads_sources_filename;

GetOptions
(
	'help'        => \$help,
	'config=s'    => \$config_filename,
	'data=s'      => \$source_directory,
	'1=s'         => \$reads_end_1_fastq,
	'2=s'         => \$reads_end_2_fastq,
	'index=s'     => \$reads_index_filename,
	'names=s'     => \$reads_names_filename,
	'sources=s'   => \$reads_sources_filename,
);

not defined $help or die @usage;
defined $config_filename or die @usage;
defined $source_directory or die @usage;
defined $reads_end_1_fastq or die @usage;
defined $reads_end_2_fastq or die @usage;
defined $reads_index_filename or die @usage;
defined $reads_names_filename or die @usage;
defined $reads_sources_filename or die @usage;

my $config = configdata->new();
$config->read($config_filename);

# Config
my $scripts_directory = $config->get_value("scripts_directory");

-d $source_directory or die "Error: Source directory $source_directory does not exist.\n";

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
	}
	else
	{
		push @data_lane_regexs, $config->get_value("$data_lane_regex");
		push @data_end_regexs, "";
		push @data_compress_regexs, $config->get_value("$data_compress_regex");
	}

	if ($config->has_value("$data_converter"))
	{
		push @data_end1_converters, $config->get_value("$data_converter");
		push @data_end2_converters, $config->get_value("$data_converter");
	}
	else
	{
		push @data_end1_converters, $config->get_value("$data_end1_converter");
		push @data_end2_converters, $config->get_value("$data_end2_converter");
	}
	
	$converter_num++;
}

# Look for a previous source list
my @previous_source_list;
if (-e $reads_sources_filename)
{
	open RSRC, $reads_sources_filename or die "Error: Unable to read $reads_sources_filename: $!\n";
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
		print "Ignoring unconvertable file $raw_filename\n";
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
			print "Ignoring unpaired file $filename\n";
			delete $convert_filename{$lane};
			delete $convert_command{$lane};
			last;
		}
	}
}

# Create a sorted list of filenames
my @source_list = sort {$a cmp $b} keys %convert_filename;

# Check that we have at least one source file
die "Error: No source data files found\n" if scalar @source_list == 0;

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
	print "No new source data files found\n";
	exit;
}

# Open both fastq files for simultaneous writing
open FQ1, ">".$reads_end_1_fastq or die "Error: Unable to write to $reads_end_1_fastq: $!\n";
open FQ2, ">".$reads_end_2_fastq or die "Error: Unable to write to $reads_end_2_fastq: $!\n";
open FQI, ">".$reads_index_filename or die "Error: Unable to open $reads_index_filename\n"; binmode(FQI);
open NAM, ">".$reads_names_filename or die "Error: Unable to open $reads_names_filename\n";

# Iterate through each lane and add fastq sequences
my $current_fragment_index = 0;
foreach my $lane (keys %convert_command)
{
	die "Error: lane $lane is unpaired\n" if not defined $convert_command{$lane}{'1'} or not defined $convert_command{$lane}{'2'};
	
	print "Converting ".$convert_filename{$lane}{'1'}." and ".$convert_filename{$lane}{'2'}."\n";
	
	open IN1, $convert_command{$lane}{"1"}." |" or die "Error: Unable to start command ".$convert_command{$lane}{"1"}."\n";
	open IN2, $convert_command{$lane}{"2"}." |" or die "Error: Unable to start command ".$convert_command{$lane}{"2"}."\n";
	
	while (1)
	{
		my $readid1 = <IN1>;
		my $sequence1 = <IN1>;
		my $comment1 = <IN1>;
		my $quality1 = <IN1>;

		last if not defined $quality1;

		chomp($readid1);
		chomp($sequence1);
		chomp($comment1);
		chomp($quality1);

		my $readid2 = <IN2>;
		my $sequence2 = <IN2>;
		my $comment2 = <IN2>;
		my $quality2 = <IN2>;

		last if not defined $quality2;

		chomp($readid2);
		chomp($sequence2);
		chomp($comment2);
		chomp($quality2);
		
		my $filepos1 = pack('q',tell(FQ1));
		my $filepos2 = pack('q',tell(FQ2));
		
		print FQI $filepos1;
		print FQI $filepos2;
		print FQ1 "\@$current_fragment_index/1\n$sequence1\n$comment1\n$quality1\n";
		print FQ2 "\@$current_fragment_index/2\n$sequence2\n$comment2\n$quality2\n";
		print NAM "$current_fragment_index\t$readid1\t$readid2\n";
		
		$current_fragment_index++;
	}
	
	close IN1;
	my $retcode1 = $? >> 8;
	
	close IN2;
	my $retcode2 = $? >> 8;

	if ($retcode1 != 0)
	{
		print STDERR "Error: the following command failed:\n".$convert_command{$lane}{"1"}."\n";
	}
	
	if ($retcode2 != 0)
	{
		print STDERR "Error: the following command failed:\n".$convert_command{$lane}{"2"}."\n";
	}
	
	if ($retcode1 != 0 or $retcode2 != 0)
	{
		die "Conversion Failed\n";
	}
}

# Check if we found reads
die "Error: No reads found\n" if $current_fragment_index == 0;
print "Imported $current_fragment_index reads\n";

# Finished writing fastq sequences
close FQ1;
close FQ2;
close FQI;
close NAM;

# Write out source list
my $reads_sources_temp = $reads_sources_filename.".tmp";
open RSRC, ">".$reads_sources_temp or die "Error: Unable to write to source list $reads_sources_temp: $!\n";
foreach my $source_filename (@source_list)
{
	print RSRC $source_filename."\n";
}
close RSRC;
rename $reads_sources_temp, $reads_sources_filename;

print "Finished Retrieval\n";

