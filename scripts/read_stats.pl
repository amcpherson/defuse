#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use List::Util qw(min max);

my $fraglength_sum = 0;
my $fraglength_sum_sq = 0;
my $fraglength_num = 0;

my %readlengths;

my $line = 0;
while (my $sam_line = <>)
{
	$line++;

	next if $sam_line =~ /^\@/;

	chomp $sam_line;
	my @sam_info_1 = split /\t/, $sam_line;

	$sam_line = <>;
	$line++;
	chomp $sam_line;
	my @sam_info_2 = split /\t/, $sam_line;
	
	my $read_name_1 = $sam_info_1[0];
	my $read_name_2 = $sam_info_2[0];

	$read_name_1 =~ s/\/[12]//;
	$read_name_2 =~ s/\/[12]//;

	$read_name_1 eq $read_name_2 or die "Error: Sam file error at line $line\n";

	my $flag_1 = $sam_info_1[1];
	my $flag_2 = $sam_info_2[1];

	next if not $flag_1 & hex('0x0002');
	next if not $flag_2 & hex('0x0002');

	my $reference_1 = $sam_info_1[2];
	my $reference_2 = $sam_info_2[2];

	my $isize_1 = abs($sam_info_1[8]);
	my $isize_2 = abs($sam_info_2[8]);

	$isize_1 == $isize_2 or die "Error: Sam file isize error at line $line\n";

	my $sequence_1 = $sam_info_1[9];
	my $sequence_2 = $sam_info_2[9];

	$reference_1 eq $reference_2 or warn "Warning: Sam file includes mismatched pair at line $line\n" and next;

	my $seq_length_1 = length($sequence_1);
	my $seq_length_2 = length($sequence_2);
	
	my $fraglength = $isize_1;

	$fraglength_sum += $fraglength;
	$fraglength_sum_sq += $fraglength ** 2;
	$fraglength_num++;
	
	$readlengths{$seq_length_1} = 1;
	$readlengths{$seq_length_2} = 1;
}

my $readlength_min = min(keys %readlengths);
my $readlength_max = max(keys %readlengths);

$readlength_min = 0 if not defined $readlength_min;
$readlength_max = 0 if not defined $readlength_max;

my $fraglength_mean = 0;
my $fraglength_variance = 0;
my $fraglength_stddev = 0;

if ($fraglength_num > 0)
{
	$fraglength_mean = $fraglength_sum / $fraglength_num;
	$fraglength_variance = $fraglength_sum_sq / $fraglength_num - $fraglength_mean ** 2;
	$fraglength_stddev = $fraglength_variance ** 0.5;
}

my $readlengths_list = join ",", keys %readlengths;

print "frag_count\tfraglength_mean\tfraglength_stddev\treadlength_min\treadlength_max\treadlengths_list\n";
print "$fraglength_num\t$fraglength_mean\t$fraglength_stddev\t$readlength_min\t$readlength_max\t$readlengths_list\n";



