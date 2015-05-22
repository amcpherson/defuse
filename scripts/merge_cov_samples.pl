#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use List::Util qw(min max);

die "Usage $0 samples1 samples2 ...\n" if scalar @ARGV == 0;

my %samples;

foreach my $samples_filename (@ARGV)
{
	open SMP, $samples_filename or die "Error: Unable to open $samples_filename\n";
	while (<SMP>)
	{
		chomp;
		my ($id,$sample) = split /\t/;
		push @{$samples{$id}}, $sample;
	}
	close SMP;
}

srand(11);

my $sum1 = 0.0;
my $sum2 = 0.0;
my $crossSum = 0.0;
my $count = 0.0;
foreach my $id (keys %samples)
{
	my $num_samples = scalar @{$samples{$id}};
	
	next if $num_samples < 2;
	
	my $sample1 = $samples{$id}->[int(rand($num_samples))];
	my $sample2 = $samples{$id}->[int(rand($num_samples))];
	
	$sum1 += $sample1;
	$sum2 += $sample2;
	$crossSum += $sample1 * $sample2;
	$count++;
}

$count > 100 or die "Error: not enough concordant read samples, set multi_exon_transcripts_stats = yes in config.txt";

my $mean = ($sum1 + $sum2) / (2.0 * $count);
my $cov = ($crossSum - $sum1 * $sum2 / $count) / $count;

print "mean\tcovariance\n";
print $mean."\t".$cov."\n";
