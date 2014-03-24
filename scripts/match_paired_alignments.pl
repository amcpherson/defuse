#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 end_1_sam end_2_sam > paired_sam\n";

my $end_1_sam_filename = shift;
my $end_2_sam_filename = shift;

defined $end_2_sam_filename or die $usage;

open END1, $end_1_sam_filename or die;
open END2, $end_2_sam_filename or die;

my $line1 = <END1>;
while (defined $line1 and $line1 =~ /^\@/)
{
	$line1 = <END1>;
}

my $line2 = <END2>;
while (defined $line2 and $line2 =~ /^\@/)
{
	$line2 = <END2>;
}

sub get_fragment_id
{
	my $line = shift;
	my @sam_info = split /\t/, $line;
	$sam_info[0] =~ /(.*)\/([12])/;
	return $1;
}

while (defined $line1 or defined $line2)
{
	if (not defined $line1)
	{
		print $line2;
		$line2 = <END2>;
		next;
	}
	
	if (not defined $line2)
	{
		print $line1;
		$line1 = <END1>;
		next;
	}
	
	my $fragment_id1 = get_fragment_id($line1);
	my $fragment_id2 = get_fragment_id($line2);
	
	if ($fragment_id1 <= $fragment_id2)
	{
		print $line1;
		$line1 = <END1>;
	}
	else
	{
		print $line2;
		$line2 = <END2>;
	}
}

close END1;
close END2;