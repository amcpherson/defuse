#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';

my $usage = "Usage: $0 paired_sam > paired_matched_sam\n";

sub get_fragment_id_read_end
{
	my $line = shift;
	my @sam_info = split /\t/, $line;
	$sam_info[0] =~ /(.*)\/([12])/;
	return ($1,$2);
}

my $current_fragment_id;
my %current_read_ends;
my @current_lines;
while (my $line = <>)
{
	my @sam_info = split /\t/, $line;
	
	$sam_info[0] =~ /(.*)\/([12])/;
	my $fragment_id = $1;
	my $read_end = $2;
	
	if (defined $current_fragment_id and $fragment_id != $current_fragment_id)
	{
		if (keys %current_read_ends == 2)
		{
			print @current_lines;
		}
		
		%current_read_ends = ();
		@current_lines = ();
	}
	
	$current_fragment_id = $fragment_id;
	$current_read_ends{$read_end} = 1;
	push @current_lines, $line;
}

if (defined $current_fragment_id)
{
	if (keys %current_read_ends == 2)
	{
		print @current_lines;
	}
}

