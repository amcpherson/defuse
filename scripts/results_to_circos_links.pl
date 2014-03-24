#!/usr/bin/perl

use strict;
use warnings;


my %found_fusions;
my $found_header = 0;
my %ind;
while (<>)
{
	chomp;
	my @fields = split /\t/;

	if (not $found_header)
	{
		foreach my $field_index (0..$#fields)
		{
			$ind{$fields[$field_index]} = $field_index;
		}

		$found_header = 1;

		next;
	}

	if ($found_header and $fields[0] eq "cluster_id")
	{
		next;
	}	

	my $id = $fields[$ind{cluster_id}];
	my $libname = $fields[$ind{library_name}];
	my $chr1 = $fields[$ind{gene_chromosome1}];
	my $chr2 = $fields[$ind{gene_chromosome2}];
	my $pos1 = $fields[$ind{genomic_break_pos1}];
	my $pos2 = $fields[$ind{genomic_break_pos2}];

	print $libname."-".$id."\ths".$chr1."\t".$pos1."\t".$pos1."\n";
	print $libname."-".$id."\ths".$chr2."\t".$pos2."\t".$pos2."\n";
}

