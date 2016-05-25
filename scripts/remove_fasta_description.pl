#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

use FindBin;
use lib "$FindBin::RealBin/../external/BioPerl-1.6.1";
use Bio::SeqIO;

my $in = Bio::SeqIO->new(-format => "fasta", -fh => \*STDIN);
my $out = Bio::SeqIO->new(-format => "fasta", -fh => \*STDOUT);

while (my $seq = $in->next_seq())
{
	$seq->description("");
	$out->write_seq($seq);
}

