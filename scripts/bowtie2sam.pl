#!/usr/bin/perl -w

# Contact: lh3
# Version: 0.1.1

use strict;
use warnings FATAL => 'all';
use Getopt::Std;

&bowtie2sam;
exit;

sub bowtie2sam {
  my %opts = ();
	  die("Usage: bowtie2sam.pl <aln.bowtie>\n") if (@ARGV == 0 && -t STDIN);
  # core loop
  my @s;
  while (<>) {
	my ($name, $nm) = &bowtie2sam_aux($_, \@s); # read_name, number of mismatches
	print join("\t", @s), "\n";
  }
}

sub bowtie2sam_aux {
  my ($line, $s) = @_;
  chomp($line);
  my @t = split("\t", $line);
  my $ret;
  @$s = ();
  # read name
  $s->[0] = $ret = $t[0];
  $s->[0] =~ s/\/[12]$//g;
  # initial flag (will be updated later)
  $s->[1] = 0;
  # read & quality
  $s->[9] = $t[4]; $s->[10] = $t[5];
  # cigar
  $s->[5] = length($s->[9]) . "M";
  # coor
  $s->[2] = $t[2]; $s->[3] = $t[3] + 1;
  $s->[1] |= 0x10 if ($t[1] eq '-');
  # mapQ
  $s->[4] = $t[6] == 0? 25 : 0;
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  my $nm = @t - 7;
  push(@$s, "NM:i:" . (@t-7));
  push(@$s, "X$nm:i:" . ($t[6]+1));
  my $md = '';
  if ($t[7]) {
	$_ = $t[7];
	my $a = 0;
	while (/(\d+):[ACGTN]>([ACGTN])/gi) {
	  my ($y, $z) = ($1, $2);
	  $md .= (int($y)-$a) . $z;
	  $a += $y - $a + 1;
	}
	$md .= length($s->[9]) - $a;
  } else {
	$md = length($s->[9]);
  }
  push(@$s, "MD:Z:$md");
  return ($ret, $nm);
}
