#!/usr/bin/perl


my $mate = shift;

usage_exit() if not defined $mate or $mate ne '1' and $mate ne '2';

sub usage_exit
{
	print "Usage: $0 [1|2]\n";
	exit;
}

while (<>)
{
	my $line = $_;

	chomp;
	next if /^#/;

	my @sam_info = split /\t/;
	
	my $qname = $sam_info[0];
	my $flag = $sam_info[1];
	my $rname = $sam_info[2];
	my $pos = $sam_info[3];
	my $mapq = $sam_info[4];
	my $cigar = $sam_info[5];
	my $mrnm = $sam_info[6];
	my $mpos = $sam_info[7];
	my $isize = $sam_info[8];
	my $seq = $sam_info[9];
	my $qual = $sam_info[10];
	my $opt = $sam_info[11];

	my $read_pair_end;
	$read_pair_end = '1' if $flag & hex('0x0040');
	$read_pair_end = '2' if $flag & hex('0x0080');
	die "SAM Format Error\n" if not defined $read_pair_end;

	print $line if $read_pair_end eq $mate;
}
