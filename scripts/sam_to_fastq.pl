#!/usr/bin/perl


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
	$read_pair_end = '/1' if $flag & hex('0x0040');
	$read_pair_end = '/2' if $flag & hex('0x0080');

	$read_pair_end = "" if not defined $read_pair_end;

	my $strand = '+';
	$strand = '-' if $flag & hex('0x0010');

	if ($strand eq '-')
	{
		$seq = reverse($seq);
		$seq =~ tr/ACGTacgt/TGCAtgca/;
		$qual = reverse($qual);
	}

	print "\@$qname$read_pair_end\n$seq\n+\n$qual\n";
}
