#!/usr/bin/perl

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];

use lib dirname($0);
use configdata;
use cmdrunner;

$| = 1;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "Split fasta and align with blat.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -r, --ref       Reference 2bit filename\n";
push @usage, "  -c, --chr       Chromosome Names, comma separated\n";
push @usage, "  -f, --fasta     Sequences Fasta\n";
push @usage, "  -o, --out       Output PSL\n";
push @usage, "  -s, --submit    Submitter Type (default: direct)\n";
push @usage, "  -p, --parallel  Maximum Number of Parallel Jobs (default: 1)\n";
push @usage, "  -b, --blat      Blat binary path (default: blat)\n";
push @usage, "  -n, --num       Number of sequences per split (default: 10000)\n";

my $help;
my $reference_2bit;
my $chromosomes;
my $sequences_fasta;
my $output_psl;
my $submitter_type;
my $max_parallel;
my $blat_bin;
my $num_sequences;

GetOptions
(
	'help'        => \$help,
	'ref=s'       => \$reference_2bit,
	'chr=s'       => \$chromosomes,
	'fasta=s'     => \$sequences_fasta,
	'out=s'       => \$output_psl,
	'submit=s'    => \$submitter_type,
	'parallel=s'  => \$max_parallel,
	'blat=s'      => \$blat_bin,
	'num=s'       => \$num_sequences,
);

not defined $help or die @usage;

defined $reference_2bit or die @usage;
defined $chromosomes or die @usage;
defined $sequences_fasta or die @usage;
defined $output_psl or die @usage;

$submitter_type = "direct" if not defined $submitter_type;
$max_parallel = 1 if $submitter_type eq "direct" and not defined $max_parallel;
$max_parallel = 200 if $submitter_type ne "direct" and not defined $max_parallel;
$blat_bin = "blat" if not defined $blat_bin;
$num_sequences = 10000 if not defined $num_sequences;

-e $reference_2bit or die "Error: Unable to find file $reference_2bit\n";
-e $sequences_fasta or die "Error: Unable to find file $sequences_fasta\n";

$reference_2bit = abs_path($reference_2bit);
$sequences_fasta = abs_path($sequences_fasta);
$output_psl = abs_path($output_psl);

my $script_directory = abs_path(dirname($0));
my $split_fasta_script = $script_directory."/split_fasta.pl";

my $log_directory = dirname($sequences_fasta)."/log";
my $log_prefix = $log_directory."/splitblatmerge.".basename($reference_2bit);

mkdir $log_directory if not -e $log_directory;

my $runner = cmdrunner->new();
$runner->name("splitblatmerge.".basename($reference_2bit));
$runner->prefix($log_prefix);
$runner->submitter($submitter_type);
$runner->maxparallel($max_parallel);
$runner->jobmem(3000000000);

print "Splitting fasta file\n";
my $split_prefix = $sequences_fasta;
my $split_catalog = $sequences_fasta.".".basename($reference_2bit).".split.catalog";
$runner->run("$split_fasta_script #<1 $num_sequences $split_prefix > #>1", [$sequences_fasta], [$split_catalog]);

my @split_filenames = `cat $split_catalog`;
chomp(@split_filenames);

print "Running blat\n";
my @split_chromosome_psls;
foreach my $chromosome (split /,/, $chromosomes)
{
	foreach my $split_filename (@split_filenames)
	{
		my $split_chromosome_psl = $split_filename.".".$chromosome.".psl";
		push @split_chromosome_psls, $split_chromosome_psl;
		
		$runner->padd("$blat_bin $reference_2bit:$chromosome #<1 #>1", [$split_filename], [$split_chromosome_psl]);
	}
}
$runner->prun();

print "Merging psl\n";
my $output_psl_tmp = $output_psl.".tmp";
unlink $output_psl_tmp;
foreach my $psl_filename (@split_chromosome_psls)
{
	$runner->run("cat $psl_filename >> $output_psl_tmp", [], []);
}
rename $output_psl_tmp, $output_psl;


