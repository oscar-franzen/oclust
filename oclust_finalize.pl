#!/usr/bin/perl

##################################################################################
#    Written by Oscar Franzen, <p.oscar.franzen@gmail.com>, Mt. Sinai, NY, USA.  #
#                                                                                #
#       "Det var inte vinsten jag ville åt, det var kampen.", A. Strindberg.     #
##################################################################################

use strict 'vars';
use Getopt::Long;
use Config;
use Cwd 'abs_path';

my $cwd = abs_path($0);

$cwd =~ /^(\S+\/)\S+/;
$cwd = $1;

sub workaround {
	my $p = abs_path($0);
	$p =~ /^(\S+\/)\S+/;
	$p = $1;

	my $path = $p . "modules";

	return $path;
}

use lib workaround();

use Bio::SeqIO;
use Memory::Usage;

# In this system equipped with LSF?
my $init = `bsub -V 2>&1`;

if ($init eq "") {
        print("This system is not supported. oclust can only run on a cluster environment equipped with the LSF scheduler.\n\n"); exit;
}

die("oclust is running on $Config{osname} ($Config{archname})\nWritten by Oscar Franzén <p.oscar.franzen\@gmail.com>, Mt. Sinai, New York

Command line arguments:
 -i [directory name]\t\tThe full path to the directory used as output in the previous step.
 -a [method]\t\t\tThe cluster algorithm to be used. Can be one of: complete (recommended), average, or single [default: complete]

 ") if (@ARGV == 0);

my $opt_dir;
my $opt_algo;

GetOptions ("i=s" => \$opt_dir,
	        "a=s" => \$opt_algo) or die("Error in command line arguments\n");

if ($opt_dir eq "") {
	print("Error: No output directory specified\n"); exit;
}

if (! -e $opt_dir) {
	print("Error: output directory does not exist\n"); exit;
}

if (! -e $opt_dir . "/submitted_jobs.txt") {
	print("Error: `oclust_finalize.pl' does not exist. You must run `oclust_finalize.pl' first.\n"); exit;
}

if ($opt_algo eq "") {
	$opt_algo = "complete";
}
elsif ($opt_algo ne "complete" && $opt_algo ne "average" && $opt_algo ne "single") {
	die("Invalid clustering algorithm specified, can be one of: complete (recommended), average or single");
}

# Check if any of the submitted jobs are pending or running
my $out = `bjobs | awk \'{print \$1}\'`;
my @items = split(/\n/,$out);

my %jobs_running;

foreach my $item (@items) {
	$jobs_running{$item} = 1;
}

my $found_running_jobs = 0;

open(fh, $opt_dir . "/submitted_jobs.txt");

while (my $line = <fh>) {
	chomp($line);

	if ($jobs_running{$line} ne "") {
		$found_running_jobs ++ ;
	}
}

close(fh);

if ($found_running_jobs > 0) {
	print("$found_running_jobs jobs are still pending or running. Please wait until these have finished.\n"); exit;
}

# Check if any of the submitted jobs failed
my @files = <$opt_dir/logs/*.log>;

my $found = 0;
my @l;

foreach my $file (@files) {
	my $out = `cat $file`;

	if ($out =~ /TERM/) {
		$found = 1;
		push(@l, $file);
	}
}

if ($found == 1) {
	print("Error: These jobs have failed:\n\n@l\n\n");
	print("Try increasing the runtime and/or memory, and rerun oclust_pipline.pl");

	exit;
}

my $mu = Memory::Usage->new();

# Record amount of memory used by current process
$mu->record('starting work');

# Create the dissimilarity matrix
my @files = <$opt_dir/alignments/*.pi>;

if (@files == 0) {
	print("No *.pi files found. Try re-running.\n\n"); exit;
}

my %matrix;

foreach my $file (@files) {
	open(fh, $file);

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+) (\S+) (\S+) /;
		my ($read1, $read2, $identity) = ($1,$2,$3);

		$matrix{$read1}{$read2} = $identity;
	}

	close(fh);
}

my @all_ids;

open(fh, $opt_dir."/targets.ss.FF.C.fa");

while (my $line = <fh>) {
	chomp($line);

	if ($line =~ /^>(\S+)/) {
		push(@all_ids, $1);
	}
}

close(fh);

open(fh, ">" . $opt_dir . "/matrix.mat");

my $first_row = " ";

foreach my $id (@all_ids) {
	$first_row .= "$id ";
}

chop($first_row);
print(fh "$first_row\n");

# Print the matrix
foreach my $id1 (@all_ids) {
	my $row_line = "$id1 ";
	
	foreach my $id2 (@all_ids) {
		my $pi = $matrix{$id1}{$id2};

		if ($id1 eq $id2) {
			$pi = 0;
		}

		if ($pi eq "") {
			print("Error: $id1 $id2\n"); exit;
		}

		$row_line .= $pi . " ";
	}

	chop($row_line);
	print(fh "$row_line\n");
}

print("Memory usage statistics:\n");
$mu->dump();

close(fh);

# Fix runscript
my $cmd = "cat $cwd/bin/R/bin/R";
my $o = `$cmd`;

$o =~ s/\/home\/foobar\/oclust\//$cwd/g;
$o =~ s/R_installed/R/g;

open(fh, ">" . $cwd."/bin/R/bin/R.fixed");
print(fh "$o\n");
close(fh);

my $cmd = "chmod +x $cwd/bin/R/bin/R.fixed";
`$cmd`;

print("Running hierarchical clustering using $opt_algo" ."-linkage hierarchical clustering\n");

my $cmd = "$cwd/bin/R/bin/R.fixed --no-save --args $opt_dir/matrix.mat $opt_dir/OTUs $opt_algo < $cwd/utils/hclust.R";

`$cmd`;

print("Finished successfully. OTU files:\n");

my $o = `ls -slht $opt_dir/*.hclust`;

print("$o\n");