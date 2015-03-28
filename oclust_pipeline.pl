#!/usr/bin/perl

###################################################################
# Oscar Franzen, <p.oscar.franzen@gmail.com>, Mt. Sinai, NY, USA. #
###################################################################

use strict 'vars';
use Getopt::Long;
use Config;
use POSIX;
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

use lib workaround(); # path to bioperl
use Bio::SeqIO;

# Can we find the binaries?
if (! -e "$cwd/bin/blastall") {
	print("Error: blastall binary not found.\n\n"); exit;
}

if (! -e "$cwd/bin/megablast") {
	print("Error: megablast binary not found.\n\n"); exit;
}

if (! -e "$cwd/bin/formatdb") {
	print("Error: formatdb binary not found.\n\n"); exit;
}

if (! -e "$cwd/bin/emboss.needle.static.linux.x86-64") {
	print("Error: emboss.needle.static.linux.x86-64 binary not found.\n\n"); exit;
}

if (! -e "$cwd/bin/uchime4.2.40_i86linux32") {
	print("Error: uchime4.2.40_i86linux32 binary not found.\n\n"); exit;
}

if (! -e "$cwd/bin/hmmscan") {
	print("Error: hmmscan binary not found.\n\n"); exit;
}

die("oclust is running on $Config{osname} ($Config{archname})\nFeedback: <p.oscar.franzen\@gmail.com>, Mount Sinai, New York, U.S.A.\n\nCommand line arguments:

	-x <method> -f <input fasta file> -o <output directory> -p 1 -minl 400 -maxl 1000

	General settings:
	-x NW or MSA               Can be NW for Needleman-Wunsch or MSA for Infernal. [MSA]
	-f [string]                Input fasta file.
	-o [string]                Name of output directory (must not exist) and use full path.
	-R HMM, BLAST, or none     Method to use for reverse complementing sequences. [HMM]
	-p [integer]               If -R is BLAST: Number of processor cores to use. [4]
	-minl [integer]            Minimum sequence length. [optional]
	-maxl [integer]            Maximum sequence length. [optional]
	-rand [integer]            Randomly sample a specified number of sequences. [optional]
	-human Y or N              If 'Y'es, then execute BLAST-based contamination
	                           screen towards the human genome. [Y]
	-chimera Y or N            Run chimera check. Can be Y or N. [Y]

	LSF settings:
	-lsf_queue [string]        Name of the LSF queue to use [scavenger].
	-lsf_account [string]      Name of the account to use.
	-lsf_time [integer]        Runtime hours per job specified as number of hours. [12]
	-lsf_memory [integer]      Requested amount of RAM in MB. [20000]
	-lsf_nb_jobs [integer]     Number of jobs. [20]

	Usage example: -f /home/foobar/long_reads.fasta -o /home/foobar/foo -p 4 -minl 700 -maxl 800\n\n") if (@ARGV == 0);

my $opt_f;
my $opt_o;
my $opt_p;
my $opt_min_length;
my $opt_max_length;
my $debug;
my $random_subset;
my $human;
my $chimera;
my $revcom_method;
my $distance;

my $lsf_nb_jobs;
my $lsf_queue;
my $lsf_time;
my $lsf_memory;
my $lsf_account;

GetOptions ("file=s" => \$opt_f,
	         "out=s" => \$opt_o,
	         "proc=i" => \$opt_p,
	         "minl=i" => \$opt_min_length,
	         "maxl=i" => \$opt_max_length,
	         "debug=i" => \$debug,
	         "rand=i" => \$random_subset,
	         "human=s" => \$human,
	         "chimera=s" => \$chimera,
	         "lsf_nb_jobs=i" => \$lsf_nb_jobs,
	         "lsf_queue=s" => \$lsf_queue,
	         "lsf_account=s" => \$lsf_account,
	         "lsf_time=i" => \$lsf_time,
	         "lsf_memory=i" => \$lsf_memory,
	         "lsf_account=s" => \$lsf_account,
	         "R=s" => \$revcom_method,
	         "x=s" => \$distance) or die("Error in command line arguments\n");

if ($opt_f eq "") {
	print("Error: No fasta input file specified\n"); exit;
}
elsif ($opt_o eq "") {
	print("Error: No output directory specified\n"); exit;
}
elsif ($opt_p eq "") {
	$opt_p = 4;
}

if ($chimera eq "") {
	$chimera = "Y";
}

if ($chimera !~ /^[YN]$/) {
	die("-chimera can take options Y or N.\n");
}

if ($lsf_nb_jobs eq "") {
	$lsf_nb_jobs = 20;
}

if ($lsf_queue eq "") {
	$lsf_queue = "scavenger";
}

if ($lsf_time eq "") {
	$lsf_time = "12";
}

if ($lsf_memory eq "") {
	$lsf_memory = 20000;
}

if ($human eq "") {
	$human = "Y";
}
elsif ($human !~ /[YN]/) {
	print("Error: The -human flag must be Y or N or not specified (switched to Y by default).\n"); exit;
}

if ($opt_o !~ /^\//) {
	my $current_path = `pwd`;
	print("Error: Specify the *full* path for the output directory, not the relative.\n\nCorrect:\n\n/home/foobar/projects/my_output_directory\n\nIncorrect: ./my_output_directory\n\nThe full path to this directory is: $current_path\n\n"); exit;
}

if ($opt_f !~ /^\//) {
	print("Error: Specify the *full* path for the input file, not the relative.\n"); exit;
}

if ($revcom_method eq "") {
	$revcom_method = "HMM";
}

if ($revcom_method ne "HMM" && $revcom_method ne "BLAST" && $revcom_method ne "none") {
	print("-R can be HMM, BLAST or none\n"); exit;
}

$distance = "MSA" if ($distance eq "");

if ($distance ne "MSA" && $distance ne "NW") {
	print("-x flag is invalid\n"); exit;
}

if ($distance eq "NW") {
	# In this system equipped with LSF?
	my $init = `bsub -V 2>&1`;

	if ($init eq "") {
		print("This system is not supported. oclust can only run on a cluster environment equipped with the LSF scheduler.\n\n"); exit;
	}
}

# Is this a fasta file?
################################################################################
my $cmd = "grep -P '^>' $opt_f";

if ($Config{osname} =~ /darwin/) {
	$cmd = "grep '^>' $opt_f";
}

my $cmd_out = `$cmd`;

if ($Config{osname} =~ /darwin/ && $debug eq "") {
	print("oclust does not support OS X at the moment."); exit;
}

if ($cmd_out eq "") {
	print("Error: Input file does not appear to be a fasta file\n"); exit;
}

# Create output directory
################################################################################
if (-e $opt_o) {
	print("Error: output directory exists\n"); exit;
}

`mkdir $opt_o 2>/dev/null`;

# Filter on sequence length
################################################################################
my $is = Bio::SeqIO->new(-file => "$opt_f", -format => 'fasta');
my $os = Bio::SeqIO->new(-file => ">$opt_o/targets.fa", -format => 'fasta');

my $count = 0;

while (my $obj = $is->next_seq()) {
	if ($opt_min_length ne "" && $opt_max_length ne "") {
		if ($obj->length > $opt_min_length && $obj->length < $opt_max_length) {
			$os->write_seq($obj);

			$count ++ ;
		}
	}
	else {
		$os->write_seq($obj);
		$count ++ ;		
	}
}

if ($count == 0) {
	print("Size filter is too strict. It filtered all sequencess.\n"); exit;
}

# Pick a random subset?
################################################################################
if ($random_subset ne "") {
	my %sequences;
	my %sequences_random;

	my $is = Bio::SeqIO->new(-file => "$opt_o/targets.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$opt_o/targets.ss.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		$sequences{$obj->display_id} = $obj;
	}

	if (keys(%sequences) < $random_subset) {
		print("Error: There are fewer sequences available after length filtering than the random sample requested\n"); exit;
	}

	while (keys(%sequences_random) < $random_subset) {
		my $r = int(rand(keys(%sequences)));
		my $obj = $sequences{(keys(%sequences))[$r]};

		$sequences_random{$obj->display_id} = $obj;
	}

	foreach my $rand (keys(%sequences_random)) {
		$os->write_seq($sequences_random{$rand});
	}

	print("Picked $random_subset random sequences of " . keys(%sequences) . ".\n");
}
else {
	`cp -v $opt_o/targets.fa $opt_o/targets.ss.fa`;
}

# Reverse complement reads
################################################################################

if ($revcom_method eq "BLAST") {
	`$cwd/bin/formatdb -pF -i $cwd/db/gg_13_5_99.20000.fasta`;
	`rm formatdb.log`;

	print("Orienting sequences in the same direction. Please be patient.\n");
	`$cwd/bin/blastall -a $opt_p -d $cwd/db/gg_13_5_99.20000.fasta -m 8 -F F -p blastn -v 1 -b 1 -e 1e-10 < $opt_o/targets.ss.fa > $opt_o/blast_screen.out`;

	print("BLAST finished. Orienting sequences.\n");

	open(fh, "$opt_o/blast_screen.out") or die("Cannot open $opt_o/blast_screen.out");

	my %list;
	my %targets;

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+)\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)\t(\S+)\t/;
		my ($id, $start, $stop) = ($1,$2,$3);

		if ($list{$id} eq "") {
			$list{$id} = 1;

			if ($start > $stop) {
				$targets{$id} = 1;
			}
		}
	}

	close(fh);

	my $is = Bio::SeqIO->new(-file => "$opt_o/targets.ss.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$opt_o/targets.ss.F.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		my $foo = $obj;

		if ($targets{$obj->display_id} ne "") {
			$foo = $foo->revcom();
		}

		$os->write_seq($foo);
	}
}
elsif ($revcom_method eq "HMM") {
	$ENV{PATH} = $cwd . "/bin/:" . $ENV{PATH};

	my $cmd = "$cwd/bin/vrevcomp/vrevcomp.pl -t $opt_o/ -h $cwd/bin/vrevcomp/HMMs/SSU/bacteria/ -c $opt_o" . "/vrevcomp.csv -o $opt_o/" . "vrevcomp.fa $opt_o/targets.ss.fa";

	print("Orienting...\n");
	`$cmd`;

	`rm -v progress.csv`;

	open(fh, "$opt_o/vrevcomp.csv") or die("Cannot open $opt_o/vrevcomp.csv");

	my %list;
	my %targets;

	<fh>;

	while (my $line = <fh>) {
		chomp($line);

		if ($line =~ /,reverse,/) {
			$line =~ /^(\S+?),/;
			my $id = $1;

			$targets{$id} = 1;
		}
	}

	my $is = Bio::SeqIO->new(-file => "$opt_o/targets.ss.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$opt_o/targets.ss.F.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		my $foo = $obj;

		if ($targets{$obj->display_id} ne "") {
			$foo = $foo->revcom();
		}

		$os->write_seq($foo);
	}

	print("Done orienting\n");
}

# Screen for human contamination
################################################################################
if ($human eq "Y") {
	if ($debug ne "") {
		print("Debug: ./bin/megablast -a $opt_p -d $cwd/db/hg19.fa -m 8 -i $opt_o/targets.ss.F.fa -o $opt_o/blast_targets.hg19.out -v 10 -b 10 -e 1e-5 -F F\n");
	}

	# First check if it's needed to download the human genome
	my $hg_path = $cwd . "/db/hg19.fa";

	if ( ! -e $hg_path) {
		print("Downloading the human genome sequence.\n");

		my $cmd = "wget -P " . $cwd . "db --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr*' 2>/dev/null >/dev/null";
		`$cmd`;

		my @files = <$cwd/db/*.gz>;

		foreach my $file (@files) {
			print("Decompressing: $file\n");
			`gunzip $file`;
		}

		print("Merging.\n");

		my $cmd = "cat $cwd" . "db/chr*.fa > $cwd" . "db/hg19.fa";
		`$cmd`;

		my $cmd = "rm -rf $cwd"."db/chr*.fa";
		`$cmd`;

		print("Formating database.\n");

		my $cmd = "formatdb -pF -i $cwd" . "db/hg19.fa";
		`$cmd`;

		my $cmd = "rm $cwd" . "formatdb.log";
		`$cmd`;
	}

	print("Screening for human contamination. This may take a while.\n");

	`$cwd/bin/megablast -a $opt_p -d $cwd/db/hg19.fa -m 8 -i $opt_o/targets.ss.F.fa -o $opt_o/blast_targets.hg19.out -v 10 -b 10 -e 1e-5 -F F`;

	open(fh, "$opt_o/blast_targets.hg19.out") or die("Cannot open $opt_o/blast_targets.hg19.out");

	my %list;
	my %targets;

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+)\t/;
		my $id = $1;
		$targets{$id} = 1;
	}

	close(fh);

	my $is = Bio::SeqIO->new(-file => "$opt_o/targets.ss.F.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$opt_o/targets.ss.FF.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		if ($targets{$obj->display_id} eq "") {
			$os->write_seq($obj);
		}
	}
}
else {
	print("Skipping human contamination check.\n");
	`cp $opt_o/targets.ss.F.fa $opt_o/targets.ss.FF.fa`;
}

# Chimera detection
################################################################################
if ($chimera eq "Y") {
	print("Running chimera check.\n");

	`$cwd/bin/uchime4.2.40_i86linux32 --minh 1.0 --quiet --db $cwd/db/gold.fa --input $opt_o/targets.ss.FF.fa --uchimeout $opt_o/targets.ss.FF.fa.uchime`;

	open(fh, "$opt_o/targets.ss.FF.fa.uchime") or die("Cannot open $opt_o/targets.ss.FF.fa.uchime");
	my %removal;

	while (my $line = <fh>) {
		chomp($line);

		if ($line =~ /\tY$/) {
			$line =~ /^\S+\t(\S+)\t/;
			my $id = $1;
			$removal{$id} = 1;
		}
	}

	print(keys(%removal) . " sequences were flaged as chimeras, these will be removed.\n");

	my $is = Bio::SeqIO->new(-file => "$opt_o/targets.ss.FF.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$opt_o/targets.ss.FF.no_chimeras.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		if ($removal{$obj->display_id} eq "") {
			$os->write_seq($obj);
		}
	}

	close(fh);
}
else {
	print("Skipping chimera check.\n\n");
	`cp $opt_o/targets.ss.FF.fa $opt_o/targets.ss.FF.no_chimeras.fa`;	
}

# Remove colon from the header (needle truncates whatever is before the colon)
################################################################################
my $cmd = "awk \'{
	if (\$0 ~ /^>/) {
		id=\$0;
		gsub(\":\",\"_\",id);
		print(id)
	}
	else {
		print(\$0)
	}
}' $opt_o/targets.ss.FF.no_chimeras.fa > $opt_o/targets.ss.FF.C.fa";

system($cmd);

# Cleanup
system("rm error.log 2> /dev/null");

my $fs = -s "$opt_o/targets.ss.FF.C.fa";

if ($fs == 0) {
	print("Error: No remaining sequences. Cannot continue.\n\n"); exit;
}
else {
	print("Output files written to: $opt_o\n\n");
}

# Write and launch job files
################################################################################
my $total_sequence_count = `grep '>' $opt_o/targets.ss.FF.C.fa | wc -l`;
chomp($total_sequence_count);
my $total_job_file_count = ceil($total_sequence_count/$lsf_nb_jobs);
my $cmd = $cwd . "utils/split_fasta_file.pl $opt_o/targets.ss.FF.C.fa $total_job_file_count";

`$cmd`;

`mkdir $opt_o/jobs`;
`mkdir $opt_o/logs`;
`mkdir $opt_o/alignments`;
`mkdir $opt_o/tmp`;

my $db = $opt_o . "/targets.ss.FF.C.fa";

my @files = <$opt_o/*.fasta>;
my $count = 0;

foreach my $file (@files) {
	if (-s $file > 0) {
		$count ++ ;

		my $job_script = "#BSUB -L /bin/bash
#BSUB -n 1
#BSUB -J oclust_".$count."
#BSUB -oo ../logs/$count.log
#BSUB -eo ../logs/$count.err
#BSUB -q $lsf_queue
#BSUB -R rusage[mem=$lsf_memory]
#BSUB -W $lsf_time" . ":00\n";

		if ($lsf_account ne "") {
			$job_script .= "#BSUB -P $lsf_account\n";
		}

		$job_script .= "cd $opt_o
$cwd";

		$job_script .= "utils/needle_runner.pl $file $db $opt_o/alignments/$count" . ".needle.out $opt_o/tmp $count\n\n";

		open(fh, ">".$opt_o."/jobs/$count".".job");
		print(fh $job_script);
		close(fh);
	}
}

print("Creating $count job file(s).\n");

# Submit jobs and record job identifiers so that we can track them
my @files = <$opt_o/jobs/*>;
my @ids;

foreach my $job (@files) {
	my $cmd = "cd $opt_o/jobs";

	$job =~ /^.+\/(\S+)$/;
	my $fn = $1;

	my $q = "$cmd; bsub < $fn";
	my $out = `$q`;

	$out =~ /^\S+ <(\S+)>/;
	my $job_id = $1;

	push(@ids, $job_id);

	print("Submitted $job\n");
}

# Write running jobs
open(fh, ">" . $opt_o . "/submitted_jobs.txt");

foreach my $id (@ids) {
	print(fh "$id\n");
}

close(fh);