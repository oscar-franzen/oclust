#!/usr/bin/perl
# -----------------------------------------------------------------------------------------
# oclust, a pipeline for clustering sequences into Operational Taxonomic Units.           #
# Oscar Franzen, <p.oscar.franzen@gmail.com>.                                             #
# -----------------------------------------------------------------------------------------

use strict 'vars';
use Getopt::Long;
use Config;
use POSIX;

# get the absolute path to the current working directory
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

use lib workaround(); # path to bundled perl modules

use Bio::SeqIO;
use Bio::AlignIO;

if ($Config{osname} !~ /linux/) {
	die("oclust must be run on a Linux system.");
}

check_file_exists($cwd."bin/blastall");
check_file_exists($cwd."bin/megablast");
check_file_exists($cwd."bin/formatdb");
check_file_exists($cwd."bin/uchime4.2.40_i86linux32");
check_file_exists($cwd."bin/hmmscan");
check_file_exists($cwd."bin/needle");
check_file_exists($cwd."bin/needle.acd");
check_file_exists($cwd."bin/codes.english");
check_file_exists($cwd."bin/knowntypes.standard");
check_file_exists($cwd."bin/EDNAFULL");
check_file_exists($cwd."bin/cmbuild");
check_file_exists($cwd."bin/cmalign");

my $info;

if (-e "/proc/cpuinfo") {
	my $o=`cat /proc/cpuinfo | grep processor | tail -n1 | awk '{print \$3}'`;
	chomp($o);
	$o++;

	$info = "oclust is running on $Config{osname} ($Config{archname}), $o CPUs are available on this system.";
}
else {
	$info = "oclust is running on $Config{osname} ($Config{archname})";
}

die("$info\nFeedback: <p.oscar.franzen\@gmail.com>, Mount Sinai, New York, U.S.A.\n\nCommand line arguments:

	-x <method> -f <input file> -o <output directory> -p <number of CPUs>

	General settings:
	-x PW or MSA               Can be PW for pairwise alignments (based on Needleman-Wunsch)
	                           or MSA for multiple sequence alignment (based on Infernal). [MSA]
	-t local or cluster        If -x is PW, should it be parallelized by running it locally
	                           on multiple cores or by submitting jobs to a cluster. [local]
	-a complete, average or    The desired clustering algorithm. [complete]
	single    
	-f [string]                Input fasta file.
	-o [string]                Name of output directory (must not exist) and use full path.
	-R HMM, BLAST, or none     Method to use for reverse complementing sequences. [HMM]
	-p [integer]               Number of processor cores to use for BLAST. [4]
	-minl [integer]            Minimum sequence length. [optional]
	-maxl [integer]            Maximum sequence length. [optional]
	-rand [integer]            Randomly sample a specified number of sequences. [optional]
	-human Y or N              If 'Y'es, then execute BLAST-based contamination
	                           screen towards the human genome. [Y]
	-chimera Y or N            Run chimera check. Can be Y or N. [Y]

	LSF settings (only valid for -x PW when -t cluster):
	-lsf_queue [string]       Name of the LSF queue to use. [scavenger]
	-lsf_account [string]     Name of the account to use. [optional]
	-lsf_time [integer]       Runtime hours per job specified as number of hours. [1]
	-lsf_memory [integer]     Requested amount of RAM in MB. [3000]
	-lsf_nb_jobs [integer]    Number of jobs. [20]

	Usage example: -x MSA -f /home/foobar/long_reads.fasta -o /home/foobar/foo -p 4 -minl 700 -maxl 800\n\n") if (@ARGV == 0);

my $setting_input_file;
my $setting_output_dir;
my $setting_cpus;
my $setting_min_seq_length;
my $setting_max_seq_length;
my $setting_debug;
my $setting_random_subset;
my $setting_human_contamination_check;
my $setting_chimera_check;
my $setting_revcom_method;
my $setting_distance_method;
my $setting_hclust_algorithm;
my $setting_parallel_type;
my $setting_lsf_nb_jobs;
my $setting_lsf_queue;
my $setting_lsf_time;
my $setting_lsf_memory;
my $setting_lsf_account;

GetOptions ("file=s" => \$setting_input_file,
	         "out=s" => \$setting_output_dir,
	         "proc=i" => \$setting_cpus,
	         "minl=i" => \$setting_min_seq_length,
	         "maxl=i" => \$setting_max_seq_length,
	         "debug=i" => \$setting_debug,
	         "rand=i" => \$setting_random_subset,
	         "human=s" => \$setting_human_contamination_check,
	         "chimera=s" => \$setting_chimera_check,
	         "lsf_nb_jobs=i" => \$setting_lsf_nb_jobs,
	         "lsf_queue=s" => \$setting_lsf_queue,
	         "lsf_account=s" => \$setting_lsf_account,
	         "lsf_time=i" => \$setting_lsf_time,
	         "lsf_memory=i" => \$setting_lsf_memory,
	         "lsf_account=s" => \$setting_lsf_account,
	         "R=s" => \$setting_revcom_method,
	         "x=s" => \$setting_distance_method,
	         "algorithm=s" => \$setting_hclust_algorithm,
	         "type=s" => \$setting_parallel_type) or die("Error in command line arguments\n");

checkSettings();

# Is this a fasta file?
################################################################################
my $cmd = "grep -P '^>' $setting_input_file";
my $cmd_out = `$cmd`;

if ($cmd_out eq "") {
	die("Error: Input file does not appear to be a fasta file\n");
}

# Create output directory
################################################################################
if (-e $setting_output_dir) {
	#die("Error: output directory exists\n");
}

`mkdir $setting_output_dir 2>/dev/null`;

# Filter on sequence length
################################################################################
my $is = Bio::SeqIO->new(-file => "$setting_input_file", -format => 'fasta');
my $os = Bio::SeqIO->new(-file => ">$setting_output_dir/targets.fa", -format => 'fasta');

my $count = 0;

while (my $obj = $is->next_seq()) {
	if ($setting_min_seq_length ne "" && $setting_max_seq_length ne "") {
		if ($obj->length > $setting_min_seq_length && $obj->length < $setting_max_seq_length) {
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
if ($setting_random_subset ne "") {
	my %sequences;
	my %sequences_random;

	my $is = Bio::SeqIO->new(-file => "$setting_output_dir/targets.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$setting_output_dir/targets.ss.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		$sequences{$obj->display_id} = $obj;
	}

	if (keys(%sequences) < $setting_random_subset) {
		print("Error: There are fewer sequences available after length filtering than the random sample requested\n"); exit;
	}

	while (keys(%sequences_random) < $setting_random_subset) {
		my $r = int(rand(keys(%sequences)));
		my $obj = $sequences{(keys(%sequences))[$r]};

		$sequences_random{$obj->display_id} = $obj;
	}

	foreach my $rand (keys(%sequences_random)) {
		$os->write_seq($sequences_random{$rand});
	}

	print("Picked $setting_random_subset random sequences of " . keys(%sequences) . ".\n");
}
else {
	`cp -v $setting_output_dir/targets.fa $setting_output_dir/targets.ss.fa`;
}

# Reverse complement reads
################################################################################

if ($setting_revcom_method eq "BLAST") {
	`$cwd/bin/formatdb -l /dev/null -pF -i $cwd/db/gg_13_5_99.20000.fasta`;

	print("Orienting sequences in the same direction. Please be patient.\n");
	`$cwd/bin/blastall -a $setting_cpus -d $cwd/db/gg_13_5_99.20000.fasta -m 8 -F F -p blastn -v 1 -b 1 -e 1e-10 < $setting_output_dir/targets.ss.fa > $setting_output_dir/blast_screen.out`;

	print("BLAST finished. Orienting sequences.\n");

	open(fh, "$setting_output_dir/blast_screen.out") or die("Cannot open $setting_output_dir/blast_screen.out");

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

	my $is = Bio::SeqIO->new(-file => "$setting_output_dir/targets.ss.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$setting_output_dir/targets.ss.F.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		my $foo = $obj;

		if ($targets{$obj->display_id} ne "") {
			$foo = $foo->revcom();
		}

		$os->write_seq($foo);
	}
}
elsif ($setting_revcom_method eq "HMM") {
	$ENV{PATH} = $cwd . "/bin/:" . $ENV{PATH};

	my $cmd = "$cwd/bin/vrevcomp/vrevcomp.pl -t $setting_output_dir/ -h $cwd/bin/vrevcomp/HMMs/SSU/bacteria/ -c $setting_output_dir" . "/vrevcomp.csv -o $setting_output_dir/" . "vrevcomp.fa $setting_output_dir/targets.ss.fa";

	print("Orienting...\n");
	`$cmd`;

	`rm -v progress.csv`;

	open(fh, "$setting_output_dir/vrevcomp.csv") or die("Cannot open $setting_output_dir/vrevcomp.csv");

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

	my $is = Bio::SeqIO->new(-file => "$setting_output_dir/targets.ss.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$setting_output_dir/targets.ss.F.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		my $foo = $obj;

		if ($targets{$obj->display_id} ne "") {
			$foo = $foo->revcom();
		}

		$os->write_seq($foo);
	}

	print("Done orienting\n");
}
else {
	my $cmd = "cp -v $setting_output_dir/targets.ss.fa $setting_output_dir/targets.ss.F.fa";
	`$cmd`;
}

# Screen for human contamination
################################################################################
if ($setting_human_contamination_check eq "Y") {
	if ($setting_debug ne "") {
		print("Debug: ./bin/megablast -a $setting_cpus -d $cwd/db/hg19.fa -m 8 -i $setting_output_dir/targets.ss.F.fa -o $setting_output_dir/blast_targets.hg19.out -v 10 -b 10 -e 1e-5 -F F\n");
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

		my $cmd = "formatdb -l /dev/null -pF -i $cwd" . "db/hg19.fa";
		`$cmd`;
	}

	print("Screening for human contamination. This may take a while.\n");

	`$cwd/bin/megablast -a $setting_cpus -d $cwd/db/hg19.fa -m 8 -i $setting_output_dir/targets.ss.F.fa -o $setting_output_dir/blast_targets.hg19.out -v 10 -b 10 -e 1e-5 -F F`;

	open(fh, "$setting_output_dir/blast_targets.hg19.out") or die("Cannot open $setting_output_dir/blast_targets.hg19.out");

	my %list;
	my %targets;

	while (my $line = <fh>) {
		chomp($line);

		$line =~ /^(\S+)\t/;
		my $id = $1;
		$targets{$id} = 1;
	}

	close(fh);

	my $is = Bio::SeqIO->new(-file => "$setting_output_dir/targets.ss.F.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$setting_output_dir/targets.ss.FF.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		if ($targets{$obj->display_id} eq "") {
			$os->write_seq($obj);
		}
	}
}
else {
	print("Skipping human contamination check.\n");
	`cp $setting_output_dir/targets.ss.F.fa $setting_output_dir/targets.ss.FF.fa`;
}

# Chimera detection
################################################################################
if ($setting_chimera_check eq "Y") {
	print("Running chimera check.\n");

	`$cwd/bin/uchime4.2.40_i86linux32 --minh 1.0 --quiet --db $cwd/db/gold.fa --input $setting_output_dir/targets.ss.FF.fa --uchimeout $setting_output_dir/targets.ss.FF.fa.uchime`;

	open(fh, "$setting_output_dir/targets.ss.FF.fa.uchime") or die("Cannot open $setting_output_dir/targets.ss.FF.fa.uchime");
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

	my $is = Bio::SeqIO->new(-file => "$setting_output_dir/targets.ss.FF.fa", -format => 'fasta');
	my $os = Bio::SeqIO->new(-file => ">$setting_output_dir/targets.ss.FF.no_chimeras.fa", -format => 'fasta');

	while (my $obj = $is->next_seq()) {
		if ($removal{$obj->display_id} eq "") {
			$os->write_seq($obj);
		}
	}

	close(fh);
}
else {
	print("Skipping chimera check.\n\n");
	`cp $setting_output_dir/targets.ss.FF.fa $setting_output_dir/targets.ss.FF.no_chimeras.fa`;	
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
}' $setting_output_dir/targets.ss.FF.no_chimeras.fa > $setting_output_dir/targets.ss.FF.C.fa";

system($cmd);

# Cleanup
system("rm error.log 2> /dev/null");

my $fs = -s "$setting_output_dir/targets.ss.FF.C.fa";

if ($fs == 0) {
	print("Error: No remaining sequences. Cannot continue.\n\n"); exit;
}
else {
	print("Output files will be written to: $setting_output_dir\n\n");
}

# Pairwise alignments
################################################################################
if ($setting_distance_method eq "PW" && $setting_parallel_type eq "local") {
	my @seqs;
	my $is = Bio::SeqIO->new(-file => $setting_output_dir."/targets.ss.FF.C.fa", -format => "fasta");

	while (my $obj = $is->next_seq()) {
		push(@seqs, $obj);
	}

	$ENV{EMBOSS_ACDROOT} = $cwd . "/bin/";

	my %done;

	for (my $i=0; $i<@seqs; $i++) {
		my $seq1 = @seqs[$i];

		# Write this sequence
		open(fh_out, ">$setting_output_dir" . "/q_" . $i . ".fa");
		print(fh_out ">". $seq1->display_id . "\n");
		print(fh_out $seq1->seq . "\n");
		close(fh_out);

		open(fh_out, ">$setting_output_dir" . "/t_" . $i . ".fa");

		# Write the db
		for (my $c=0; $c<@seqs; $c++) {
			my $seq2 = @seqs[$c];

			if ($i != $c) {
				if ($done{$seq1->display_id . "-" . $seq2->display_id} eq "") {
					print(fh_out ">". $seq2->display_id . "\n");
					print(fh_out $seq2->seq . "\n");

					# Mark comparison as done
					$done{$seq1->display_id . "-" . $seq2->display_id} = 1;
					$done{$seq2->display_id . "-" . $seq1->display_id} = 1;
				}
			}
		}

		close(fh_out);
	}

	print("Running $setting_cpus threads.\n");

	# Bundle them together into jobs
	my $partition = int( @seqs / $setting_cpus );
	my $file_suffix = 1;
	my $count = 0;
	my $out_suffix = 0;

	open(my $fh_out, ">$setting_output_dir" . "/partition_" . $file_suffix . ".job");

	my @files = <$setting_output_dir/t_*.fa>;

	fisher_yates_shuffle( \@files );

	foreach my $file (@files) {
		my $size = -s $file;
		$count ++ ;

		if ($size > 0) {
			my $query = $file;
			my $db = $file;

			$query =~ s/\/t_/\/q_/;

			if ($count > $partition) {
				$file_suffix ++ ;
				$count = 0;

				close($fh_out);
				open($fh_out, ">$setting_output_dir" . "/partition_" . $file_suffix . ".job");
			}

			$out_suffix ++ ;

			my $cmd = $cwd . "bin/needle -asequence $query -bsequence $db -datafile " . $cwd . "bin/EDNAFULL -auto -stdout -aformat3 fasta > $setting_output_dir/" . "needle_" . $out_suffix . ".aln";

			print($fh_out "$cmd\n");
		}
	}

	close($fh_out);

	my @files = <$setting_output_dir/*.job>;

	foreach my $job (@files) {
		`chmod +x $job`;
	}

	my @files = <$setting_output_dir/partition*.job>;

	open(fh_out, ">$setting_output_dir/" . "run_pw");
	print(fh_out "cd $setting_output_dir\n");

	foreach my $file (@files) {
		print(fh_out "$file &\nsleep 0.2s\n");
	}

	print(fh_out "wait");

	close(fh_out);

	print("Running alignments. This may take a while.\n");

	my $p = $setting_output_dir . "/run_pw";
	#system("chmod +x $p; $p");

	finish();
}
elsif ($setting_distance_method eq "PW" && $setting_parallel_type eq "cluster") {
	my @seqs;
	my $is = Bio::SeqIO->new(-file => $setting_output_dir."/targets.ss.FF.C.fa", -format => "fasta");

	while (my $obj = $is->next_seq()) {
		push(@seqs, $obj);
	}

	# Number of alignments per CPU
	# Divided by two, because performing the alignment A:B is the same as B:A.
	my $partition = int( ((@seqs * @seqs)/2) / $setting_lsf_nb_jobs);

	print("Submitting $partition pairwise alignments per job...\n");

	my $count = 0;
	my $file_suffix = 1;

	open(my $fh_out, ">$setting_output_dir" . "/partition_" . $file_suffix . ".fa");

	my %done;

	for (my $i=0; $i<@seqs; $i++) {
		my $seq1 = @seqs[$i];

		for (my $c=0; $c<@seqs; $c++) {
			my $seq2 = @seqs[$c];

			if ($i != $c) {
				if ($done{$c.":".$i} eq "") {
					print($fh_out ">".$seq1->display_id . "\n");
					print($fh_out $seq1->seq . "\n");

					print($fh_out ">".$seq2->display_id . "\n");
					print($fh_out $seq2->seq . "\n");

					$count ++ ;

					if ($count > $partition) {
						$file_suffix ++ ;
						$count = 0;
						open($fh_out, ">$setting_output_dir" . "/partition_" . $file_suffix . ".fa");
					}

					$done{$c.":".$i} = 1;
					$done{$i.":".$c} = 1;
				}
			}
		}
	}

	`mkdir $setting_output_dir/jobs`;
	`mkdir $setting_output_dir/logs`;

	for (my $i=1; $i<=$file_suffix; $i++) {
		my $job_script = "#BSUB -L /bin/bash
#BSUB -n 1
#BSUB -J oclust_".$i."
#BSUB -oo ../logs/$i.log
#BSUB -eo ../logs/$i.err
#BSUB -q $setting_lsf_queue
#BSUB -R rusage[mem=$setting_lsf_memory]
#BSUB -W $setting_lsf_time" . ":00\n";

		if ($setting_lsf_account ne "") {
			$job_script .= "#BSUB -P $setting_lsf_account\n";
		}

		$job_script .= "cd $setting_output_dir\n";
		$job_script .= $cwd . "bin/needleman_wunsch --printfasta --file " . $setting_output_dir . "/partition_" . $i . ".fa > " . $setting_output_dir . "/partition_" . $i . ".fa.fas\n";

		open(fh, ">".$setting_output_dir."/jobs/$i".".job");
		print(fh $job_script);
		close(fh);
	}

	my @files = <$setting_output_dir/jobs/*>;
	my @ids;

	foreach my $job (@files) {
		my $cmd = "cd $setting_output_dir/jobs";

		$job =~ /^.+\/(\S+)$/;
		my $fn = $1;

		my $q = "$cmd; bsub < $fn";
		my $out = `$q`;

		$out =~ /^\S+ <(\S+)>/;
		my $job_id = $1;

		push(@ids, $job_id);

		print("Submitted $job\n");
	}

	print(@files . " jobs have been submitted to the cluster. Now waiting for jobs to finish.\n");

	while (1) {
		# Check if the number of completed log files correspond to the number of submitted jobs
		my @logs = <$setting_output_dir/logs/*.log>;
		my $found_completed_logs = 0;

		foreach my $log (@logs) {
			open(fh, $log);

			while (my $line = <fh>) {
				chomp($line);
				
				if ($line =~ /Successfully completed./) {
					$found_completed_logs ++ ;
					last;
				}
			}

			close(fh);
		}

		if ($found_completed_logs == $setting_lsf_nb_jobs) {
			print("All submitted jobs have completed.\n");
			last;
		}

		sleep(60);
	}

	finish();
}
else {
	# Build the covariance model
	my $dir_path = "$cwd" . "/bin/RDPinfernalTraindata";

	if (! -d $dir_path) {
		my $cmd = "unzip -d $cwd" . "bin $cwd" . "bin/RDPinfernalTraindata.zip";
		`$cmd`;
	}

	my $f = $cwd . "bin/RDPinfernalTraindata/bacteria16S_508_mod5.stk.cm";

	if (! -e $f) {
		print("Running cmbuild.\n");
		my $cmd = $cwd . "bin/cmbuild --ere 1.4 $cwd"."bin/RDPinfernalTraindata/bacteria16S_508_mod5.stk.cm $cwd"."bin/RDPinfernalTraindata/bacteria16S_508_mod5.stk";

		`$cmd`;
	}

	# Infernal-based
	print("Running cmalign.\n");
	my $f = $setting_output_dir . "/targets.ss.FF.C.fa";
	my $cmd = $cwd."bin/cmalign --cpu $setting_cpus -o $setting_output_dir" . "/infernal.sto $cwd" . "bin/RDPinfernalTraindata/bacteria16S_508_mod5.stk.cm $f 2>/dev/null >/dev/null";

	`$cmd`;

	# Convert to fasta and remove sequences without any homologous positions
	my $in = Bio::AlignIO->new(-file => $setting_output_dir."/infernal.sto", '-format' => 'stockholm');
	my @alignments;

	while ( my $aln = $in->next_aln() ) {
		foreach my $seq ($aln->each_seq()) {
			push(@alignments, $seq);
		}
	}

	print("Checking alignment for consistency.\n");

	my %res;

	for (my $i=0; $i<@alignments; $i++) {
		my $seq1 = @alignments[$i];

		for (my $k=0; $k<@alignments; $k++) {
			if ($k != $i) {
				my $seq2 = @alignments[$k];
				my $str1 = $seq1->seq;
				my $str2 = $seq2->seq;

				my $count = 0;

				for (my $q=0; $q<length($str1); $q++) {
					my $char1 = substr($str1, $q, 1);
					my $char2 = substr($str2, $q, 1);

					if ($char1 =~ /[ATGC]/i && $char2 =~ /[ATGC]/i) {
						$count ++ ;
					}
				}

				if ($count == 0) {
					$res{$i} ++ ;
				}
			}
		}
	}

	my $n = @alignments;
	$n -- ;
	my %removal;

	foreach my $index (keys(%res)) {
		if ($res{$index} == $n) {
			#print("$index $res{$index} " . @alignments[$index]->display_id . "\n");
			$removal{@alignments[$index]->display_id} = 1;
		}
	}

	my $os = Bio::SeqIO->new(-file => ">" . $setting_output_dir . "/infernal.F.fasta", -format => "fasta");

	print(keys(%removal) . " sequences removed\n");

	foreach my $item (@alignments) {
		if ($removal{$item->display_id} eq "") {
			$os->write_seq($item);
		}
	}

	my $cmd = "cat $cwd/bin/R/bin/R";
	my $o = `$cmd`;

	$o =~ s/\/home\/foobar\/oclust\//$cwd/g;
	$o =~ s/R_installed/R/g;

	open(fh, ">" . $cwd."/bin/R/bin/R.fixed");
	print(fh "$o\n");
	close(fh);

	my $cmd = "chmod +x $cwd/bin/R/bin/R.fixed";
	`$cmd`;

	my $cmd = "$cwd/bin/R/bin/R.fixed --no-save --no-restore --args $cwd $setting_output_dir"."/infernal.F.fasta $setting_output_dir $setting_hclust_algorithm MSA < $cwd/utils/hclust_fr_aln.R";
	`$cmd`;

	print("*** oclust running in MSA-mode has finished. *** \nResults are in:\n$setting_output_dir\n");
}

sub finish {
	print("Alignment has finished. Writing distance matrix.\n");

	my @files = <$setting_output_dir/*.aln>;
	my %distances;
	my %all_ids;

	foreach my $file (@files) {
		my $is = Bio::SeqIO->new(-file => $file, -format => "fasta");

		while (my $obj1 = $is->next_seq()) {
			my $obj2 = $is->next_seq();

			my $id1 = $obj1->display_id;
			my $id2 = $obj2->display_id;

			$all_ids{$id1} = 1;
			$all_ids{$id2} = 1;

			my $seq1 = $obj1->seq;
			my $seq2 = $obj2->seq;

			# Remove terminal gaps
			$seq1 = uc($seq1);
			$seq2 = uc($seq2);

			if ($seq1 =~ /^\-/) {
				$seq1 =~ /^(-+)(\S+)/;
				my $gaps_str = $1;
				$seq1 = $2;
				$seq2 = substr($seq2, length($gaps_str), length($seq2));
			}

			if ($seq2 =~ /^\-/) {
				$seq2 =~ /^(-+)(\S+)/;
				my $gaps_str = $1;
				$seq2 = $2;
				$seq1 = substr($seq1, length($gaps_str), length($seq1));
			}

			if ($seq1 =~ /\-+$/) {
				$seq1 =~ /(\-+)$/;
				my $gaps_str = $1;

				$seq1 = substr($seq1, 0, length($seq1) - length($gaps_str));
				$seq2 = substr($seq2, 0, length($seq2) - length($gaps_str));
			}

			if ($seq2 =~ /\-+$/) {
				$seq2 =~ /(\-+)$/;
				my $gaps_str = $1;

				$seq2 = substr($seq2, 0, length($seq2) - length($gaps_str));
				$seq1 = substr($seq1, 0, length($seq1) - length($gaps_str));
			}

			# Delete single base indels
			my @positions_to_delete;

			while ($seq2 =~ m/[ATGC](-)[ATGC]/g) {
				my $p = pos($seq2);
				push(@positions_to_delete, $p - 2);
			}

			for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($seq1, $item, 1) = "";
				substr($seq2, $item, 1) = "";
			}

			my @positions_to_delete;

			while ($seq1 =~ m/[ATGC](-)[ATGC]/g) {
				my $p = pos($seq1);
				push(@positions_to_delete, $p - 2);
			}

			for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($seq1, $item, 1) = "";
				substr($seq2, $item, 1) = "";
			}

			# Delete indels larger than 4
			my @positions_to_delete;

			while ($seq1 =~ m/[ATGC](-{5,})[ATGC]/g) {
				my $del_start = pos($seq1) - length($1);
				push(@positions_to_delete, { del_start => $del_start, del_len => length($1) });
			}

			 for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($seq1, $item->{del_start} - 1, $item->{del_len}) = "";
				substr($seq2, $item->{del_start} - 1, $item->{del_len}) = "";
			}

			my @positions_to_delete;

			while ($seq2 =~ m/[ATGC](-{5,})[ATGC]/g) {
				my $del_start = pos($seq2) - length($1);
				push(@positions_to_delete, { del_start => $del_start, del_len => length($1) });
			}

			 for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($seq1, $item->{del_start} - 1, $item->{del_len}) = "";
				substr($seq2, $item->{del_start} - 1, $item->{del_len}) = "";
			}

			# Calc identity
			my @t1 = split(//, $seq1);
			my @t2 = split(//, $seq2);

			my $differences = 0;
			my $total = 0;

			for (my $i=0; $i<@t1; $i++) {
				if (@t1[$i] ne "N" && @t2[$i] ne "N") {

					# There are some IUPAC-encoded nucleotides in the miseq data - IGNORE!
					if (@t1[$i] =~ /[ATGC\-]/ && @t2[$i] =~ /[ATGC\-]/) {
						if (@t1[$i] ne @t2[$i]) {
							$differences ++ ;
						}

						$total ++ ;
					}
				}
			}

			my $fraction_differences = $differences / $total;
			$fraction_differences = sprintf("%.5f", $fraction_differences);

			$distances{$id1}{$id2} = $fraction_differences;
		}
	}

	open(fh_dist, ">$setting_output_dir" . "/dist.mat");

	my $first_row;

	foreach my $id (keys(%all_ids)) {
		$first_row .= "$id ";
	}

	chop($first_row);
	print(fh_dist "$first_row\n");

	foreach my $id1 (keys(%all_ids)) {
		my $l = "$id1 ";

		foreach my $id2 (keys(%all_ids)) {
			if ($id1 eq $id2) {
				$l .= "0 ";
			}
			elsif ($distances{$id1}{$id2} ne "") {
				$l .= $distances{$id1}{$id2} . " ";
			}
			elsif ($distances{$id2}{$id1} ne "") {
				$l .= $distances{$id2}{$id1} . " ";
			}
		}

		chop($l);

		print(fh_dist "$l\n");
	}

	close(fh_dist);

	my $cmd = "cat $cwd/bin/R/bin/R";
	my $o = `$cmd`;

	$o =~ s/\/home\/foobar\/oclust\//$cwd/g;
	$o =~ s/R_installed/R/g;

	open(fh, ">" . $cwd."/bin/R/bin/R.fixed");
	print(fh "$o\n");
	close(fh);

	my $cmd = "chmod +x $cwd/bin/R/bin/R.fixed";
	`$cmd`;

	my $cmd = "$cwd/bin/R/bin/R.fixed --no-save --no-restore --args $cwd $setting_output_dir"."/dist.mat $setting_output_dir $setting_hclust_algorithm PW < $cwd/utils/hclust_fr_aln.R";
	`$cmd`;

	print("*** oclust running in PW-mode has finished. ***\n\n Results are in:\n$setting_output_dir\n");
}

sub check_file_exists {
	my $file = shift;

	if (! -e $file) {
		print("Error: $file could not be found. Terminating.\n"); exit;
	}
}

sub checkSettings() {
	$setting_parallel_type = "local" if ($setting_parallel_type eq "");

	die("-t must be local or cluster\n") if ($setting_parallel_type ne "local" && $setting_parallel_type ne "cluster");

	$setting_hclust_algorithm = "complete" if ($setting_hclust_algorithm eq "");

	die("-a can be one of: complete, single, or average\n") if ($setting_hclust_algorithm ne "complete" && $setting_hclust_algorithm ne "single" && $setting_hclust_algorithm ne "average");

	die("Error: No fasta input file specified\n") if ($setting_input_file eq "");
	
	die("Error: No output directory specified\n") if ($setting_output_dir eq "");
	
	$setting_cpus = 4 if ($setting_cpus eq "");

	$setting_chimera_check = "Y" if ($setting_chimera_check eq "");

	die("-chimera can take options Y or N.\n") if ($setting_chimera_check !~ /^[YN]$/);

	$setting_lsf_nb_jobs = 20 if ($setting_lsf_nb_jobs eq "");

	$setting_lsf_queue = "scavenger" if ($setting_lsf_queue eq "");

	$setting_lsf_time = "1" if ($setting_lsf_time eq "");

	$setting_lsf_memory = 3000 if ($setting_lsf_memory eq "");

	$setting_human_contamination_check = "Y" if ($setting_human_contamination_check eq "");

	die("Error: The -human flag must be Y or N or not specified (switched to Y by default).\n") if ($setting_human_contamination_check !~ /[YN]/);

	if ($setting_output_dir !~ /^\//) {
		my $current_path = `pwd`;
		die("Error: Specify the *full* path for the output directory, not the relative.\n\nCorrect:\n\n/home/foobar/projects/my_output_directory\n\nIncorrect: ./my_output_directory\n\nThe full path to this directory is: $current_path\n\n");
	}

	die("Error: Specify the *full* path for the input file, not the relative.\n") if ($setting_input_file !~ /^\//);

	$setting_revcom_method = "HMM" if ($setting_revcom_method eq "");

	die("-R can be HMM, BLAST or none\n") if ($setting_revcom_method ne "HMM" && $setting_revcom_method ne "BLAST" && $setting_revcom_method ne "none");

	$setting_distance_method = "MSA" if ($setting_distance_method eq "");

	print("-x must be PW or MSA.\n") if ($setting_distance_method ne "MSA" && $setting_distance_method ne "PW");

	if ($setting_distance_method eq "PW" && $setting_parallel_type eq "cluster") {
		# In this system equipped with LSF?
		my $init = `bsub -V 2>&1`;

		if ($init eq "") {
			die("-x PW -t cluster can only run on a cluster environment equipped with the LSF scheduler. Try -x MSA instead or -x PW -t local.\n\n");
		}
	}
}

sub fisher_yates_shuffle {
    my $deck = shift;
    my $i = @$deck;
    while ($i--) {
        my $j = int rand ($i+1);
        @$deck[$i,$j] = @$deck[$j,$i];
    }
}
