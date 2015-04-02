#!/usr/bin/perl
# oclust, a pipeline for clustering long PacBio CCS reads into operational taxonomic units.
# Oscar Franzen, <p.oscar.franzen@gmail.com>, Mt. Sinai, NY, USA.
# -----------------------------------------------------------------------------------------

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
use Bio::AlignIO;
use Parallel::ForkManager;

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
	-lsf_queue [string]        Name of the LSF queue to use [scavenger].
	-lsf_account [string]      Name of the account to use.
	-lsf_time [integer]        Runtime hours per job specified as number of hours. [12]
	-lsf_memory [integer]      Requested amount of RAM in MB. [20000]
	-lsf_nb_jobs [integer]     Number of jobs. [20]

	Usage example: -x MSA -f /home/foobar/long_reads.fasta -o /home/foobar/foo -p 4 -minl 700 -maxl 800\n\n") if (@ARGV == 0);

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
my $hclust_algorithm;
my $parallel_type;

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
	         "x=s" => \$distance,
	         "algorithm=s" => \$hclust_algorithm,
	         "type=s" => \$parallel_type) or die("Error in command line arguments\n");

if ($parallel_type eq "") {
	$parallel_type = "local";
}

if ($parallel_type ne "local" && $parallel_type ne "cluster") {
	print("-t must be local or cluster\n"); exit;
}

if ($hclust_algorithm eq "") {
	$hclust_algorithm = "complete";
}

if ($hclust_algorithm ne "complete" && $hclust_algorithm ne "single" && $hclust_algorithm ne "average") {
	print("-a can be one of: complete, single, or average\n"); exit;
}

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

if ($distance ne "MSA" && $distance ne "PW") {
	print("-x must be PW or MSA.\n"); exit;
}

if ($distance eq "PW") {
	# In this system equipped with LSF?
	my $init = `bsub -V 2>&1`;

	if ($init eq "") {
		print("-x PW can only run on a cluster environment equipped with the LSF scheduler. Try -x MSA instead.\n\n"); exit;
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
	#print("Error: output directory exists\n"); exit;
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
else {
	my $cmd = "cp -v $opt_o/targets.ss.fa $opt_o/targets.ss.F.fa";
	`$cmd`;
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
	print("!\n");

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
	print("Output files will be written to: $opt_o\n\n");
}

# Write and launch job files
################################################################################
if ($distance eq "PW") {
	my @seqs;
	my $is = Bio::SeqIO->new(-file => $opt_o."/targets.ss.FF.C.fa", -format => "fasta");

	while (my $obj = $is->next_seq()) {
		push(@seqs, $obj);
	}

	# Number of alignments per CPU
	my $partition = int((@seqs * @seqs) / $opt_p);

	print("Performing $partition pairwise alignments per CPU...\n");

	my $count = 0;
	my $file_suffix = 1;

	open(my $fh_out, ">$opt_o" . "/partition_" . $file_suffix . ".fa");

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
						open($fh_out, ">$opt_o" . "/partition_" . $file_suffix . ".fa");
					}

					$done{$c.":".$i} = 1;
					$done{$i.":".$c} = 1;
				}
			}
		}
	}

	my $pm = Parallel::ForkManager->new($opt_p);

	for (my $i=1; $i<=$file_suffix; $i++) {
		my $pid = $pm->start and next;

		# Inside the child process
		my $cmd = $cwd . "/bin/needleman_wunsch --printfasta --file " . $opt_o . "/partition_" . $i . ".fa > " . $opt_o . "/partition_" . $i . ".fa.fas";

		#system($cmd);

		$pm->finish; # Terminates the child process
	}

	print("Running alignments. This may take a while.\n");

	$pm->wait_all_children();

	print("Alignments finished. Writing distance matrix.\n");

	my @files = <$opt_o/*.fas>;
	my %distances;
	my %all_ids;

	foreach my $file (@files) {
		open(fh, "$file\n");

		while (my $id1 = <fh>) {
			my $str1 = <fh>;

			my $id2 = <fh>;
			my $str2 = <fh>;

			chomp($id1);
			chomp($id2);

			$all_ids{$id1} = 1;
			$all_ids{$id2} = 1;

			chomp($str1);
			chomp($str2);

			<fh>; # blank line

			# Remove terminal gaps
			$str1 = uc($str1);
			$str2 = uc($str2);

			if ($str1 =~ /^\-/) {
				$str1 =~ /^(-+)(\S+)/;
				my $gaps_str = $1;
				$str1 = $2;
				$str2 = substr($str2, length($gaps_str), length($str2));
			}

			if ($str2 =~ /^\-/) {
				$str2 =~ /^(-+)(\S+)/;
				my $gaps_str = $1;
				$str2 = $2;
				$str1 = substr($str1, length($gaps_str), length($str1));
			}

			if ($str1 =~ /\-+$/) {
				$str1 =~ /(\-+)$/;
				my $gaps_str = $1;

				$str1 = substr($str1, 0, length($str1) - length($gaps_str));
				$str2 = substr($str2, 0, length($str2) - length($gaps_str));
			}

			if ($str2 =~ /\-+$/) {
				$str2 =~ /(\-+)$/;
				my $gaps_str = $1;

				$str2 = substr($str2, 0, length($str2) - length($gaps_str));
				$str1 = substr($str1, 0, length($str1) - length($gaps_str));
			}

			# Delete single base indels
			my @positions_to_delete;

			while ($str2 =~ m/[ATGC](-)[ATGC]/g) {
				my $p = pos($str2);
				push(@positions_to_delete, $p - 2);
			}

			for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($str1, $item, 1) = "";
				substr($str2, $item, 1) = "";
			}

			my @positions_to_delete;

			while ($str1 =~ m/[ATGC](-)[ATGC]/g) {
				my $p = pos($str1);
				push(@positions_to_delete, $p - 2);
			}

			for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($str1, $item, 1) = "";
				substr($str2, $item, 1) = "";
			}

			# Delete indels larger than 4
			my @positions_to_delete;

			while ($str1 =~ m/[ATGC](-{5,})[ATGC]/g) {
				my $del_start = pos($str1) - length($1);
				push(@positions_to_delete, { del_start => $del_start, del_len => length($1) });
			}

			 for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($str1, $item->{del_start} - 1, $item->{del_len}) = "";
				substr($str2, $item->{del_start} - 1, $item->{del_len}) = "";
			}

			my @positions_to_delete;

			while ($str2 =~ m/[ATGC](-{5,})[ATGC]/g) {
				my $del_start = pos($str2) - length($1);
				push(@positions_to_delete, { del_start => $del_start, del_len => length($1) });
			}

			 for (my $i=@positions_to_delete-1; $i>=0; $i--) {
				my $item = @positions_to_delete[$i];
				
				substr($str1, $item->{del_start} - 1, $item->{del_len}) = "";
				substr($str2, $item->{del_start} - 1, $item->{del_len}) = "";
			}

			# Calc identity
			my @t1 = split(//, $str1);
			my @t2 = split(//, $str2);

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

		close(fh);
	}

	open(fh_dist, ">$opt_o" . "/dist.mat");

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

	#`mkdir $opt_o/jobs`;
	#`mkdir $opt_o/logs`;

	#my $total_sequence_count = `grep '>' $opt_o/targets.ss.FF.C.fa | wc -l`;
	#chomp($total_sequence_count);

	# my $total_job_file_count = ceil($total_sequence_count/$lsf_nb_jobs);
	# my $cmd = $cwd . "utils/split_fasta_file.pl $opt_o/targets.ss.FF.C.fa $total_job_file_count";

	# `$cmd`;

	# `mkdir $opt_o/tmp`;

	# my $db = $opt_o . "/targets.ss.FF.C.fa";

	# my @files = <$opt_o/*.fasta>;
	# my $count = 0;

	# foreach my $file (@files) {
	# 	if (-s $file > 0) {
	# 		$count ++ ;

	# 		my $job_script = "#BSUB -L /bin/bash
	# #BSUB -n 1
	# #BSUB -J oclust_".$count."
	# #BSUB -oo ../logs/$count.log
	# #BSUB -eo ../logs/$count.err
	# #BSUB -q $lsf_queue
	# #BSUB -R rusage[mem=$lsf_memory]
	# #BSUB -W $lsf_time" . ":00\n";

	# 		if ($lsf_account ne "") {
	# 			$job_script .= "#BSUB -P $lsf_account\n";
	# 		}

	# 		$job_script .= "cd $opt_o
	# $cwd";

	# 		$job_script .= "utils/needle_runner.pl $file $db $opt_o/alignments/$count" . ".needle.out $opt_o/tmp $count\n\n";

	# 		open(fh, ">".$opt_o."/jobs/$count".".job");
	# 		print(fh $job_script);
	# 		close(fh);
	# 	}
	# }

	# print("Creating $count job file(s).\n");

	# # Submit jobs and record job identifiers so that we can track them
	# my @files = <$opt_o/jobs/*>;
	# my @ids;

	# foreach my $job (@files) {
	# 	my $cmd = "cd $opt_o/jobs";

	# 	$job =~ /^.+\/(\S+)$/;
	# 	my $fn = $1;

	# 	my $q = "$cmd; bsub < $fn";
	# 	my $out = `$q`;

	# 	$out =~ /^\S+ <(\S+)>/;
	# 	my $job_id = $1;

	# 	push(@ids, $job_id);

	# 	print("Submitted $job\n");
	# }

	# # Write running jobs
	# open(fh, ">" . $opt_o . "/submitted_jobs.txt");

	# foreach my $id (@ids) {
	# 	print(fh "$id\n");
	# }

	# close(fh);
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
	my $f = $opt_o . "/targets.ss.FF.C.fa";
	my $cmd = $cwd."bin/cmalign --cpu $opt_p -o $opt_o" . "/infernal.sto $cwd" . "bin/RDPinfernalTraindata/bacteria16S_508_mod5.stk.cm $f 2>/dev/null >/dev/null";

	`$cmd`;

	# Convert to fasta and remove sequences without any homologous positions
	my $in = Bio::AlignIO->new(-file => $opt_o."/infernal.sto", '-format' => 'stockholm');
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

	my $os = Bio::SeqIO->new(-file => ">" . $opt_o . "/infernal.F.fasta", -format => "fasta");

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

	my $cmd = "$cwd/bin/R/bin/R.fixed --no-save --no-restore --args $cwd $opt_o"."/infernal.F.fasta $opt_o $hclust_algorithm MSA < $cwd/utils/hclust_fr_aln.R";
	`$cmd`;

	print("oclust running in MSA-mode has finished. Results are in:\n$opt_o\n");
}