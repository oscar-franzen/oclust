#!/usr/bin/perl -w
# Program name: V-RevComp
# Version: 1.1
# Purpose: Determining orientation of SSU and LSU rRNA gene sequences
# Copyright (c) Hartmann et al. 2010
# Contact: Martin Hartmann (contact(at)microbiome.ch)
# Swiss Federal Research Institute WSL, Forest Soils and Biogeochemistry, 
# 8903 Birmensdorf, Switzerland
# Programmer: Charles Howes (vrevcomp(at)ch.pkts.ca)
# 
# Citation: Hartmann M, Howes CG, Veldre V, Schneider S, Vaishampayan PA,
# Yannarell AC, Quince C, Johansson P, Björkroth, Abarenkov K, Hallam SJ,
# Mohn WW, Nilsson RH (2011). V-RevComp: Automated high-throughput detection
# of reverse complementary 16S ribosomal RNA gene sequences in large
# environmental and taxonomic datasets. FEMS Microbiology Letters 319(2):
# 140-145.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have receive a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>
#
# Description:
# Hidden Markov Models (HMM) can be used to detect conserved regions
# on SSU and LSU rRNA sequences and detection frequency can give a reliable
# measure of their orientation (i.e. sense or reverse complement). Bacterial,
# archaeal, and fungal SSU and LSU HMMs were created previously (see Hartmann
# et al. 2010, Journal of Microbiological Methods 83(2): 250-253 for details).
# These HMMs are located on the query sequence by screening the sequence first
# in its input orientation and then in its reverse complemenatary orientation.
# The ratio of the detection frequency determines the orientation of the
# sequence. The output fasta file contains the input sequences in the correct
# orientation. The comma/separated value file cotnains the detection statistics
# for each sequence and flags sequences showing ambiguous results.

print STDERR "V-RevComp v. 1.1. Copyright (c) Hartmann et al. 2010.\n";
print STDERR "\n";

use strict;
use Getopt::Long;
use POSIX qw(strftime);

# Record the full command line for logging purposes:
my $COMMAND=join(" ",$0,@ARGV);

# Parse command line:
my $OUTFILE;
my $EVALUE=0.01;
my $SCORE=0;
my $HMMDIR="HMMs/bacteria";
my $OUTCSV;
my $INCLUDEHMM;
my $ALPHABET;
my $USEBITSCORE;
my $DEBUG;
my $INCLUDEV9R;
my $REGIONS;
my $CHUNKSIZE=1000;
my $TMPDIR="/tmp"; # Default on linux?

my $result=GetOptions(
  "bitscore"=>\$USEBITSCORE,
  "csv:s"=>\$OUTCSV,
  "regions:s"=>\$REGIONS,
  "debug:i"=>\$DEBUG,
  "evalue:f"=>\$EVALUE,
  "hmmdir:s"=>\$HMMDIR,
  "output=s"=>\$OUTFILE,
  "includeV9r"=>\$INCLUDEV9R,
  "score:f"=>\$SCORE,
  "zize:i"=>\$CHUNKSIZE,
  "tempdir:s"=>\$TMPDIR,
);
my @INFILES=@ARGV;

my @TODELETE;

# Checking for errors and usage message:
if (!@INFILES or !defined $OUTFILE or !$result) {
  die("Usage: perl vrevcomp.pl [-h hmmdirectory] [-b] [-d] [-e evalue] [-s score]
              [-t tempdir] [-i] [-z chunksize] [-h hmmdirectory] [-c csvoutput]
              [-o outputfile] inputfile [inputfiles...]

This program will examine SSU or LSU rRNA gene datasets for the presence of
reverse complementray entries and report and reorient them. It processes
full-length as well as partial SSU and LSU rRNA gene sequences.

  Options:
    -o outputfile: write corrected sequences to FASTA file
    -c csvoutput: write detection information to CSV file
    -h hmmdirectory: the directory containing HMM files named
	V[1-x]leftlong.HMM
	V[1-x]rightlong.HMM
    -r V1,V2,V7 (SSU) or V01,V09,V15 (LSU): only check these regions
    -t tempdir: use this directory for temporary files (default: $TMPDIR)
	the program runs faster when tempdir is a ramdisk.
    -e evalue: set the global evalue threshold (default: $EVALUE)
    -s score: set the global score threshold (default: $SCORE)
    -i include V9right (SSU only) in the calculations (default: not included)
	note: V9right (SSU only) is slightly more prone to false positives
    -z batchsize: run hmmscan on batches of this many sequences (default: $CHUNKSIZE)
    -b use bitscore instead of evalue threshold (only use on or the other)
    -d write hmm regions to debug.csv

    Example:
    perl vrevcomp.pl -h HMMs/SSU/bacteria/ -o out.fasta -c out.csv -r V1,V2 in.fasta
    --  this will screen the orientation of bacterial SSU rRNA gene sequences from the
	file in.fasta only based on the regions V1 and V2, and save the orientation-
	corrected sequences, as well as the other sequences to out.fasta.
	File out.csv will contain detection statistics for further evaluation where necessary.
");
}

# Look for the specified input files:
if (my @m=grep {!-f} @INFILES) {
  if (@m==1) {die("$m[0]: File not found\n");}
  die("These files were not found: @m\n");
}

# Find size of all input files, for the progress bar later:
my $INSIZE=0;
map {(-e $_) && ($INSIZE+= -s $_)} @INFILES;

# Check to see if an output file was specified:
if (!defined $OUTCSV and !defined $OUTFILE) {
  die("Either -c or -o must be specified, or else there's no point in running this program.\n");
}

# Looking for the hmmscan program in the path:
my $HMMSCAN=join(",,,",grep {-f} map {"$_/hmmscan"} split(/:/,$ENV{'PATH'}.":."));
$HMMSCAN=~s/,,,.*//s;
if ($HMMSCAN eq "") {die("Couldn't find 'hmmscan' in your path; is it installed?\n");}

# Getting the list of long.HMM files as a hash:
my $HMMcount=0;
opendir(DIR,$HMMDIR) or die("$HMMDIR: $!");
my %HMMS=map {$HMMcount++;$_=>1} grep {-f && m/long[.]HMM$/i} map {"$HMMDIR/$_"} readdir(DIR);
closedir(DIR);

# Check the file names of the hmm files:
my $badness=0;
foreach my $n (sort keys %HMMS) {
  if ($n!~ m/^(\S+)(left|right)long[.]HMM$/i) {
    warn("The filename '$n' doesn't match the pattern '???(left|right)long.HMM'\n");
    $badness++;
  }
}
if ($badness) {exit 2;}
if ($HMMcount==0) {die("No files named '*.HMM' were found in $HMMDIR");}

# Get the e-values and scores from each HMM file.
# The hmmscan program can be given an e-value and a score cutoff on
# the command line with -e and -T; we're storing this information
# inside the HMM files in a DESC line.  Format:
# DESC  evalue=0.005  score=0
my %HMMevalue=map {$_=>getevalue($_)} keys %HMMS;
my %HMMscore=map {$_=>getscore($_)} keys %HMMS;

# Turn the HMM filenames into more useful names by stripping the
# directory and suffix (in HMMbasenames) and also the left/right
# long/short words (in HMMnames).  Each hash contains the full
# filename(s) of the unstripped files.
my %HMMbasenames=map {basename($_,".HMM")=>$_} keys %HMMS;
my %HMMnames;
map {$HMMnames{basename($_,"(left|right)long[.]HMM")}.=" ".$_} keys %HMMS;

# Error checking to make sure that each region has a left and right
# long hmm:
foreach my $x (keys %HMMnames) {
  my %set=map {basename($_,"[.]HMM")=>1} grep {/./} split(/ /,$HMMnames{$x});
  if (!defined $set{"${x}leftlong"}) {die("${x}leftlong.HMM was not found\n");}
  if (!defined $set{"${x}rightlong"}) {die("${x}rightlong.HMM was not found\n");}
}

# Optimization: sort the keys only once, and store them for later use:
my @HMMS=sort keys %HMMS;

# Exclude V9right unless told to include it:
if (! defined $INCLUDEV9R) { @HMMS=grep {! /V9right/} @HMMS; }
if (defined $REGIONS) {
  $REGIONS=~s/,/|/g;
  $REGIONS="($REGIONS)(left|right)(short|long)";
  @HMMS=grep {/$REGIONS/} @HMMS;
}

# Final summary:
my %summary;
my %CACHE;

# Prepare the output files:
if (defined $OUTCSV) { open(OUTCSV,">$OUTCSV") or die("$OUTCSV: $!"); print OUTCSV "Name,Orientation,Forward count,Reverse count\n"; }
if (defined $OUTFILE) {open(OUTFILE,">$OUTFILE") or die("$OUTFILE: $!");}

# Progress bar variables:
my $last=time();
my $start=time();
my $offset=0;
$last=progress($offset,$INSIZE,$start,$last);
open(PROGRESS,">progress.csv") or die("progress.csv: $!\n"); # CGH debug

# Main loop:
foreach my $f (@INFILES) {
  print "Processing file $f\n";
  open(IN2,"<$f") or die("$f: $!\n");
  open(OUT,">temp.chunk.$$") or die("temp.chunk.$$: $!\n");
  my $zcount=0;
  do {
    local $/="\n>";
    while (my $s=<IN2>) {
      $s=~s/^>?/>/s; # Fix start
      $s=~s/\n>/\n/s; # Fix end
      print OUT $s;
      $zcount++; # Number of sequences written
      $last=progress($offset+tell(IN2),$INSIZE,$start,$last);
      if ($zcount>=$CHUNKSIZE) {
        $zcount=0;
        close(OUT);
	precalculate("temp.chunk.$$",\@HMMS);
	process_file("temp.chunk.$$");
	%CACHE=(); # Empty cache
        open(OUT,">temp.chunk.$$") or die("temp.chunk.$$: $!\n");
      }
    }
  };
  close(OUT);
  if (-e "temp.chunk.$$" and !-z "temp.chunk.$$") {
    precalculate("temp.chunk.$$",\@HMMS);
    process_file("temp.chunk.$$");
    $last=progress($offset+tell(IN2),$INSIZE,$start,$last);
  }
  close(IN2);
  unlink("temp.chunk.$$");
}
close(PROGRESS); # CGH debug

# Close everything that's open:
if (defined $OUTFILE) {close(OUTFILE);}
if (defined $OUTCSV) {close(OUTCSV);}

# Blank line at the end, for progress bar cleanup:
print STDERR "\n";

# Summary:
print "Summary:\n";
printf "%-20s  %s\n","Orientation","Count";
map {printf "%-20s  %s\n",$_,$summary{$_}} sort keys %summary;

# End:
exit 0;

#----------------------------------------------------------------------
# Do the HMM searches ahead of time, as a batch:
#  
# Arguments: precalculate($file,$hmmsref)
#   $file: which file to scan
#   $hmmsref: ref to array of HMM files to use
sub precalculate {
  my ($filename,$hmmsref)=@_;

  # Adapt to newline problems:
  open(MAININPUT,"<$filename") or die("$filename: $!\n");
  local $/="\n>"; # read input one sequence at a time, not one line at at time

  # Create reverse file:
  my $rev="$TMPDIR/temp.rev.$$";
  open(OUTREV,">$rev") or die("$rev: $!\n");
  while(my $sequence=<MAININPUT>) {
    $sequence=~s/>$//s; $sequence=~s/^>?/>/s; # Fix '>' symbol at start and end
    if ($sequence!~m/^(>[^\n]+\n)(.*)$/s) {
      die("Bad sequence found in $filename, sequence: $sequence\n");
    }
    print OUTREV $1,revcomp($2),"\n";
  }
  close(OUTREV);
  close(MAININPUT);

  # Scan them:
  scanfile($filename,"fwd",$hmmsref);
  scanfile($rev,"rev",$hmmsref); unlink($rev);
}

# This runs hmscan on the given input file.
#
# Arguments: scanfile ($file,$label,$hmmsref)
#   $file: which file to scan; should contain $CHUNKSIZE sequences
#   $label: 'fwd' or 'rev' to indicate direction
#   $hmmsref: a reference to a list of hmm files to run
sub scanfile {
  my ($file,$label,$hmmsref)=@_;

  # Check each file with each hmm:
  for (my $hcount=0;$hcount<@$hmmsref;$hcount++) {
    my $hmm=${$hmmsref}[$hcount];
    # Run hmmpfam:
    # --domE thresh = evalue threshold
    my $criteria="--domE $HMMevalue{$hmm}";
    if (defined $USEBITSCORE) {
      # --domT score = bit score threshold (can be negative)
      my $criteria="--domT $HMMscore{$hmm}";
    }

    # Run hmmscan on the normal input files:
    my $cmd="$HMMSCAN --max $criteria $hmm $file";
    my $tempout="$TMPDIR/temp.out.$$.".rand();
    my $output=`$cmd 2>&1 > $tempout`;
    if (!defined $output) {die("$cmd: $!\n");}
    if ($output =~ m/cannot execute binary file|Exec format error|command not found/) {
      print "Error!  cmd=$cmd\n";
      print "output:\n$output";
      exit 2;
    }

    # Read the output from hmmscan:
    open(IN,"<$tempout") or die("$tempout: $!\n");
    local $/="\n//\n";
    while(my $aa=<IN>) {
      if ($aa!~m/Query: +(\S+)/s) {die("$hmm $label: $aa");}
      my $k="$hmm|$1|$label";
      $k=~s/^.*\///;
      $CACHE{$k}=process_query($aa);
    }
    close(IN);
    unlink($tempout);
  }
}

# Process a file:
#
# Arguments: process_file($file)
#   $file: The file to process
sub process_file {
  my $file=$_[0];
  if (!defined $file) {die("process_file called with undefined filename");}
  if (!-e $file) {die("$file doesn't exist.\n");}
  if (!-f $file) {die("$file isn't a file.\n");}
  open(IN,"<$file") or die("$file: $!\n");
  do {
    local $/="\n>";  # Record separator is '\n>'
    while (<IN>) {
      my $v=$_;
      $v=~s/^>?/>/s; $v=~s/\n>?$/\n/s; # Fixup endings
      process_sequence($v);
    };
  };
  close(IN);
}

# Process a sequence:
#
# Arguments: process_sequence($seq)
#   $seq: a header+sequence from an input fasta file
sub process_sequence {
  my ($seq)=@_;

  # For each query sequence:
  #
  # (1) check the default orientation for HMMER matches against all 18(?) long
  # HMMs of V-Xtractor. Count the number of such hits (0-18?).
  #
  # (2) compute the reverse complement of the sequence - trivial
  #
  # (3) check the recomputed orientation for HMMER matches against all 18 (?)
  # long HMMs of V-Xtractor. Count the number of such hits (0-18?).
  #
  # (4) decision: if only one of the two orientations have matches, then you
  # know which orientation the query comes in
  #
  # (5) decision: if you have hits on both strands, then mark the query as
  # uncertain. You could use the two counts to find out the most probable
  # direction and report it to the user.

  # Clean and separate input:
  $seq=~s/\r/\n/g;
  $seq=~s/^>?([^\n]+)\n//; # Remove name/description from sequence
  my $nameseq=$1;             # Store it here
  $nameseq=~s/ .*//;       # Delete spaces and any other comments
  $seq=~tr/\n\t //d;       # Delete newlines, spaces and tabs

  # (1) check the default orientation for HMMER matches against all 18(?) long
  # HMMs of V-Xtractor. Count the number of such hits (0-18?).
  my %fresult=map {$_=>matchhmm($_,"$nameseq|fwd")} @HMMS;
  my $fcount=grep {!/notfound/} values %fresult;

  # (2) deleted

  # (3) check the reverse orientation for HMMER matches against all 18 (?)
  # long HMMs of V-Xtractor. Count the number of such hits (0-18?).
  my %rresult=map {$_=>matchhmm($_,"$nameseq|rev")} @HMMS;
  my $rcount=grep {!/notfound/} values %rresult;
  #print "$nameseq|rev: ",map {"$_=>$rresult{$_}\n"} keys %rresult;

  # (4) decision: if only one of the two orientations have matches, then you
  # know which orientation the query comes in
  #
  # (5) decision: if you have hits on both strands, then mark the query as
  # uncertain. You could use the two counts to find out the most probable
  # direction and report it to the user.

  my $orientation="notfound";
  if ($fcount==0 and $rcount!=0) { $orientation="reverse"; }
  if ($fcount!=0 and $rcount==0) { $orientation="forward"; }
  if ($fcount!=0 and $rcount!=0) {
    $orientation="uncertain-eitherway?";
    if ($fcount>$rcount) { $orientation="uncertain-forward?"; }
    if ($fcount<$rcount) { $orientation="uncertain-reverse?"; }
  }

  #if ($orientation eq "notfound") {print "notfound: fcount=$fcount rcount=$rcount\n";}
  
  # Summarize the results in various ways and provide a "ready-to-go" FASTA
  # file. It would operate on 16S sequences of any size: pyro up to full-length.
  $summary{$orientation}++;
  if (defined $OUTCSV) {print OUTCSV "$nameseq,$orientation,$fcount,$rcount\n";}
  if (defined $OUTFILE) {
    if ($orientation !~ m/reverse/) {
      print OUTFILE ">$nameseq\n$seq\n";
    } else {
      # compute the reverse complement of the sequence - trivial
      my $rseq=revcomp($seq);
      print OUTFILE ">$nameseq-reverse\n$rseq\n";
    }
  }
}

# Retrieve an HMM search.
sub matchhmm { # Hmm=filename, Name=sequencename
  my ($hmm,$name)=@_;
  $name=~s/ .*//;

  my $k="$hmm|$name"; # The "|fwd" or "|rev" was added to the name earlier
  $k=~s/^.*\///;
  my $output=$CACHE{$k};
  if (!defined $output) {
    die("matchhmm($hmm,$name)=undefined");
  }
  return $output;
}

# Parse the output of a single hmm against a single sequence
#
# Arguments: process_query($output)
#   $output: the output of hmmscan
sub process_query { 
  my $output=$_[0];

  my $oldout=$output;

  # Error checking:
  if (!defined $output) {die("process_query called with no output");}

  if ($output=~m/(FATAL:[^\n]*)/s) {die("\nprocess_query: \n$1\n");}
  if ($output=~m/No individual domains that satisfy reporting thresholds/) {return "notfound";}
  if ($output=~m/No hits detected that satisfy reporting thresholds/) {return "notfound";}
  if ($output=~m/No targets detected that satisfy reporting thresholds/) {return "notfound";}

  # Data extraction:
  $output=~s/\r//g; # Delete carriage returns (fine for windows and unix)

  # Hmmer2:
  #$output=~s/.*\nParsed for domains:\nModel[^\n]*\n----[- ]*\n//s; # Cut
  #$output=~s/\nAlignments of top-scoring domains:.*//s; # Cut
  #my ($start,$end)=$output=~m/^\S+\s+\S+\s+(\d+)\s+(\d+)/s; # Match it

  # Hmmer3:
  # Always present, even when no matches found:
  $output=~s/.*\nDomain annotation for each model \(and alignments\):\n//s;
  $output=~s/\nInternal pipeline statistics summary:.*//s;
  # Only present when a match found:
  $output=~s/>>[^\n]*\n   #[^\n]*\n ---[- ]*\n//s;
  $output=~s/\n  Alignments for each domain:.*//s;

  # Take the 'alifrom/ali to' values from the best domain:
  my @outlist=grep {/./} split(/\n/,$output);
  my $bestoutput=$outlist[0];
  my $bestev=10;
  for (my $x=0;$x<@outlist;$x++) {
    my $ev=substr($outlist[$x],31,9);
    $ev=~s/\s+//g;
    $ev=$ev+0;
    if ($ev !~ m/^[\d.eE-]+$/) {
      die("\nMangled output: ev=($ev)\n".join("\n",@outlist)."\n\nFrom:\n$oldout");
    }
    if ($ev+0<$bestev) {
      $bestev=$ev;$bestoutput=$outlist[$x]
    }
  }
#   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
# ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#   1 !   31.3   0.0   1.1e-11   1.1e-11       4      38 ..       5      39 ..       2      41 .. 0.93
  my $values=substr($bestoutput,78,15); # Doing envfrom-env to
  my ($start,$end)=$values=~m/^\s*(\d+)\s+(\d+)$/s; # Match it
  if (defined $DEBUG) {
    $values=substr($bestoutput,20,9);
    my ($evalue)=$values=~m/^\s*(\S+)\s*$/s; # Match it
    $values=substr($bestoutput,7,6);
    my ($bitscore)=$values=~m/^\s*(\S+)\s*$/s; # Match it
    open(DEBUG,">>debug.csv");
    print DEBUG "$evalue,$bitscore,$start,$end\n".join("\n",@outlist)."\n";
    close(DEBUG);
  }

  # Return something if not found:
  if (!defined $start) {return "notfound";}
  # Shouldn't happen:
  if ($start>=$end) {die("\nhmmpfam reported a match starting at $start but ENDING at $end (wrong order)");}
  # Success!
  return "$start-$end"; # = the region that was found
}

# Extract the e-value that's been hardcoded into the HMMs:
sub getevalue {
  my $hmm=$_[0];
  open(HMMIN,"<$hmm") or die("$hmm: $!");
  my @l=grep {/^DESC\s(.*\s)?evalue=[0-9.e-]+/i} <HMMIN>;
  close(HMMIN);
  if (@l) {return $l[0]=~m/\sevalue=([0-9.e-]+)/i;} # First evalue
  return $EVALUE;  # Default value
}

# Extract the score that's been hardcoded into the HMMs:
sub getscore {
  my $hmm=$_[0];
  open(HMMIN,"<$hmm") or die("$hmm: $!");
  my @l=grep {/^DESC\s(.*\s)?score=[0-9.-]+/i} <HMMIN>;
  close(HMMIN);
  if (@l) {return $l[0]=~m/\sscore=([0-9.e-]+)/i;} # First score
  return $SCORE;  # Default value
}

# Remove the directory name from the string:
sub basename {
  (my $aa=$_[0])=~s/.*\///;
  if (defined $_[1]) { $aa=~s/${_[1]}$//i; } # Trim extension if given
  return $aa;
}

# Progress bar routine:
sub progress {
  my ($pos,$total,$stime,$ltime)=@_;
  if ($ltime==time()) {return time();} # No change
  if ($stime==$ltime) { # At start
    print STDERR "Initializing progress bar...";
    return time();
  }
  $ltime=time();
  my $done=$pos/($total||$pos||1); # Avoid div by zero
  print PROGRESS time(),",",$pos,",",$done,"\n"; # CGH debug
  my $ttime=($ltime-$stime)/$done;
  print STDERR "\r";
  print STDERR strftime("Start:%T",localtime($stime))." ";
  print STDERR sprintf("[%-34s]",("="x(int($done*30))).int($done*100)."%");
  print STDERR strftime(" End:%T",localtime($stime+$ttime))."  ";
  if ($pos!=$total) {
    print STDERR "Left:".nicetime($stime+$ttime-time());
  } else {
    print STDERR "Total:".nicetime($ttime)."\n";
  }
  return $ltime;
}

# Format a number of seconds as hour:minute:second.  Can't handle
# negative numbers, and doesn't care about days.
sub nicetime {
  my $d=$_[0];
  my $hour=int($d/3600); $d-=$hour*3600;
  my $min=int($d/60); $d-=$min*60;
  my $sec=int($d);
  return sprintf("%d:%02d:%02d",$hour,$min,$sec);
}

# Reverse complement a nucleotide sequence with or without IUPAC ambiguity codes.
sub revcomp {
  my $r=reverse($_[0]);
  $r=~s/\s//gs;
  $r=~y/abcdghkmrstvwnABCDGHKMRSTVWN[]/tvghcdmkysabwnTVGHCDMKYSABWN][/;  # IUPAC ambiguity codes
  return $r;
}
