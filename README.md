# About oclust.
A pipeline for clustering long 16S rRNA sequencing reads, or any sequences, into
[Operational Taxonomic Units](http://en.wikipedia.org/wiki/Operational_taxonomic_unit).

# Requirements
* Linux v.2.6.x
* Perl v.5.10.1
* R (should be in path)
  * package _seqinr_ should be installed:
```
> install.packages("seqinr")
```


# Input files
The only input file to oclust is a file in FASTA format containing the sequencing reads to be
clustered.

FASTQ files can be converted to FASTA:

```
   $ cd utils
   $ chmod +x fastq_to_fasta.pl
   $ ./fastq_to_fasta.pl file.fastq > file.fasta
```

# Installation
1. Get the repository:

   `$ git clone https://github.com/oscar-franzen/oclust.git oclust`

2. Make executable (might not be necessary):

   ```
   $ cd oclust
   $ chmod +x *.pl
   ```

3. Decide if you want to compute distances based on Needleman-Wunsch or Infernal. The latter will
   be substantially faster.

   First time executed, `oclust_pipeline.pl` will download the human genome sequence and
   format it.

```
   $ ./oclust_pipeline.pl -x <method> -f <input file> -o <output directory> -p <number of CPUs>

   General settings:
   -x PW or MSA               Can be PW for pairwise alignments (based on Needleman-Wunsch)
                               or MSA for multiple sequence alignment (based on
                               Infernal). [MSA]
   -t local or cluster        If -x is PW, should it be parallelized by running it locally
                               on multiple cores or by submitting jobs to a cluster
                               (requires a system with the LSF scheduler). [local]
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
```

# Dependencies
The oclust pipeline bundles together the following open source/public domain software:

* [R](http://www.r-project.org), compiled with: `$ ./configure --prefix=~/R/ --enable-static=yes --with-x=no --with-tcltk=no`
* [The seqinr R package](http://cran.r-project.org/web/packages/seqinr/)
* [Perl](http://www.perl.org) and [BioPerl](http://www.bioperl.org)
* [Parallel::ForkManager](http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm)
* [Memory::Usage](http://search.cpan.org/~doneill/Memory-Usage-0.201/lib/Memory/Usage.pm)
* [NCBI BLAST](http://ftp.ncbi.nlm.nih.gov/blast/executables/release/)
* [uchime (public domain version)](http://drive5.com/uchime/uchime_download.html)
* [HMMER (hmmscan)](http://hmmer.janelia.org/)
* [vrevcomp](http://www.microbiome.ch/web/Tools.html)
* [infernal](http://infernal.janelia.org/)
* [EMBOSS Needleman-Wunsch implementation (needle)](http://emboss.sourceforge.net/), compiled with: `$ ./configure --prefix=~/e/ --disable-shared --without-mysql --without-postgresql --without-axis2c --without-hpdf --without-x --without-pngdriver`

# Contact
* p.oscar.franzen at gmail.com
