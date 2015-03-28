# Requirements
* Linux
* Perl 5.x
* A computing cluster running LSF if you plan to run the Needleman-Wunsch pairwise alignments.

# Input files
The only input file to oclust is a file in FASTA format containing the sequencing reads to be
clustered. There is no support for FASTQ format, etc. If your data is in FASTQ, you can to FASTA
for example using the 'fastq_to_fasta' script:

```
   $ cd utils
   $ chmod +x fastq_to_fasta
   $ fastq_to_fasta my_fastq.fq my_fasta.fasta
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

   ```
   $ ./oclust_pipeline.pl -x <method> -f <input fasta file> -o <output directory> -p 1 -minl 400 -maxl 1000

   General settings:
   -x NW or MSA               Can be NW for Needleman-Wunsch or MSA for Infernal. [MSA]
   -f [string]                Input fasta file.
   -o [string]                Name of output directory (must not exist) and use full path.
   -R hmm or BLAST            Method to use for reverse complementing sequences. [hmm]
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
   ```

4. Run the last step in the pipeline:

   ```
   $ ./oclust_finalize.pl -i <input directory>

   Settings:
   -i [string]                  Name of the output directory of `oclust_pipeline.pl'.
   -a [clustering algorithm]    Can be one of: complete, average, or single. [complete]
   ```

5. If step 5 fails this is likely due to insufficient memory on the node. Try requesting more memory and run it again.

6. If everything went fine there should now be files with the extension `hclust` in the specified directory. These files are space-delimited and contain two columns (with header). The first column contains the identifier of the sequencing read and the second column contains the OTU (cluster) designation. Clustering is performed on four distances (0.01 to 0.04). If other distances are desired this can be achieved by modifying `hclust.R` in the `utils` directory.

# Dependencies
The oclust pipeline glues together the following open source components:

* EMBOSS, version 6.6.0 was downloaded from ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz,
and compiled using: `./configure --prefix=~/local/emboss --disable-shared --without-mysql`.
* R, compiled with: `./configure --prefix=~/R/ --enable-static=yes --with-x=no --with-tcltk=no`.
* Perl.
* BioPerl.
* NCBI BLAST.
* amos.
* uchime (public domain version).
* hmmer (hmmscan)
* vrevcomp
* beer.

# Help/suggestions
* p.oscar.franzen at gmail.com
