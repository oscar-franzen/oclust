# Requirements
* Linux
* Perl 5.x
* A computing cluster running LSF if you plan to run the pairwise alignments.

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

   First time executed, `oclust_pipeline.pl` will download the human genome sequence and
   format it.

   ```
   $ ./oclust_pipeline.pl -x <method> -f <input fasta file> -o <output directory> -p 1 -minl 400 -maxl 1000

   General settings:
   -x PW or MSA               Can be PW for pairwise alignments (based on Ngila)
                               or MSA for multiple sequence alignment (based on
                               Infernal). [MSA]
   -a complete, average or    The desired clustering algorithm. [complete]
      single    
   -f [string]                Input fasta file.
   -o [string]                Name of output directory (must not exist) and use full path.
   -R HMM, BLAST, or none     Method to use for reverse complementing sequences. [HMM]
   -p [integer]               Number of processor cores to use for -R BLAST and -x MSA. [4]
   -minl [integer]            Minimum sequence length. [optional]
   -maxl [integer]            Maximum sequence length. [optional]
   -rand [integer]            Randomly sample a specified number of sequences. [optional]
   -human Y or N              If 'Y'es, then execute BLAST-based contamination
                               screen towards the human genome. [Y]
   -chimera Y or N            Run chimera check. Can be Y or N. [Y]

   LSF settings (only valid for -x PW):
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
   ```

5. If step 5 fails this is likely due to insufficient memory on the node. Try requesting more memory and run it again.

6. If everything went fine there should now be files with the extension `hclust` in the specified directory. These files are space-delimited and contain two columns (with header). The first column contains the identifier of the sequencing read and the second column contains the OTU (cluster) designation. Clustering is performed on four distances (0.01 to 0.04). If other distances are desired this can be achieved by modifying `hclust.R` in the `utils` directory.

# Dependencies
The oclust pipeline bundles together the following open source/public domain programs:

* R [http://www.r-project.org], compiled with: `./configure --prefix=~/R/ --enable-static=yes --with-x=no --with-tcltk=no`
* The seqinr R package [http://cran.r-project.org/web/packages/seqinr/]
* Perl and BioPerl [http://www.bioperl.org]
* NCBI BLAST [ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/]
* uchime (public domain version) [http://drive5.com/uchime/uchime_download.html]
* HMMER (hmmscan) [http://hmmer.janelia.org/]
* vrevcomp [http://www.microbiome.ch/web/Tools.html]
* infernal [http://infernal.janelia.org/]
* Ngila [http://scit.us/projects/ngila/]
* beer [http://en.wikipedia.org/wiki/India_Pale_Ale]

# Help/suggestions
* p.oscar.franzen at gmail.com
