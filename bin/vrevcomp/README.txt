V-RevComp: An open-source, high-throughput software tool to detect and re-orient reverse complementary small and large subunit ribosomal RNA gene sequences of bacterial, archaeal, and fungal origin

Source code available at:
http://www.microbiome.ch/web/Tools.html

Version:	1.1
Purpose:	Detect and re-orient reverse complementary SSU and LSU rRNA sequences
Copyright:	Hartmann et al. 2010
Contact:	Martin Hartmann (contact(at)microbiome.ch)
Affiliation:	Swiss Federal Research Institute WSL, Forest Soils and Biogeochemistry, 8903 Birmensdorf, Switzerland
Programmer: 	Charles Howes (vrevcomp(at)ch.pkts.ca)

The full installation instructions can be found in Users' Guide
Quick installation instructions can be found below.


1) Perl and HMMER

Perl is available at http://www.perl.org/

Installing HMMER version 3
note: V-RevComp 1.1 requires HMMER version 3 and will not work with version 2.

 - download			http://hmmer.janelia.org/#download  
 - uncompress/unpack		tar zxf hmmer-3.0.tar.gz
 - move to new directory	cd hmmer-3.0
 - configure			./configure
 - build			make
 - automated tests		make check
 - automated install		make install

Precompiled binaries are available for some operating systems. They can be found in the "binaries/" directory after the unpacking the tarball. Installing the package takes nothing more than moving these binaries wherever you want them (e.g. /usr/local/bin).


2) V-RevComp

 - download	at http://www.microbiome.ch/web/Tools.html
 - unpack	"unzip vrevcomp"


3) Input file and running V-RevComp

Use sequences in FASTA formatted input file. Copy FASTA file to directory containing vrevcomp.pl. Move into this directory with "cd path/". Typing "perl vrevcomp.pl" lists all options of the program.

Example: perl vrevcomp.pl -h HMMs/SSU/bacteria/ -o out.fasta  -c out.csv -r V1,V2 in.fasta
 -- this will screen the orientation of bacterial 16S sequences from the file in.fasta only based on the regions V1 and V2, and save the orientation-corrected sequences to out.fasta. File out.csv will contain detection statistics for further evaluation where necessary.
