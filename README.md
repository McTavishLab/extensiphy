# RapUp
### Overview

RapUp is a pipeline that assembles homologous loci by aligning reads to a reference from a multiple sequence alignment, calls consensus, adds to the existing alignment, and places the new lineages in a phylogeny using EPA in RAxML.

RapUp takes an alignment, a tree, and sets of sequencing reads from query taxa.

### Setup and Use

RapUp now takes inputs in the commandline without requiring a config file.
RapUp allows for control over both how many RapUp runs happen
in parallel and how many threads are allocated to each RapUp run
Make sure you dont ask your computer to work too hard by adding more runs and threads than your computer can handle
find out how many cores you have available and calculate (cores * rapup_runs) you wish to run as the same time
if you have 8 cores available, consider starting 2 runs with 3 threads available to each,
then adjust to your optimum setting.

### First Run
Once you've cloned this repo and installed all dependencies to your PATH, begin here. Dependencies are outlined at the bottom of this readme.

If you only plan on using RapUp to add data to an existing alignment and tree, use the following command:

```bash
./multi_map.sh -a ./testdata/combo.fas -d ./testdata
```

If you plan to generate a starting alignment and tree that you wish to add sequences to, test gon_phyling with this command:

```bash
./gon_phyling.sh -d ./gon_phy_testdata
```

#### IMPORTANT!!

RapUp requires that you limit the loci you include for updating to sequences with lengths of 1000bp or above. This is to protect the read mapping and basecall accuracy.

### Controls and Flags For Use:

 (-a) alignment in fasta format,
 (-d) directory of paired end fastq read files for all query taxa,
 (-t) tree in Newick format produced from the input alignment that you wish to update with new sequences or specify NONE to perform new inference (DEFAULT: NONE),
 (-m) alignment type (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA) (DEFAULT: CONCAT_MSA),
 (-o) directory name to hold results (DEFAULT: creates rapup_run),
 (-r) Selected a reference sequence from the alignment file for read mapping or leave as default and a random reference will be chosen (DEFAULT: RANDOM),
 (-p) number of taxa to process in parallel,
 (-c) number of threads per taxon being processed,
 (-1, -2) suffix (ex: R1.fastq or R2.fastq) for both sets of paired end files (DEFAULTS: R1.fq and R2.fq),
 (-g) output format (CONCAT_MSA or SINGLE_LOCUS_FILES) (DEFAULT: CONCAT_MSA),
 (-s) specify the suffix (.fa, .fasta, etc) (DEFAULT: .fasta),
 (-b) bootstrapping tree ON or OFF (DEFAULT: OFF)

 if using single locus MSA files as input,
 (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv),

### Output Files!
- concatenated file: found in your output folder [OUTDIR]/combine_and_infer/extended.aln
- Phylogeny in newick file format: found in your output folder [OUTDIR]/combine_and_infer/RAxML_bestTree.consensusFULL
- taxon specific intermediate files: found in your output folder [OUTDIR]/[TAXON_NAME]. .sam, .bam and .vcf files can be found in here for any additional analyses.

## If all you have is raw reads and you need to create a starting tree:
Creating a starting tree!
You need a tree and alignment with any number of taxa in order to update these with new taxa.
gon_phyling.sh is a simple pipeline to de novo assemble reads before using parsnp for loci selection and finally phylogenetic inference.

1. Move some fraction of your reads to a new directory for assembly and starting tree inference.
2. run: ./gon_phyling.sh -d [PATH/TO/NEW/READ/DIRECTORY]
3. Use the produced alignment file, tree file and the rest of the reads as the inputs for a full RapUp run by running multi_map.sh -a [PATH/TO/ALIGNMENT/FILE] -d [PATH/TO/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2].

### Dependencies

**Using RapUp is limited to Linux at the moment.** Using Ubuntu will ensure the smoothest performance. If you want to use another distro, you'll have to make sure you install analogous one-liners and all that. You have been warned.

Dependencies (Separate programs you'll need to install):

1. [Python 3](https://www.python.org/)
2. [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
3. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)
4. [Seqtk](https://github.com/lh3/seqtk)
5. [Samtools](http://www.htslib.org/)
6. [Bcftools](http://www.htslib.org/)
7. [Fastx](http://hannonlab.cshl.edu/fastx_toolkit/download.html)
8. [Dendropy](https://dendropy.org/)

Additionally, RapUp comes with an additional pipeline for generating a phylogenetic tree from scratch: **Gon\_phyling**. These programs are not required for running Phycorder itself but Gon\ling can be useful if you have a lot of data and aren't interested in hand selecting the loci/genes you include in your alignment. Gon\_phyling's dependencies are as follows:

1. [PARSNP](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
2. [Spades](https://github.com/ablab/spades)
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
4. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)

### Quick dependency install
Use conda for fastest dependency install.

Add appropriate channels to your conda install:

```bash
conda config --prepend channels conda-forge
conda config --prepend channels bioconda
```

Run this command to add the necessary dependencies to your conda environment:

```bash
conda create -n phycorder samtools hisat2 seqtk bcftools fastx-toolkit dendropy raxml
```

Activate your installation

```bash
conda activate phycorder
```

Conda install recipe on the way.





PARALLEL BRANCH command for multimap data test

You must set up the config file for use after you have tested your install.
The PARALLEL branch allows for control over both how many phycorder runs happen
in parallel (hence the name) and how many threads are allocated to each phycorder run_name
Make sure you dont ask your computer to work too hard by adding more runs and threads than your computer can handle
find out how many cores you have available and calculate (cores*phycorder runs you wish to run as the same time)
if you have 8 cores available, consider starting 2 runs with 3 threads available to each,
then adjust to your optimum setting.
First, use the following command as Phycorder should be able to find the included test datafiles

./multi_map.sh ./sample_phycorder.cfg

It is recommended that you leave the sample_phycorder.cfg file alone so you always have a reference
Make a copy and then alter that for your analyses

When you're ready to load your own data, adjust the variable values in the new config file

IMPORTANT!!
Before you do anything else, make a copy of your read data and move that copy into an empty directory
Run the name_parser.py program on that data with the following options:
name_parser.py -d --newtaxa_dir [PATH/TO/DIRECTORY/OF/READS]

This will rename your reads in a way that is easily parsed by Phycorder

Some specifics about the values to change:
The alignment, tree and read directory options should be pretty obvious.
- The alignment is the one used to generate the tree you're plugging into the program as a starting tree
- The tree is the starting tree that you wish to add taxa to.
- the read directory is the directory of paired end reads that you wish to add to your tree.
If you are renaming your reads with name_parser.py (which you absolutely should. File name schemes are generally terrible)
then leave the r1_tail and r2_tail options alone so they function with the name_parser.py outputs
- phycorder runs dictates how many taxa you want to map at a single time.
  - Say if you have 10 taxa to add to your tree, you can tell phycorder that you want to run 5 mappings at a time in order to be efficient
- Threads should be obvious but interacts with [Phycorder runs] by dictating how many threads are given to each program IN AN INDIVIDUAL Run
  - So if you had 5 runs going at a time and you assigned 4 threads per run, this requires 20 threads in order to run.

Creating a starting tree!
You need a tree and alignment with any number of taxa in order to update these with new taxa.
gon_phyling.sh is a simple pipeline to de novo assemble reads before using parsnp for loci selection and finally phylogenetic inference.

To test gon_phyling.sh, run:
./gon_phyling.sh ./sample_gon_phyling.cfg

!!! If all you have is raw reads and you need to create a starting tree:
1. Run name_parser.py on your reads.
2. Move some fraction of your reads to a new folder for assembly and starting tree inference.
3. Make a copy of the sample_gon_phyling.cfg file and alter the variable for the read directory to be assembled.
4. run: ./gon_phyling.sh ./[GON_PHYLING_CONFIG_FILE] in the Phycorder directory.
5. Use the produced alignment file, tree file and the rest of the reads as the inputs for a full Phycorder run.


!!! Installation of programs necessary for Phycorder:
Python packages:
    	Dendropy 4.0 (pip install dendropy)
Software in path for multi_map.sh rapid-updating:
        hisat2 https://ccb.jhu.edu/software/hisat2/index.shtml
        fastx  http://hannonlab.cshl.edu/fastx_toolkit/download.html
        raxmlHPC http://sco.h-its.org/exelixis/web/software/raxml/index.html
        seqtk https://github.com/lh3/seqtk
        samtools / bcftools
        NOTE: requires samtools and bcftools 1.0 - not currently avail via apt-get. Install from http://www.htslib.org/
        Installs nicely but to /usr/local unlike apt-get - make sure paths are correct!
Software in path for gon_phyling.sh:
	parsnp https://harvest.readthedocs.io/en/latest/content/parsnp.html
	spades http://spades.bioinf.spbau.ru/
	bbmap https://sourceforge.net/projects/bbmap/
	raxmlHPC http://sco.h-its.org/exelixis/web/software/raxml/index.html




Old README info
---------------------------------------------------------------------

Pipline that takes an alignment, a tree, and set of sequencing reads form a query taxon.


Assembles  homologous locus, if possible, using closest match as reference,
calls consensus, adds to existing alignment, and places in the tree using EPA in RAxML.


Example run:

    ./map_to_align.sh -a example.aln -t tree.tre -p SRR610374_1.fastq -e SRR610374_2.fastq

NOTE: This example requires downloading SRR610374_1.fastq and SRR610374_2.fastq from  
http://www.ebi.ac.uk/ena/data/view/SRR610374&display=html  
or   
http://www.ncbi.nlm.nih.gov/sra/?term=SRR610374  


##Arguments:
###required
 -a alignment in DNA fasta format. Use 'preprocessing.py' if not available as fasta  (required)  
 -t tree in newick format. Tip lables must match alignment labels, and polytomies must be resolved. Use 'preprocessing.py' to do so if necessary.
 (required)  
 -p paired end query reads fastq (read 1)
 -e paired end query reads fastq (read 2)
    or  
 -s query reads (stub of single end read names. File should be named stub.fastq)  
 (required)  
###optional arguments   
 -o output directory. Optional default is phycorder_run  
 -n run_name.  Optional, default is QUERY.  
 -r Boolean. Align and place reads. Default is 0, set to 1 if you want to align and place reads. Can be SLOWWWW if lots map.  
 -m Boolean. Map reads. Default is 1, set to 0 only if you have already mapped reads.  
 -b Boolean. Align reads to best reference in alignment, call and place consensus sequence. Default is 1, set to 0 if you only want to align and place reads.  
 -w Boolean. Try creating consensus from non optimal (worse) reference locus. Useful for investigating effects of refence choice.  

##Output files:
 (with -n QUERY)  
  ref_nogap.fas : The refence alignement will all gap characters removed, used as potential mapping reference loci  
  {}.bt : bowtie2 refences  
  full_alignment.sam : reads mapped to all loci in refence alignement (same full_sorted.bam, full_sorted.bam.bai)  
  mapping_info : listing of how many reads mapped to each locus. Used in determining best refence locus  
  matches : list of reads that mapped to any locus  
  matches.fq : fastq of all reads mapped to any locus  
  matches.fa : fasta of all reads mapped to any locus  
  matches_unique.fa : matches.fa with duplicate reads removed  
  {alignname}.phy : reference alignement in phylip format for PaPaRa  
  papara_{}.reads : PaPaRa read alignment output files  
  papara_alignment.reads : Alignment of reads in fasta  
  RAxML_{}.QUERY_reads_EPA : RAxML reads output files  


##Requirements:

runs on linux, probably not anywhere else  
(Too many probably...)   
Python packages:
    Dendropy 4.0 (pip install dendropy)  
Software in path:
	bowtie2  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml  
	fastx  http://hannonlab.cshl.edu/fastx_toolkit/download.html  
	raxmlHPC http://sco.h-its.org/exelixis/web/software/raxml/index.html  
	seqtk https://github.com/lh3/seqtk  
	samtools / bcftools  
	NOTE: requires samtools and bcftools 1.0 - not currently avail via apt-get. Install from http://www.htslib.org/  
	Installs nicely but to /usr/local unlike apt-get - make sure paths are correct!  
