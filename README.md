# Phycorder

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
