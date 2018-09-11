# Phycorder

DEV BRANCH command for multimap data test

./multi_map.sh -a /shared/phycorder/testdata/phycorder_testdata.fasta -t /shared/phycorder/testdata/RAxML_bestTree.phycorder_testdata.out -p /shared/phycorder/testdata/ -c 4 -o tmp_o1
## NOTE this currently does not run with -c 2

./multi_map.sh -a /shared/phycorder/example.aln -t /shared/phycorder/tree.tre -p /shared/phycorder/tmp_reads/ -c 4 -o tmp_ncbi

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
