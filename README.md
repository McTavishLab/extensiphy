# RapUp
### Overview

RapUp assembles homologous loci by mapping to a reference from the alignment, calls consensus, adds to existing alignment, and places in the tree using EPA in RAxML.

RapUp takes an alignment, a tree, and sets of sequencing reads from query taxa.

### Setup and Use

Currently you must set up the config file for use after you have tested your install.
RapUp allows for control over both how many RapUp runs happen
in parallel and how many threads are allocated to each RapUp run
Make sure you dont ask your computer to work too hard by adding more runs and threads than your computer can handle
find out how many cores you have available and calculate (cores*RapUp runs you wish to run as the same time)
if you have 8 cores available, consider starting 2 runs with 3 threads available to each,
then adjust to your optimum setting.

### First Run

If you plan to generate a starting alignment and tree that you wish to add sequences to, test gon_phyling with this command:

```bash
./gon_phyling.sh ./sample_gon_phyling.cfg
```

If you only plan on using RapUp to add data to an existing alignment and tree, use the following command as RapUp should be able to find the included test datafiles

```bash
./multi_map.sh ./sample_rapup.cfg
```

It is recommended that you leave the sample_RapUp.cfg file alone so you always have a reference
Make a copy and then alter that for your analyses

When you're ready to load your own data, adjust the variable values in the new config file

#### IMPORTANT!!
Before you do anything else, make a copy of your read data and move that copy into an empty directory
Run the name_parser.py program on that data with the following options:

```bash
name_parser.py -d --newtaxa_dir [PATH/TO/DIRECTORY/OF/READS]
```

Additionally, limit the loci you include for updating to sequences with lengths of 1000bp or above. This is to protect the read mapping and basecall accuracy.

This will rename your reads in a way that is easily parsed by RapUp

### Some specifics about the values to change:

- align: The alignment is the one used to generate the tree you're plugging into the program as a starting tree. Depending on the type of input, you'll need to adjust the align_type option.
- align_type: set the input type for the alignment. Can be a parsnp.xmfa file output by parsnp, a directory of seperate loci multiple sequence alignment files in the fasta format or a single concatenated multiple sequence alignment in the fasta format.
- tree: The tree is the starting tree that you wish to add taxa to. Set it to "PATH/TO/TREE" if you desire to update the tree. Set to "NONE" if you want a new tree inference.
- read_dir: the read directory is the directory of paired end reads that you wish to add to your tree.
If you are renaming your reads with name_parser.py (which you absolutely should. File name schemes are generally terrible) then leave the r1_tail and r2_tail options alone so they function with the name_parser.py outputs.
- output_type: Dictates if you want to only output a concatenated MSA fasta file or if you want to also output updated seperate loci fasta files as well. CURRENTLY ONLY USE CONCATENATED OPTION!
- loci_positions: a csv file that delimits loci lengths in a concatenated fasta MSA. This will produced for your automatically if you input single locus MSAs or you can make one for your concatenated MSA fasta by using the format found in example_positions.csv in /testdata
- single_locus_suffix: for use when inputting sepereate loci files.
- outdir: the directory name to create and store your results in
- bootstrapping: dictates whether you want a bootstrapped or single best tree. greatly affects speed. 
- RapUp_runs dictates how many taxa you want to map at a single time.
  - Say if you have 10 taxa to add to your tree, you can tell RapUp that you want to run 5 mappings at a time in order to be efficient
- Threads should be obvious but interacts with [RapUp runs] by dictating how many threads are given to each program IN AN INDIVIDUAL Run
  - So if you had 5 runs going at a time and you assigned 4 threads per run, this requires 20 threads in order to run.


### Output Files!
- concatenated file: found in your output folder [OUTDIR]/combine_and_infer/extended.aln
- Phylogeny in newick file format: found in your output folder [OUTDIR]/combine_and_infer/RAxML_bestTree.consensusFULL
- taxon specific intermediate files: found in your output folder [OUTDIR]/[TAXON_NAME]. .sam, .bam and .vcf files can be found in here for any additional analyses.

Creating a starting tree!
You need a tree and alignment with any number of taxa in order to update these with new taxa.
gon_phyling.sh is a simple pipeline to de novo assemble reads before using parsnp for loci selection and finally phylogenetic inference.

To test gon_phyling.sh, run:
```bash
./gon_phyling.sh ./sample_gon_phyling.cfg
```

## If all you have is raw reads and you need to create a starting tree:
1. Run name_parser.py on your reads.
2. Move some fraction of your reads to a new folder for assembly and starting tree inference.
3. Make a copy of the sample_gon_phyling.cfg file and alter the variable for the read directory to be assembled.
4. run: ./gon_phyling.sh ./[GON_PHYLING_CONFIG_FILE] in the RapUp directory.
5. Use the produced alignment file, tree file and the rest of the reads as the inputs for a full RapUp run.


### Dependencies

Unfortunately, RapUp requires some dependencies. You know what they say about not reinventing the wheel. We'll walk through the basics of installation and adding the installed programs to your path so that RapUp can use them.

Using RapUp is limited to Linux at the moment. Using Ubuntu will ensure the smoothest performance. If you want to use another distro, you'll have to make sure you install analogous one-liners and all that. You have been warned.

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
