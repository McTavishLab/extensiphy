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
./multi_map.sh -a ./testdata/combo.fas -t ./testdata/combo.tre -d ./testdata
```

If you plan to generate a starting alignment and tree that you wish to add sequences to, test gon_phyling with this command:

```bash
./gon_phyling.sh -d ./gon_phy_testdata
```

#### IMPORTANT!!

RapUp requires that you limit the loci you include for updating to sequences with lengths of 1000bp or above. This is to protect the read mapping and basecall accuracy.

### RapUp Controls and Flags For Use:

- (-a) alignment in fasta format,
- (-d) directory of paired end fastq read files for all query taxa,
- (-t) tree in Newick format produced from the input alignment that you wish to update with new sequences or specify NONE to perform new inference (DEFAULT: NONE),
- (-m) alignment type (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA) (DEFAULT: CONCAT_MSA),
- (-o) directory name to hold results (DEFAULT: creates rapup_run),
- (-r) Selected a reference sequence from the alignment file for read mapping or leave as default and a random reference will be chosen (DEFAULT: RANDOM),
- (-p) number of taxa to process in parallel,
- (-c) number of threads per taxon being processed,
- (-1, -2) suffix (ex: R1.fastq or R2.fastq) for both sets of paired end files (DEFAULTS: R1.fq and R2.fq),
- (-g) output format (CONCAT_MSA or SINGLE_LOCUS_FILES) (DEFAULT: CONCAT_MSA),
- (-s) specify the suffix (.fa, .fasta, etc) (DEFAULT: .fasta),
- (-b) bootstrapping tree ON or OFF (DEFAULT: OFF)
- (-i) set whether to clean up intermediate output files to save disk space)(KEEP, CLEAN)(DEFAULT: KEEP)

 if using single locus MSA files as input,
- (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv),

### Output Files!
- concatenated file: found in your output folder [OUTDIR]/combine_and_infer/extended.aln
- Phylogeny in newick file format: found in your output folder [OUTDIR]/combine_and_infer/RAxML_bestTree.consensusFULL
- taxon specific intermediate files: found in your output folder [OUTDIR]/[TAXON_NAME]. .sam, .bam and .vcf files can be found in here for any additional analyses.

### Gon_phyling Controls and Flags For Use

INPUT OPTIONS:
- (-d) directory of paired end reads. All output folders and files will be contained here
- (-g) the name of the genome you wish to use as a reference during loci selection (if any)(DEFAULT: NONE)
- (-1, -2) suffixes of paired-end input files in read directory (DEFAULT: -1 R1.fastq -2 R2.fastq)

 OUTPUT
- (-b) bootstrapping setting. Do you want to perform 100 boostrap replicates and add the support values to the best tree? (DEFAULT: OFF)
- (-o) output type. Output either a concatenated multiple sequence alignment only or also output separate loci alignment files (DEFAULT: LOCI) (OPTIONS: LOCI, LOCUS)
- (-l) Locus position file. Use if selecting -o LOCUS. Outputs a csv file tracking the loci names and their positions within the concatenated MSA (DEFAULT: gon_phy_locus_positions.csv)

 RUNNING PROGRAM
- (-r) gon_phyling runs. This is the number of genomes assembled at a single time (DEFAULT: 2)
- (-c) Threads for each gon_phyling run. Figure out how many cores you have available and input [# of threads x # of parrallel genome assemblies] = cores you can allocate. (DEFAULT: 2)

## If all you have is raw reads and you need to create a starting tree:
Creating a starting tree!
You need a tree and alignment with any number of taxa in order to update these with new taxa.
gon_phyling.sh is a simple pipeline to de novo assemble reads before using parsnp for loci selection and finally phylogenetic inference.

1. Move some fraction of your reads to a new directory for assembly and starting tree inference.
2. run: 
```bash
./gon_phyling.sh -d [PATH/TO/NEW/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2]
```
3. Use the produced alignment file, tree file and the rest of the reads as the inputs for a full RapUp run by running:
```bash
 multi_map.sh -a [PATH/TO/ALIGNMENT/FILE] -d [PATH/TO/READ/DIRECTORY] -t [PATH/TO/TREE/FILE] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2].
```

### Tutorial

For a more indepth walkthrough of how to install dependencies for use with RapUp and how to run RapUp using different data types and options, try the tutorial in the tutorial folder. You can copy code snippets into your terminal window.

### Dependencies

**Using RapUp is limited to Linux at the moment.** Using Ubuntu will ensure the smoothest performance. If you want to use another distro, you'll have to make sure you install analogous one-liners and all that. You have been warned.

Dependencies (Separate programs you'll need to install):

1. [Python 3](https://www.python.org/)
2. [Hisat2](https://daehwankimlab.github.io/hisat2/download/)
3. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)
4. [Seqtk](https://github.com/lh3/seqtk)
5. [Samtools](http://www.htslib.org/)
6. [Bcftools](http://www.htslib.org/)
7. [Fastx](http://hannonlab.cshl.edu/fastx_toolkit/download.html)
8. [Dendropy](https://dendropy.org/)

Additionally, RapUp comes with an additional pipeline for generating a phylogenetic tree from scratch: **Gon\_phyling**. These programs are not required for running RapUp itself but Gon\_ling can be useful if you have a lot of data and aren't interested in hand selecting the loci/genes you include in your alignment. Gon\_phyling's dependencies are as follows:

1. [PARSNP](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
2. [Spades](https://github.com/ablab/spades)
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
4. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)


### Apt-get dependency install
Almost all programs for running RapUp are available with apt-get. Hisat2 is not available with apt-get. Run the commands found below to install:

```bash
apt-get install raxml
apt-get install seqtk
apt-get install samtools
apt-get install bcftools
apt-get install fastx-toolkit
pip install dendropy
```

Installs with apt-get for Gon\_phyling are not currently available. You will have to install these programs manually or with conda.


### Conda dependency install
Use conda for fastest dependency install.

Add appropriate channels to your conda install:

```bash
conda config --prepend channels conda-forge
conda config --prepend channels bioconda
```

Run this command to add the necessary dependencies to your conda environment:

```bash
conda create -n rapup samtools hisat2 seqtk bcftools fastx-toolkit dendropy raxml
```

Activate your installation

```bash
conda activate rapup
```

Conda install recipe on the way.
