# Extensiphy
## Overview

Extensiphy is a pipeline that assembles homologous loci by aligning reads to a reference from a multiple sequence alignment, calls consensus, adds to the existing alignment. Homologous loci may be kept concatenated or split back into individual alignments prior to phylogenetic estimation.

![Extensiphy_worlflow](./doc/images/EP_workflow_2.png)

Extensiphy takes an alignment and sets of sequencing reads from query taxa (a). Reads are aligned to a reference sequence and a consensus sequence is called (b). The new sequence is added to the alignment and the updated alignment is used to estimate a phylogeny (c).

## Setup and Use

### Docker
The simplest and most hassle free way to run Extensiphy is using Docker.
the Quick Install and Run section (LINK) will review the docker installation instructions.

### Anaconda
You can also install the dependencies of Extensiphy using Anaconda. This link will take you to the Anaconda installation guide.

### Advanced
*For advanced users of Linux* If you're comfortable installing programs by hand, this section is for you.

## Quick Install and Run
### Building and testing your own Extensiphy Docker image
First we'll building the Docker image and a container to test your Extensiphy installation. Then we'll connect your data to a new container so you can begin updating your own alignments!

1. Make sure you have [Docker desktop installed](https://www.docker.com/products/docker-desktop). Then, download the [extensiphy repository](https://github.com/McTavishLab/extensiphy) to your computer. One way is to type from the terminal:

```bash
git clone https://github.com/McTavishLab/extensiphy.git
cd extensiphy
```

2. To build your Docker installation of Extensiphy, we'll need to build the Docker image.

```bash
docker build --tag ep_image .
```

3. We'll build your Extensiphy Docker container using this command.
The `-i` flag will make the container interactive and allow you to run Extensiphy
within the container.
```bash
docker run --name ep_container -i -t ep_image bash
```

4. Your command line prompt should change to indicate that you are now working
inside your Extensiphy container. To test your installation, run this command:
```bash
./multi_map.sh -a ./testdata/combo.fas -d ./testdata
```

Once Extensiphy has finished running on the test data, you should see a line saying:
```bash
Alignment file is: /usr/src/app/extensiphy/EP_output/outputs/extended.aln
```
If you did not get this message, you'll have to check output log `ep_dev_log.txt`
to learn more about the issue before proceeding.
For a deeper walkthrough of what has actually happened, you can try the tutorial.
To get right down to business and update your own alignment, continue with the walkthrough.

### Using Extensiphy on your own data

5. Next, you'll need to move the data you want to use to a directory we can link to a new container.
First, let's create a new directory and move the data we want to use into the new directory:
We'll use brackets `[]` to indicate variables you should replace with your own files or paths
```bash
mkdir new_data_dir
mv [/path/to/your/alignment_file] [/path/to/new_data_dir]
mv [/path/to/your/raw_read_files] [/path/to/new_data_dir]
```

6. We'll build a new Extensiphy Docker container and connect the directory containing your data to the container.
Replace the `[stuff inside the brackets]` with the appropriate paths and folder names you've used so far.
```bash
docker run --name ep_container -i -t -v [/path/to/new_data_dir]:/usr/src/app/linked_data ep_image bash
```
7. Now you can run the same command as earlier but we'll specify that the `data`
as where your data is located. The output will be an updated sequence alignment.
You will also need to specify the suffixes of your read files using the
`-1` and `-2` flags.
The `-o` flag lets you specify the name of the output folder.
```bash
./multi_map.sh -a /usr/src/app/linked_data/[alignment_file] -d /usr/src/app/linked_data -1 [suffix_1] -2 [suffix_2] -o [output_dir_name]
```

Once the Extensiphy run is finished, you can check the `outputs` directory
for the updated alignment file.

## Extensiphy Controls and Flags For Use:

### Required flags
- (-a) alignment in fasta format,
- (-d) directory of paired end fastq read files for all query taxa,

### Optional flags
- (-u) produce only an updated alignment or perform full phylogenetic estimation (ALIGN or PHYLO) (DEFAULT: ALIGN),
- (-t) tree in Newick format produced from the input alignment that you wish to update with new sequences or specify NONE to perform new inference (DEFAULT: NONE),
- (-1, -2) suffix (ex: R1.fastq or R2.fastq) for both sets of paired end files. Required if suffix is different than default (DEFAULTS: R1.fq and R2.fq),
- (-m) alignment type (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA) (DEFAULT: CONCAT_MSA),
- (-o) directory name to hold results (DEFAULT: creates rapup_run),
- (-r) Selected a reference sequence from the alignment file for read mapping or leave as default and a random reference will be chosen (DEFAULT: RANDOM),
- (-p) number of taxa to process in parallel,
- (-c) number of threads per taxon being processed,
- (-e) set read-type as single end (SE) or pair-end (PE) (DEFAULT: PE),
- (-g) output format (CONCAT_MSA or SINGLE_LOCUS_FILES) (DEFAULT: CONCAT_MSA),
- (-s) specify the suffix (.fa, .fasta, etc) (DEFAULT: .fasta),
- (-b) bootstrapping tree ON or OFF (DEFAULT: OFF)
- (-i) set whether to clean up intermediate output files to save disk space)(KEEP, CLEAN)(DEFAULT: KEEP)

 if using single locus MSA files as input,
- (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv),
- (-n) Set size of loci size cutoff used as input or output (Options: int number)(DEFAULT: 700)

## Output Files!
- Concatenated alignment file: found in your output folder [OUTDIR]/outputs/extended.aln
- Phylogeny in newick file format (if you selected to output a phylogeny): found in your output folder [OUTDIR]/combine_and_infer/RAxML_bestTree.consensusFULL
- Taxon specific intermediate files (if you kept intermediate files): found in your output folder [OUTDIR]/[TAXON_NAME]. .sam, .bam and .vcf files can be found in here for any additional analyses.

### Gon_phyling Controls and Flags For Use
Additionally, Extensiphy comes with an additional pipeline for generating a
phylogenetic tree from scratch: **Gon\_phyling**.
These programs are not required for running Extensiphy itself but Gon\_ling
can be useful if you have a lot of data and aren't interested in hand selecting
the loci/genes you include in your alignment.

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
3. Use the produced alignment file, tree file and the rest of the reads as the inputs for a full Extensiphy run by running:
```bash
 multi_map.sh -a [PATH/TO/ALIGNMENT/FILE] -d [PATH/TO/READ/DIRECTORY] -t [PATH/TO/TREE/FILE] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2].
```

## Tutorial

We recommend you run through the tutorial before using Extensiphy on your own dataset. The tutorial will walk through how to install program dependencies for use with Extensiphy and how to run Extensiphy using different data types and options. You can find the tutorial file in the tutorial folder (link). You can copy code snippets into your terminal window.

## Installation Methods

### Anaconda Installation
You can install the dependencies of Extensiphy using the Anaconda package manager.
 Install instructions for Anaconda can be found here (link).
Once Anaconda has been installed, use this command to create an environment with
 all of the Extensiphy dependencies added to it.

1. Add appropriate channels to your conda install:

```bash
conda config --prepend channels conda-forge
conda config --prepend channels bioconda
```

2. Run this command to create a new environment (extensiphy_env) and add the necessary dependencies:

```bash
conda create -n extensiphy_env samtools bwa-mem2 seqtk bcftools fastx-toolkit dendropy raxml
```

3. Activate your installation

```bash
conda activate extensiphy_env
```

4. Once all installation processes are complete and you've activated your environment,
you can clone the Extensiphy repository.

```bash
git clone https://github.com/McTavishLab/extensiphy.git
cd extensiphy
```

5. Finally, you can test your installation by running the following command.

```bash
./multi_map.sh -a ./testdata/combo.fas -d ./testdata
```

Once Extensiphy has finished running on the test data, you should see a line saying:
```bash
Alignment file is: [path/to]/EP_output/outputs/extended.aln
```

This indicates that the run completed successfully and you can find the updated
alignment file in the `output` directory.
If you did not get this message, you'll have to check output log `ep_dev_log.txt`
to learn more about the issue before proceeding.


### Advanced Installation Methods

**Using Extensiphy is limited to Linux at the moment.** Using Ubuntu will ensure
the smoothest performance. If you want to use another distro,
you'll have to make sure you install analogous one-liners and all that.
You have been warned.


#### Requirements

You can forgo installing dependencies with Conda or Docker and
instead install everything by hand if you feel comfortable with computer pathing.

Dependencies (Separate programs you'll need to install):

1. [Python 3](https://www.python.org/)
2. [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
3. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)
4. [Seqtk](https://github.com/lh3/seqtk)
5. [Samtools](http://www.htslib.org/)
6. [Bcftools](http://www.htslib.org/)
7. [Fastx](http://hannonlab.cshl.edu/fastx_toolkit/download.html)
8. [Dendropy](https://dendropy.org/)

Additionally, Extensiphy comes with an additional pipeline for generating a
phylogenetic tree from scratch: **Gon\_phyling**.
These programs are not required for running Extensiphy itself but Gon\_ling
can be useful if you have a lot of data and aren't interested in hand selecting
the loci/genes you include in your alignment. Gon\_phyling's dependencies are as
follows:

1. [PARSNP](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
2. [Spades](https://github.com/ablab/spades)
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (BBDUK.sh)
4. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)


### Apt-get dependency install
Many programs for running Extensiphy are available with apt-get.
Run the commands found below to install:

```bash
apt-get install raxml
apt-get install seqtk
apt-get install samtools
apt-get install bcftools
pip install dendropy
```

Installs with apt-get for Gon\_phyling are not currently available.
You will have to install these programs manually or with conda.
