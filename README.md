# Extensiphy

Extensiphy is a pipeline that assembles homologous loci by aligning reads to a reference from a multiple sequence alignment, calls consensus, adds to the existing alignment. Homologous loci may be kept concatenated or split back into individual alignments prior to phylogenetic estimation.

![Extensiphy_worlflow](./doc/images/EP_workflow_2.png)

Extensiphy takes an alignment and sets of sequencing reads from query taxa (a). Reads are aligned to a reference sequence and a consensus sequence is called (b). The new sequence is added to the alignment and the updated alignment is used to estimate a phylogeny (c).

[1. Setup and Use](#1-setup-and-use)

[2. Quick Install and Run](#2-quick-install-and-run)

[3. Extensiphy Controls and Flags](#3-extensiphy-controls-and-flags)

[4. Output Files](#4-output-files)

[5. Gon_phyling](#5-gon-phyling)

[6. Starting from Raw Reads](#6-starting-from-raw-reads)

[7. Anaconda Installation](#7-anaconda-installation)

[8. Advanced Installation Methods](#8-advanced-installation-methods)

## 1. Setup and Use

### Docker
The simplest and most hassle free way to run Extensiphy is using Docker.
the **Quick Install and Run** section will review the docker installation instructions.

### Anaconda
You can also install the dependencies of Extensiphy using Anaconda. The **Anaconda Installation** section of this readme will walk through this process in more detail.

### Advanced
*For advanced users of Linux* If you're comfortable installing programs by hand, the **Advanced Installation Methods** section is for you. Extensiphy dependencies are also found here.

### Tutorial
We recommend you run through the [tutorial](https://github.com/McTavishLab/extensiphy/blob/dev/tutorial/extensiphy_tutoria.md) for a more in-depth walkthrough of Extensiphy's features. The tutorial will walk through different installation methods and how to run Extensiphy using different data types and options. You can copy code snippets into your terminal window.

## 2. Quick Install and Run
### Building and testing your own Extensiphy Docker image
First we'll building the Docker image and a container to test your Extensiphy installation. Then we'll connect your data to a new container so you can begin updating your own alignments!

1. Make sure you have [Docker installed](https://www.docker.com/products/docker-desktop) according to your operating system.

2. To pull the Docker build of Extensiphy, run this command.

```bash
docker pull mctavishlab/extensiphy
```

3. We'll build your Extensiphy Docker container using this command.
* `-i` makes the container interactive.
* `-t` specifies the image to use as a template.
* `--name` specifies the container name.
```bash
docker run --name ep_container -i -t extensiphy bash
```
Your command line prompt should change to indicate that you are now working
inside your Extensiphy container.


4. Ok, Lets run Extensiphy and test our installation!  
We want to produce an updated alignment and estimate a phylogeny from that alignment.  
You'll need to input some information using flags:
* `-a` passes Extensiphy the alignment file you wish to update.
* `-d` passes the folder containing the fastq files.
* `-1` and `-2` pass the suffixes of your fastq reads (assuming paired-end files!).
* `-u PHYLO` tells Extensiphy to estimate a phylogeny from the updated alignment.
* `-o` passes the name of the output directory.  
To test your installation, run this command:
```bash
./multi_map.sh -u PHYLO -a ./testdata/combo.fas -d ./testdata -1 _R1.fq -2 _R2.fq -o ep_output
```

Once Extensiphy has finished running on the test data, you should see a line saying:
```
Alignment file is: /usr/src/app/extensiphy/ep_output/outputs/extended.aln

Tree file is: /usr/src/app/extensiphy/ep_output/outputs/RAxML_bestTree.consensusFULL
```
* If you did not get this message, you'll have to check output log `ep_dev_log.txt`
to learn more about the issue before proceeding.

* For a deeper walk through of what has actually happened, take a look through the [tutorial](https://github.com/McTavishLab/extensiphy/blob/dev/tutorial/extensiphy_tutoria.md).

* To get right down to business and update your own alignment, continue to the next section.


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
* `-v` specifies the directory your linking to the container and where in the container your linking it.
```bash
docker run --name ep_container -i -v [/path/to/new_data_dir]:/usr/src/app/linked_data -t extensiphy bash
```

7. Now you can run the same command as earlier but you'll specify the directory where your data is located and your file suffixes.
```bash
./multi_map.sh -a /usr/src/app/linked_data/[alignment_file] -d /usr/src/app/linked_data -1 [suffix_1] -2 [suffix_2] -o [output_dir_name]
```

Once the Extensiphy run is finished, you can check the find your updated alignment in:
```
/usr/src/app/extensiphy/EP_out/outputs/extended.aln
```

## 3. Extensiphy Controls and Flags:

### Required flags
```
- (-a) alignment in fasta format,
- (-d) directory of paired end fastq read files for all query taxa,
```

### Optional flags
```
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
```

## 4. Output Files
* Concatenated alignment file: found in your output folder ```[OUTDIR]/outputs/extended.aln```
* Phylogeny in newick file format (if you selected to output a phylogeny): found in your output folder ```[OUTDIR]/outputs/RAxML_bestTree.consensusFULL```
* Taxon specific intermediate files (if you kept intermediate files): found in your output folder ```[OUTDIR]/[TAXON_NAME]```. .sam, .bam and .vcf files can be found in here for any additional analyses.

## 5. Gon phyling
Additionally, Extensiphy comes with an pipeline for generating a
phylogenetic tree from scratch: **Gon\_phyling**.
This program is not required for running Extensiphy itself but Gon\_phyling
can be useful if you have a lot of data and aren't interested in hand selecting
the loci/genes you include in your alignment.

#### Input Options:
```
- (-d) directory of paired end reads. All output folders and files will be contained here
- (-g) the name of the genome you wish to use as a reference during loci selection (if any)(DEFAULT: NONE)
- (-1, -2) suffixes of paired-end input files in read directory (DEFAULT: -1 R1.fastq -2 R2.fastq)
```

#### Output Options
```
- (-b) bootstrapping setting. Do you want to perform 100 boostrap replicates and add the support values to the best tree? (DEFAULT: OFF)
- (-o) output type. Output either a concatenated multiple sequence alignment only or also output separate loci alignment files (DEFAULT: LOCI) (OPTIONS: LOCI, LOCUS)
- (-l) Locus position file. Use if selecting -o LOCUS. Outputs a csv file tracking the loci names and their positions within the concatenated MSA (DEFAULT: gon_phy_locus_positions.csv)
```

#### Performance Options
```
- (-r) gon_phyling runs. This is the number of genomes assembled at a single time (DEFAULT: 2)
- (-c) Threads for each gon_phyling run. Figure out how many cores you have available and input [# of threads x # of parrallel genome assemblies] = cores you can allocate. (DEFAULT: 2)
```

## 6. Starting from Raw Reads
Creating a starting alignment!
You need a alignment with any number of taxa in order to update with new taxa.
The commands below will use Gon_phyling to assemble a starting alignment that can then be built-upon with Extensiphy.

1. Move some fraction or subset of your read files to a new directory for assembly and starting alignment construction.
2. run:
```bash
./gon_phyling.sh -d [PATH/TO/NEW/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2]
```
3. Use the produced alignment file and the rest of the reads as the inputs for a full Extensiphy run by running:
```bash
 multi_map.sh -a [PATH/TO/ALIGNMENT/FILE] -d [PATH/TO/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2].
```

## 7. Anaconda Installation
You can install the dependencies of Extensiphy using the Anaconda package manager.
 Install instructions for Anaconda can be found [here](https://docs.conda.io/en/latest/miniconda.html).
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


## 8. Advanced Installation Methods

**Using Extensiphy is limited to Linux at the moment.** Using Ubuntu will ensure
the smoothest performance. If you want to use another distro,
you'll have to make sure you install analogous one-liners and all that.
You have been warned.


### Requirements

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

**Gon\_phyling**.
Gon\_phyling's dependencies are as follows:

1. [PARSNP](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
2. [Spades](https://github.com/ablab/spades)
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (BBDUK.sh is the program you need from this package)
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
