# Extensiphy

Extensiphy is a pipeline that assembles homologous loci by aligning reads to a reference from a multiple sequence alignment, calls consensus, adds to the existing alignment. Homologous loci may be kept concatenated or split back into individual alignments prior to phylogenetic estimation.

![Extensiphy_worlflow](./doc/images/EP_workflow_2.png)

Extensiphy takes an alignment and sets of sequencing reads from query taxa (a). Reads are aligned to a reference sequence and a consensus sequence is called (b). The new sequence is added to the alignment and the updated alignment is used to estimate a phylogeny (c).

[Setup and Use](#setup-and-use)

[Building and testing your own Extensiphy Docker image](#building-and-testing-your-own-extensiphy-docker-image)

[Extensiphy Controls and Flags](#extensiphy-controls-and-flags)

[Output Files](#output-files)

[Gon_phyling](#gon-phyling)

[Starting from Raw Reads](#starting-from-raw-reads)

[Anaconda Installation](#anaconda-installation)

[Advanced Installation Methods](#advanced-installation-methods)

## Setup and Use

### Docker
The simplest and most hassle free way to run Extensiphy is using Docker.
the **Quick Install and Run** section will review the docker installation instructions.

### Anaconda
You can also install the dependencies of Extensiphy using Anaconda. The **Anaconda Installation** section of this readme will walk through this process in more detail.

### Advanced
If you're comfortable installing programs by hand, the **Advanced Installation Methods** section is for you. Extensiphy dependencies are also found here.

### Tutorial
We recommend you run through the [tutorial](https://github.com/McTavishLab/extensiphy/blob/dev/tutorial/extensiphy_tutoria.md) for a more in-depth walkthrough of Extensiphy's features. The tutorial will walk through different installation methods and how to run Extensiphy using different data types and options. You can copy code snippets into your terminal window.

## Building and testing your own Extensiphy Docker image
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
docker run --name ep_container -i -t mctavishlab/extensiphy bash
```
Your command line prompt should change to indicate that you are now working
inside your Extensiphy container.

You can exit the docker container by typing `exit`.

To restart it and return to interactive analyses, run:

```bash
docker container restart ep_container
docker exec -it ep_container bash
```


### Quick test run
If you have followed one of the install approaches above, you are now ready to try a test run!  
We'll use the `combo.fas` alignment file as our starting alignment. `combo.fas` can be found in:
```
/extensiphy/testdata/combo.fas
```
When you examine the alignment file with a viewer like [Seaview](http://doua.prabi.fr/software/seaview), you should see something like ![this](tutorial/images/starting_alignment.png?raw=true)

Now, either from the docker container, your anaconda env, or from the directory where you installed Extensiphy, run:

```bash
./multi_map.sh -a ./testdata/combo.fas -d ./testdata -1 _R1.fq -2 _R2.fq -u PHYLO -o EP_output
```
This is a simple run on three paired end read samples, which are found in the directory `extensiphy/testdata`
* The `-a` flag provides the path to the existing alignment to update.
* The `-d` flag provides the path to your directory of fastq files.
* The `-1` and `-2` flags specify the filename endings for each of the readfiles. (defaults are `_R1.fq` and `_R2.fq` )
* The `-u` flag specfies what analysis to run. Here we are buildinga phylogeny. (default is `ALIGN`, building an alignment only.)
* The `-o` flag specifies the output directory. (default is `EP_output`)

Once Extensiphy has finished running on the test data, you should see a lines saying:
```
Alignment file is: /project/extensiphy/EP_output/outputs/extended.aln

Tree file is: /project/extensiphy/EP_output/outputs/RAxML_bestTree.consensusFULL

```
* If you did not get this message, you'll have to check output log `ep_dev_log.txt`
to learn more about the issue before proceeding.  

The alignment file should now 3 new sequences ![this](tutorial/images/updated_alignment.png?raw=true)  
and you should when you examine the phylogeny (`RAxML_bestTree.consensusFULL`), you should see ![this](tutorial/images/tree_image_2.png?raw=true)  

We just added 3 new taxa to a starting multiple sequence alignment and obtained a tree that includes these new taxa. Notice that the new sequences we wanted to add (taxon_30, taxon_31 and taxon_32) have been added to the clade highlighted in the red box.  


* If you are using docker - exit the container by typing
```
exit
```

* You can copy the extended tree to your local directory using:

```
docker cp ep_container:/project/extensiphy/EP_output/outputs/RAxML_bestTree.consensusFULL .
```

* For a deeper walk through, take a look through the [tutorial](https://github.com/McTavishLab/extensiphy/blob/dev/tutorial/extensiphy_tutoria.md).

* To get right down to business and update your own alignment, continue to the next section.


### Using Extensiphy on your own data.

We'll use brackets `[]` to indicate variables you should replace with your own files or paths.
Replace the `[stuff inside the brackets]` with the appropriate paths and folder names you've used so far.


If you have installed Extensiphy locally, you can just pass in the paths to your data, and run the analysis.

````
./multi_map.sh -a [path to your_input_alignment] -d [path to your_directory_of_reads] -1 [r1_suffix] -2 [r2_suffix] -u [either PHYLO or ALIGN, depending on if you want a phylogeny or just and alignment] -o [your_output_dir]
````

If you are using docker, it is simplest to link your data directory to a new container.

Put the input alignment and raw reads you want to align in a directory. e.g. [my_data_dir]

We'll build a new Extensiphy Docker container and connect the directory containing your data to the container.

```bash
docker run --name ep_container_link -i -t -v [/path/to/my_data_dir]:/project/linked_data mctavishlab/extensiphy bash
```

This shares the 'my_data_dir' folder between your operating system and the docker container. (In this example it is named "my_data_dir" locally and "linked_data" in your docker container, but you can name them the same thing in both places if you prefer.)

Now you can run multi_map as earlier but we'll specify the directory where your data is located.

```bash
./multi_map.sh -a /project/linked_data/[alignment_file] -d /project/linked_data -1 [suffix_1] -2 [suffix_2] -o linked_data/[output_dir_name]
```

By putting the outputs into the linked directory, you can access them directly through your operating system without having to copy them.


## Extensiphy Controls and Flags:

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

## Output Files
* Concatenated alignment file: found in your output folder ```[OUTDIR]/outputs/extended.aln```
* Phylogeny in newick file format (if you selected to output a phylogeny): found in your output folder ```[OUTDIR]/outputs/RAxML_bestTree.consensusFULL```
* Taxon specific intermediate files (if you kept intermediate files): found in your output folder ```[OUTDIR]/[TAXON_NAME]```. .sam, .bam and .vcf files can be found in here for any additional analyses.

## Gon phyling
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

## Starting from Raw Reads
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

## Anaconda Installation
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


## Advanced Installation Methods

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
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (BBDUK.sh and repair.sh are the programs you need from this package)
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
