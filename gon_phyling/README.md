# Gon Phyling

Extensiphy comes with an pipeline for generating an alignment and a
phylogenetic tree from scratch: **Gon\_phyling**.
This program is not required for running Extensiphy itself but Gon\_phyling
can be useful if you have a lot of data and aren't interested in hand selecting
the loci/genes you include in your alignment. Gon\_phyling is more for internal use
but if you've dug into everything Extensiphy has to offer, you might find it useful.

## Installation and Use
Gon\_phyling also relies on dependency programs.
Below are some instructions for installing Gon\_phyling using Docker.
Dependency programs are found at the end of the tutorial if you'd rather install each program by hand.

### Building and testing your own Extensiphy Docker image
First we'll building the Docker image and a container to test your Gon\_phyling installation. Then we'll connect your data to a new container so you can begin building your own alignments!

1. Make sure you have [Docker installed](https://www.docker.com/products/docker-desktop) according to your operating system.

2. To pull the Docker build of Gon\_phyling, run this command.

```bash
docker pull mctavishlab/gon_phyling
```

3. We'll build your Gon\_phyling Docker container using this command.
* `-i` makes the container interactive.
* `-t` specifies the image to use as a template.
* `--name` specifies the container name.
```bash
docker run --name gon_phy_container -i -t mctavishlab/gon_phyling bash
```
Your command line prompt should change to indicate that you are now working
inside your Gon\_phyling container.

You can exit the docker container by typing `exit`.

andTo restart it and return to interactive analyses, run:

```bash
docker container restart gon_phy_container
docker exec -it gon_phy_container bash
```


### Quick test run
If you have followed the install approaches above, you are now ready to try a test run!  
We'll use the raw read files in folder:
```
~/extensiphy/gon_phyling/gon_phy_testdata
```
You'll see that there are a number of files with a `.fastq` suffix (for information on suffixes, check out this tutorial (LINK)). If we examine those files, we'll see that they all contain data in the same format like this:

```
@taxonxon_66_R1_sim_taxonxon_109_R1-2040/1
CGCGTGGCAGGAAACCAGCCACGCCATCCAACGCCTGCGCGACAACCCTGCCTGCGCCGACAGCGAGCTCGCCCTGATTGGCGACAACGAACGCAGCGCA
+
;C0FFFDFHHHGAAEJJJIIJJJDJDJJ4JHBIHJ(J*IIIC)JFCICJBIJJH1JIJJIIDJDIHCFH'CD?BBBCFEC7DB((BFD?DDD9C@D?DDE
@taxonxon_66_R1_sim_taxonxon_109_R1-2038/1
ACGGTAGGATGAATGTGAATGAGCCGACCGAGCATCGGCATCTGTCGCGCTCTGCTGTGTCGCACCATTTAAAAATCATGCTGCAAGCCGGAGCGGTGGC
+
BCCDFFDDHHHHHIGIJGJJII@JJIF<GJIJID4G'EII?2FJHIFGJE@HJDDIJGI=FD'8F0DEC0<D9DCC,DIH)=DECDD@ACBDDBDBCADD
```

Now, from the starting folder of your Docker container, run:

```bash
./gon_phyling/gon_phyling.sh -d ./gon_phyling/gon_phy_testdata -1 _R1.fastq -2 _R2.fastq
```

This simple run uses a couple of flags (You might be familiar with some of them if you've used Extensiphy already).
* The `-d` flag provides the path to your directory of fastq files.
* The `-1` and `-2` flags specify the filename endings for each of the readfiles. (defaults are `_R1.fq` and `_R2.fq` )

Once Gon\_phyling has finished running on the test data, you should see a lines saying:
```
Newly constructed alignment file can be found here: gon_phy_testdata/outputs/RAxML_bestTree.core_genome_run.out
Newly estimated phylogeny file can be found here: gon_phy_testdata/outputs/combo.fas

Genome assembly, alignment construction and phylogenetic estimation complete.

```
* The tree file will look like the one found in the Extensiphy tutorial but with different taxa names ![Extensiphy tutorial](tutorial/images/tree_image_2.png?raw=true).

* The alignment file will look like the one found in the Extensiphy tutorial but with different taxa names ![alignment image](tutorial/images/starting_alignment.png?raw=true).  

* If you are using docker - exit the container by typing
```
exit
```

* You can copy the alignment to your local directory using:

```
docker cp gon_phy_container:/project/extensiphy/gon_phyling/gon_phy_testdata/outputs/combo.fas .
```

* To get right down to business and build your own alignment, continue to the next section.


### Using Gon\_phyling on your own data.

We'll use brackets `[]` to indicate variables you should replace with your own files or paths.
Replace the `[stuff inside the brackets]` with the appropriate paths and folder names you've used so far.

If you have installed Gon_phyling locally, you can just pass in the paths to your data, and run the analysis.

````
./gon_phyling/gon_phyling.sh -d [path to your_directory_of_reads] -1 [r1_suffix] -2 [r2_suffix]
````

If you are using docker, it is simplest to link your data directory to a new container.

Put the input alignment and raw reads you want to align in a directory. e.g. ``[my_data_dir]``

We'll build a new Gon\_phyling Docker container and connect the directory containing your data to the container.

```bash
docker run --name gon_phy_container_link -i -t -v [/path/to/my_data_dir]:/project/linked_data mctavishlab/gon_phyling bash
```

This shares the `my_data_dir` folder between your operating system and the docker container. (In this example it is named `my_data_dir` locally and `linked_data` in your docker container, but you can name them the same thing in both places if you prefer.)

Now you can run Gon\_phyling as earlier but we'll specify the directory where your data is located.

```bash
./gon_phyling/gon_phyling.sh -d /project/linked_data -1 [suffix_1] -2 [suffix_2]
```

By putting the outputs into the linked directory, you can access them directly through your operating system without having to copy them.

* The output files can be found in the `outputs` directory that is created in the input directory.

## Program Options
These are the various flags you can use with Gon\_phyling.

### Input Options:
```
- (-d) directory of paired end reads. All output folders and files will be contained here
- (-g) the name of the genome you wish to use as a reference during loci selection (if any)(DEFAULT: NONE)
- (-1, -2) suffixes of paired-end input files in read directory (DEFAULT: -1 R1.fastq -2 R2.fastq)
```

### Output Options
```
- (-b) bootstrapping setting. Do you want to perform 100 boostrap replicates and add the support values to the best tree? (DEFAULT: OFF)
- (-o) output type. Output either a concatenated multiple sequence alignment only or also output separate loci alignment files (DEFAULT: LOCI) (OPTIONS: LOCI, LOCUS)
- (-l) Locus position file. Use if selecting -o LOCUS. Outputs a csv file tracking the loci names and their positions within the concatenated MSA (DEFAULT: gon_phy_locus_positions.csv)
```

### Performance Options
```
- (-r) gon_phyling runs. This is the number of genomes assembled at a single time (DEFAULT: 2)
- (-c) Threads for each gon_phyling run. Figure out how many cores you have available and input [# of threads x # of parrallel genome assemblies] = cores you can allocate. (DEFAULT: 2)
```

## Starting from Raw Reads
Creating a starting alignment!
You need a alignment with any number of taxa in order to update with new taxa.
The commands below will use Gon_phyling to assemble a starting alignment that can then be built-upon with Extensiphy.

1. Move some fraction or subset of your read files to a new directory for assembly and starting alignment construction.
2. run from the `extensiphy directory`, replacing the paths and suffix information:
```bash
./gon_phyling/gon_phyling.sh -d [PATH/TO/NEW/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2]
```
3. Use the produced alignment file and the rest of the reads as the inputs for a full Extensiphy run by running:
```bash
 ./multi_map.sh -a [PATH/TO/ALIGNMENT/FILE] -d [PATH/TO/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2].
```
## Dependencies
Gon\_phyling's dependencies are as follows:

1. [PARSNP](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
2. [Spades](https://github.com/ablab/spades)
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (**BBDUK.sh and repair.sh are the programs you need from this package**)
4. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)
