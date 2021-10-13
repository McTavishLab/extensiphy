---

title: gon_phyling Tutorial

author: Jasper Toscani Field

output: html_document

---
# Gon Phyling tutorial

Extensiphy comes with an pipeline for generating an alignment and a
phylogenetic tree from scratch: **Gon\_phyling**.
This program is not required for running Extensiphy itself but Gon\_phyling
can be useful if you have a lot of data and aren't interested in hand selecting
the loci/genes you include in your alignment.

## Installation
Gon\_phyling also relies on dependency programs.
Below are some instructions for installing Gon\_phyling using Docker.
Dependency programs are found at the end of the tutorial if you'd rather install each program by hand.




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
2. run:
```bash
./gon_phyling.sh -d [PATH/TO/NEW/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2]
```
3. Use the produced alignment file and the rest of the reads as the inputs for a full Extensiphy run by running:
```bash
 multi_map.sh -a [PATH/TO/ALIGNMENT/FILE] -d [PATH/TO/READ/DIRECTORY] -1 [READ SUFFIX 1] -2 [READ SUFFIX 2].
```
## Dependencies
Gon\_phyling's dependencies are as follows:

1. [PARSNP](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
2. [Spades](https://github.com/ablab/spades)
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) (**BBDUK.sh and repair.sh are the programs you need from this package**)
4. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)
