---

title: RapUp Tutorial

author: Jasper Toscani Field

output: html_document

---
# RapUp tutorial

Hello!

This tutorial will walk you through installing and using **RapUp**, a program for rapidly updating phylogenies and multiple sequence alignments.

## Description

RapUp is a tool for updating an existing multiple sequence alignment and a phylogeny inferred from that alignment. Say you built a phylogeny for a group of bacteria during an outbreak and then received some new sequencing data that you wish to quickly incorporate into the phylogeny. RapUp makes it convenient to do this while also ensuring you can easily do this again any time you acquire new data.

### Use cases and required files

You can use RapUp in multiple ways. Depending on the input and output you wish to receive, you will need different data in different formats.Here are some example use cases:

##### You want an alignment with new taxa added and a tree from that alignment
You need:
* Alignment file (fasta format)
* Directory of paired-end read files
* (Optional) a tree generated from input alignment

##### You want multiple single locus files and a tree from all of those alignment files
You need:
* Multiple seperate locus multiple sequence alignment files (fasta format) **OR**
  * A single concatenated multiple sequence alignment file 
  * A CSV file illustrating locus lengths and positions
* Directory of paired-end read files
* (Optional) a tree generated from combining all the loci


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

Additionally, RapUp comes with an additional pipeline for generating a phylogenetic tree from scratch: **Gon\_phyling**. These programs are not required for running RapUp itself but Gon\ling can be useful if you have a lot of data and aren't interested in hand selecting the loci/genes you include in your alignment. Gon\_phyling's dependencies are as follows:

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
conda create -n rapup samtools hisat2 seqtk bcftools fastx-toolkit dendropy raxml
```

Activate your installation

```bash
conda activate rapup
```

Conda install recipe on the way.



### Pathing

RapUp will need to automatically look for these programs in your computers PATH. If you're new to the inner workings of computers, think of your PATH as a set of programs or locations on your computer that your computer automatically knows the location of. The following is a basic tutorial on adding programs to your PATH.

First you'll need to find the ```.bash_profile``` file on your computer. This file lives in your home directory when you first open a terminal. Here's how you find it:

```bash
open terminal
$ls -a
.bash_profile
other_files_in_your_home_directory

$echo $PATH
home/your_name/bin:home/your_name/other_folders/in_your/path:

$vim .bash_profile
```
At this point you'll need to add the programs you've downloaded and unpacked to your .bash_profile in the following format. Add all of this as one line, filling in the appropriate absolute pathing. If you don't see this file when using `ls -a` then you need to create it.

```bash

export PATH="/home/your_name/path/to/program:$PATH"

```

An example of what I see is:

```bash

export PATH="/home/jasper/src/SPAdes-3.13.0-Linux:$PATH"

```

Once you've added all your programs to your PATH, close the terminal window and reopen it. Then type use the `which` command to see if your computer knows where the program is located

```bash

$which hisat2

/home/your_name/hisat2_folder/hisat2

```

I installed hisat2 in a bin folder so I see:

```bash

/home/jasper/bin/hisat2

```



## Running RapUp

### RapUp Help Menu
RapUp (on branch overhaul_dev) takes command line arguments to update a phylogenetic tree with new taxa sequences. Lets look at the options used by RapUp. RapUp use revolves around calling the 

```bash
multi_map.sh
```

command followed by flags (dashes next to a letter corresponding the command you wish to use or input).

```bash
multi_map.sh -h
```

has the following output:

```bash
RapUp is a program for quickly adding genomic sequence data to multiple sequence alignments and phylogenies. View the README for more specific information. Inputs are generally a multiple sequence file in .fasta format and a directory of .fastq paired-end read sequences.


 EXAMPLE COMMAND:

 /path/to/multi_map.sh -a /path/to/alignment_file -d /path/to/directory_of_reads [any other options]

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
 (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv)
```

RapUp has a number of default settings for these so you will not always have to explicitly use all of these options for every run. The use of these flags depends on the input you wish to use and the output you desire to have at the end of a run.

First lets try a basic test case. Within the RapUp folder is a folder called 

```bash
testdata
```

### A Basic RapUp Run
This folder contains a variety of files for testing your installation and use of RapUp.

Lets use RapUp to update a multiple sequence alignment with some new data. Run this command in the RapUp folder:

```bash

./multi_map.sh -a ./testdata/combo.fas -d ./testdata

```

This command takes in an alignment file (combo.fas) with the -a flag and a directory containing some paired-end read files with the -d flag. The other flags use their default values in this case. The result of this is that RapUp will add the taxa sequences from the read files to combo.fas and will infer a new phylogenetic tree. You can examine the results by looking at these files:

```bash
~/rapup_run/combine_and_infer/extended.aln
~/rapup_run/combine_and_infer/RAxML_bestTree.consensusFULL
```

Primarily, RapUp should be used for adding new sequences to an alignment AND a tree produced from that alignment. Lets run a command to take as intput the same alignment and a tree that was produced from that alignment before adding new sequences. Run the following command:

First, remove the folder with the output we just generated. In the rapup folder, run:

```bash
$rm -r ./rapup_run
```

Then run:

```bash
./multi_map.sh -a ./testdata/combo.fas -t ./testdata/combo.tre -d ./testdata
```

The -t flag indicates that you're assigning a tree file as input that corresponds with the alignment file you indicated. The tree file is then used as a starting tree when performing the new, full maximum likelihood search instead of a randomly generated tree.

If either run completed successfully you'll see a full phylogenetic tree that looks like this: ![this](images/tree_image_1.png?raw=true)

We just added 3 new taxa to a starting multiple sequence alignment and obtained a tree that includes these new taxa.



### Starting with Multiple Single Locus Alignments
Lets try starting with multiple single locus files instead of a single, concatenated sequence file. Use the following commands to look at our single locus files:

```bash
$ls /testdata/single_locus_align_dir/
/testdata/single_locus_align_dir/single_locus_1_.fasta
/testdata/single_locus_align_dir/single_locus_2_.fasta
```

Both of these files contain the homologous sequence for a number of taxa. The first one is a rather big sequence of over 10,000 bases. The second one is a smaller sequence of 4 bases. Take a look at them like so:

```bash
$head -2 /testdata/single_locus_align_dir/single_locus_1_.fasta
$head -2 /testdata/single_locus_align_dir/single_locus_2_.fasta
```

You'll see that one file's sequence is indeed very large while the second file's sequence is only a few letters. This is deliberate to display a function of RapUp when selecting which data to input. RapUp should only be used with loci over 1,000 bases long. If taking these files as input, the smaller sequence file will be identified and removed before construction of a concatenated sequence. Lets run a new analyses.

First, remove the folder with the output we just generated. In the rapup folder, run:

```bash
$rm -r ./rapup_run
```

Now, enter the following command:

```bash
./multi_map.sh -a ./testdata/single_locus_align_dir -d ./testdata -m SINGLE_LOCUS_FILES
```

This run will take the single loci MSA files, check for loci longer than the cut-off of 1,000 nucleotides and construct a concatenated alignment of those loci. A file capturing the length and positions of those loci can be found in the

```bash
loci_positions.csv
```

file. This file is a comma delimited file capturing the loci's position in the concatenated alignment, the loci's file name and the loci's length. This will be useful if you decide to split your concatenated multiple sequence alignment back into single locus (gene) alignment files.

Lets do that now!
Run the following command to split your files back into individual locus alignment files that now include your newly added taxa. We'll also give our new RapUp run output folder a new name so we can distinguish it from our old run that contains only concatenated alignemnts.

```bash
./multi_map.sh -a ./testdata/single_locus_align_dir -d ./testdata -t ./testdata/combo.tre -m SINGLE_LOCUS_FILES -o new_rapup_run

```

The -m flag allows you to specify a number of alignment input options (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA). The -o flag allows you to specify the name and location of your new output folder.

```bash
./multi_map.sh -a ./testdata/single_locus_align_dir -d ./testdata -t ./testdata/combo.tre -m SINGLE_LOCUS_FILES -g SINGLE_LOCUS_FILES -o locus_out_rapup_run

```

The -g flag allows you to specify the format you wish your output alignments to take (CONCAT_MSA or SINGLE_LOCUS_FILES). When selecting SINGLE_LOCUS_FILES, a concatenated alignment file is produced ASWELL as a folder that contains the separated single locus alignment files in fasta format. Lets take a look at them.

From the RapUp folder, run:

```bash
ls ./locus_out_rapup_run/combine_and_infer/updated_single_loci/
```

You'll notice that there is only a single file here corresponding to a single locus alignment. This is because one of the loci we tried to input into RapUp was only 4 nucleotides long, far short of the 1000 nucleotide cutoff for using RapUp!

As an additional use-case, RapUp can take Parsnp output files as inputs. This time, we'll use a raw Parsnp file (parsnp.xmfa) file as our input alignment. Parsnp is an automatic homologous locus selection tool for bacterial sequences. The output is not in a standard .fasta format so there is normally some additional processing necessary. Alternatively, RapUp will convert the .xmfa file to a fasta file when you add new sequences to the alignment. Use this command:

```bash
./multi_map.sh -a ./testdata/parsnp.xmfa -d ./testdata -t ./testdata/combo.tre -m PARSNP_XMFA -o parsnp_rapup_run

```


