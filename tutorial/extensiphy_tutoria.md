---

title: Extensiphy Tutorial

author: Jasper Toscani Field

output: html_document

---
# Extensiphy tutorial

Hello!

This tutorial will walk you through installing and using **Extensiphy**, a program for rapidly updating multiple sequence alignments with new sequences. I try not to assume too much knowledge of programming and program use standards so this tutorial should be accessible for beginning and advanced bioinformaticians.


## Description

Extensiphy is a tool for updating an existing multiple sequence alignment. Say you built a sequence alignment and phylogeny for a group of bacteria during an outbreak and then received some new sequencing data that you wish to quickly incorporate into the phylogeny. Extensiphy makes it convenient to do this while also ensuring you can easily do this again any time you acquire new data.

### Use cases and required files

You can use Extensiphy in multiple ways. Depending on the input and output you wish to receive, you will need different data in different formats. Here are some example use cases:

#### You want an alignment with new taxa added and a tree from that alignment
You need:
* Alignment file (fasta format)
* Directory of paired-end read files
* (Optional) a phylogeny generated from the input alignment

#### You want multiple single locus files and a phylogeny from all of those alignment files
You need:
* Multiple seperate locus multiple sequence alignment files (fasta format) **OR**
  * A single concatenated multiple sequence alignment file
  * A CSV file illustrating locus lengths and positions (formatting explained later)
* Directory of paired-end read files
* (Optional) a tree generated from combining all the loci


### A few notes on code examples and using this tutorial

Using Extensiphy requires you to use the command line. If you are new to the command line, we've included some additional tutorials to help you get to grips with the new environment you'll have to navigate. Check out the links below for these tutorials.

* [Intro to Command Line Tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/command_line_tutorial.md)

* [Intro to File Suffixes Tutorial](https://github.com/McTavishLab/extensiphy/blob/main/tutorial/suffix_tutorial.md)

* You'll also see `$` (dollar signs) in front of code snippets. These should not be included when you copy and paste the code into your terminal window. These dollar signs are only included to indicate the line with the actual code snippet and not the results of the command.


## Installation and testing

To install Extensiphy and its dependency programs, you should use one of the [methods outlined in the README](https://github.com/McTavishLab/extensiphy#setup-and-use). Once you've followed the installation instructions of your choice, come back here and continue to running Extensiphy!


## Running Extensiphy

Ok, you've installed Extensiphy and all of its dependency programs.
Regardless of the method of installation you used, we're going to use the example dataset in the `testdata` directory to explore some of Extensiphy's functionality. All of the commands will work regardless of the installation method.
Lets get to it!

### Extensiphy Help Menu
Extensiphy takes command line arguments to update a sequence alignment with new taxa sequences. Lets look at the options used by Extensiphy. Extensiphy use revolves around calling the

```bash
$ ./extensiphy.sh
```

command followed by flags (dashes next to a letter corresponding the option you wish to use or input).
Calling the help menu with the following command:

```bash
$ ./extensiphy.sh -h
```

has the following output:

```
Extensiphy is a program for quickly adding genomic sequence data to multiple sequence alignments and phylogenies.
     View the README for more specific information.
     Inputs are generally a multiple sequence file in fasta format and a directory of
     Fastq paired-end read sequences.     


 EXAMPLE COMMAND:     

 /path/to/extensiphy.sh -u ALIGN -a /path/to/alignment_file -d /path/to/directory_of_reads [any other options]     

 REQUIRED FLAGS     
 (-a) alignment in fasta format,     
 (-d) directory of paired end fastq read files for all query taxa,     
 (-u) produce only an updated alignment or perform full phylogenetic estimation (ALIGN or PHYLO) (DEFAULT: ALIGN)


 OPTIONAL FLAGS     
 (-t) tree in Newick format produced from the input alignment that you wish to update with new sequences or specify NONE to perform new inference (DEFAULT: NONE),     
 (-m) alignment type (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA) (DEFAULT: CONCAT_MSA),     
 (-o) directory name to hold results (DEFAULT: creates EP_output),     
 (-i) clean up intermediate output files to save HD space (Options: CLEAN, KEEP)(DEFAULT: KEEP),     
 (-r) Selected a reference sequence from the alignment file for read mapping or leave as default and a random reference will be chosen (DEFAULT: RANDOM),     
 (-p) number of taxa to process in parallel,     
 (-c) number of threads per taxon being processed,     
 (-e) set read-type as single end (SE) or pair-end (PE) (DEFAULT: PE)     
 (-1, -2) suffix (ex: R1.fastq or R2.fastq) for both sets of paired end files (DEFAULTS: R1.fq and R2.fq),     
 (-g) output format (CONCAT_MSA or SINGLE_LOCUS_FILES) (DEFAULT: CONCAT_MSA),     
 (-s) specify the suffix (.fa, .fasta, etc) (DEFAULT: .fasta),     
 (-b) bootstrapping tree ON or OFF (DEFAULT: OFF)     


 if using single locus MSA files as input,     
 (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv),     
 (-n) Set size of loci size cutoff used as input or output (Options: int number)(DEFAULT: 700)
 ```

Extensiphy has a number of default settings for these so you will not always have to explicitly use all of these options for every run. The use of these flags depends on the input you wish to use and the output you desire to have at the end of a run.


### Examining our input data
First lets try a basic test case. Within the Extensiphy folder is a folder called

```
testdata
```

This folder contains a variety of files for testing your installation and use of Extensiphy.

Lets take a look at some of the files in the testdata folder.

```bash
$ grep ">" ./testdata/combo.fas
```
This command should output:

```
>taxon_12
>taxon_16.ref
>taxon_22
>taxon_25
>taxon_10
>taxon_17
>taxon_20
>taxon_23
>taxon_21
>taxon_18
>taxon_13
>taxon_15
>taxon_14
>taxon_27
>taxon_11
>taxon_1
>taxon_28
>taxon_26
>taxon_19
>taxon_24
```

A quick search of our `combo.fas` sequence alignment shows that we have a number of sequences, presumably in fasta format because the sequence names start with a `>`. Lets count our sequences so we can later assess if our Extensiphy runs are actually adding sequences to the alignment.

```bash
$ grep -c ">" ./testdata/combo.fas
```
You should see there are this many sequences in the sequence alignment.

```
20
```

We can also examine some of the sequences directly. This process is easier with a program like [seaview](http://doua.prabi.fr/software/seaview) but if you want to skip installing another program, you can run the following command to look at the first 4 lines of the alignment file.

```bash
$ head -4 ./testdata/combo.fas
```
The output of this command will be rather intense. **DONT PANIC**.
You should see a few header lines starting with a `>` and the name of the taxon followed by another long line of DNA sequence.
There should be 2 sequences total output by this command. This displays that there is indeed aligned DNA sequences in our alignment file.

Finally, lets examine the size of our alignment file.
This will give us tool to show that our alignment file is expanding with new sequences as we run Extensiphy.

```bash
$ ls -lh ./testdata/combo.fas
```

This command should output something similar to the following information.

```
-rw-rw-r-- 1 your_user_name your_user_name 200K Oct  4 11:43 ./testdata/combo.fas
```

This tells us that the alignment file is 200 kilobytes in size.
Ok, we've examined the alignment file, lets take a look at the fastq files we'll be working with.
Run the following command to show all of the read files we'll use in this tutorial.

```bash
$ ls ./testdata/*.fq
```

Our output should look like this:

```
./testdata/taxon_30_R1.fq  ./testdata/taxon_31_R1.fq  ./testdata/taxon_32_R1.fq
./testdata/taxon_30_R2.fq  ./testdata/taxon_31_R2.fq  ./testdata/taxon_32_R2.fq
```

Examining these results, we can see that our files are paired-end fastq format files, denoted by the suffixes (`_R1.fq` and `_R2.fq`).
We have read files for 3 taxa.
Lets take a look inside to confirm that the read files look like real fastq files.
Run the following command.

```bash
$ for i in $(ls ./testdata/*.fq ); do head $i; done
```

The output of this command should look like repetitions of the following:

```
@taxonxon_66_R1_sim_taxonxon_127_R1-2040/1
AGACTCCAGCCTGCCTGACGTGGTCAGCCACGGCGAAGGCCGCGCCGACTTCGCGCTTCACGGCGGCAGTATTTCCGCCGACTTGGGCGTCGCGCTGCAA
+
CC?FDAFAGH2HHJIJJ?@GE@GFJJHJ=<IJIIHJJE.?J8IGJ*FI>GJEH@HDJJ>J.HHHB?F(DF=CFFECE7@FD.A8@DC9DD5D=DDCE:CC
@taxonxon_66_R1_sim_taxonxon_127_R1-2038/1
ATGACCGTCGCCGGTCGCCGTGCCGACTACGGCAAACGGGCAGCGTTCGCGTTCGCAGAGGGCGCGGAAGGTGTACAAATCTTTTTCCAAAATCGACAAG
+
@;@8FF?FFGHGFFIHDJI-IGGJJ?IJHHHJIJJ?FJJICJJH*HJBHEC9IHJFHHE((JIFGEGJDEDECD:A0&DDCEDFCDBF>DE>DDD?@C<D
@taxonxon_66_R1_sim_taxonxon_127_R1-2036/1
```

Bioinformatics data is often formatted in distressing ways. Again, DONT PANIC.
This is what we want to see. Its not critical that you recognize all of the lines of a fastq file.
A cursory glance will show you lines that have the taxon's name and lines that have the DNA sequence information.
Great! With all of our files containing short reads, we're ready to run Extensiphy and update our alignment with these new taxa.


### A Basic Extensiphy Run
Lets use Extensiphy to update a multiple sequence alignment with some new data. Run this command in the Extensiphy folder:

```bash

$ ./extensiphy.sh -a ./testdata/combo.fas -d ./testdata -o first_extensiphy_run

```

This command takes in an alignment file (combo.fas) with the -a flag and a directory containing some paired-end read files with the -d flag. The -o flag specifies the name of our output folder. The other flags use their default values in this case. The result of this is that Extensiphy will add the taxa sequences from the read files to combo.fas. You can examine the results by looking at the updated alignment file:

```
~/first_extensiphy_run/outputs/extended.aln
```

Lets take a look at our new alignment file that SHOULD have 3 new taxa/sequences added to it.

```bash
$ grep ">" ./first_extensiphy_run/outputs/extended.aln
```

The output you should see is:

```
>taxon_30_
>taxon_31_
>taxon_32_
>taxon_12
>taxon_16.ref
>taxon_22
>taxon_25
>taxon_10
>taxon_17
>taxon_20
>taxon_23
>taxon_21
>taxon_18
>taxon_13
>taxon_15
>taxon_14
>taxon_27
>taxon_11
>taxon_1
>taxon_28
>taxon_26
>taxon_19
>taxon_24
```

Follow that up with a quick count of the number of sequences in the alignment file.

```bash
$ grep -c ">" ./first_extensiphy_run/outputs/extended.aln
```

This command should return
```
23
```

Awesome! You've successfully added 3 new sequences to a sequence alignment without having to realign the entire file.
Lets play with some of the optional functionality Extensiphy provides.

### Building a phylogeny

Extensiphy can also build a phylogeny for you based on the updated sequence alignment using RAxML.
Adding the `-u PHYLO` option to the command will instruct Extensiphy to estimate a phylogeny once the sequence alignment has been updated.
Here is our new command.

```bash

$ ./extensiphy.sh -a ./testdata/combo.fas -d ./testdata -u PHYLO -o second_extensiphy_run

```

The outputs will all appear in the `outputs` folder. When looking at your `outputs` folder, you'll see all the RAxML output files.

```
extended.aln          RAxML_bestTree.consensusFULL  RAxML_log.consensusFULL            RAxML_result.consensusFULL
extended.aln.reduced  RAxML_info.consensusFULL      RAxML_parsimonyTree.consensusFULL
```

The tree file, `RAxML_bestTree.consensusFULL` is the maximum likelihood most likely phylogeny based on the input data.

### Updating an existing phylogeny

In some cases, you might already have a phylogeny and you want to add your new sequences to the alignment and phylogeny in one command.
Run the following command:

```bash
$ ./extensiphy.sh -u PHYLO -a ./testdata/combo.fas -t ./testdata/combo.tre -d ./testdata -o third_extensiphy_run
```
All of the output files will be the same as generating a phylogeny from scratch.
The -t flag indicates that you're assigning a tree file as input that corresponds with the alignment file you indicated. The tree file is then used as a starting tree when performing the new, full maximum likelihood search instead of a de novo generated tree.

Lets look at the original tree that DOESN'T included our new sequences: ![sequences](images/tree_image_3.png?raw=true)

Now, if either run completed successfully you'll see a full phylogenetic tree that looks like this: ![this](images/tree_image_2.png?raw=true)

We just added 3 new taxa to a starting multiple sequence alignment and obtained a tree that includes these new taxa. Notice that the new sequences we wanted to add (taxon_30, taxon_31 and taxon_32) have been added to the clade highlighted in the red box.

### Bootstrapping an updated phylogeny

By using the `-b` flag, we can toggle bootstrapping on or off when estimating a phylogeny.
You'll need to specify `-b ON` to turn on bootstrapping. Bootstrapping is off by default.
Run this command to run Extensiphy on our test data with bootstrapping on for our phylogenetic estimation.

```bash
$ ./extensiphy.sh -u PHYLO -b ON -a ./testdata/combo.fas -t ./testdata/combo.tre -d ./testdata -o fourth_extensiphy_run
```
One the run completes, you should see a few additional files compared to our previous outputs.
The complete `outputs` folder should look something like this:

```
extended.aln                                                      RAxML_info.consensusFULL
extended.aln.reduced                                              RAxML_info.consensusFULL_bootstrap
RAxML_bestTree.consensusFULL                                      RAxML_info.majority_rule_bootstrap_consensus
RAxML_bipartitionsBranchLabels.majority_rule_bootstrap_consensus  RAxML_log.consensusFULL
RAxML_bipartitions.majority_rule_bootstrap_consensus              RAxML_result.consensusFULL
RAxML_bootstrap.consensusFULL_bootstrap
```

These are standard [RAxML](https://github.com/stamatak/standard-RAxML) outputs.
Visit the linked github repo and brush up on what each file contians.


### Choosing a specific reference sequence

Selecting a particular reference from the alignment may be important to a particular analysis. You can select a reference by using the `-r` flag followed by the name of the taxon/sequence. First, lets look at a list of our included taxa. Run the following command in the Extensiphy directory:

```bash
$ grep ">" ./testdata/combo.fas
```

You should see something like this:

```
>taxon_12
>taxon_16.ref
>taxon_22
>taxon_25
>taxon_10
>taxon_17
>taxon_20
>taxon_23
>taxon_21
>taxon_18
>taxon_13
>taxon_15
>taxon_14
>taxon_27
>taxon_11
>taxon_1
>taxon_28
>taxon_26
>taxon_19
>taxon_24
```

Lets take one of our taxa and use that sequence as the reference. I chose taxon_11 because its on a long branch far from where we expect our sequences to placed in the tree. Make sure you leave off the carrot (>) from the taxon name. Run this command:

```bash
$ ./extensiphy.sh -a ./testdata/combo.fas -u PHYLO -d ./testdata -t ./testdata/combo.tre -o fifth_Extensiphy_run -r taxon_11
```

Did you get the same tree as the output of our original run?



### Starting with Multiple Single Locus Alignments

Extensiphy will accept multiple, single locus files as input instead of a single, concatenated sequence file.
This is useful if you would also like to output single locus files that have been updated with new sequences for each taxon.
Lets try starting with multiple single locus files instead of a single, concatenated sequence file. Use the following commands to look at our single locus files:

```bash
$ ls ./testdata/single_locus_align_dir/
```
This command should produce this result:
```
/testdata/single_locus_align_dir/single_locus_1_.fasta
/testdata/single_locus_align_dir/single_locus_2_.fasta
```

Both of these files contain homologous sequences for a number of taxa. The first one is a rather big sequence of over 10,000 bases. The second one is a smaller sequence of 4 bases. Take a look at them like so:

```bash
$ head -2 /testdata/single_locus_align_dir/single_locus_1_.fasta
$ head -2 /testdata/single_locus_align_dir/single_locus_2_.fasta
```

You'll see that one file's sequence is indeed very large while the second file's sequence is only a few letters.
This is deliberate to display a function of Extensiphy when selecting which data to input.
Extensiphy should only be used with loci over 700 bases long.
The read alignment software used in Extensiphy produces the best results when the reads are aligned to longer sequences.
We've set this default length to 700 bases but we'll show you how to adjust it in a moment.
If taking these files as input, the smaller sequence file will be identified and removed before construction of a concatenated sequence.
Lets run a new analyses.

Now, enter the following command:

```bash
$ ./extensiphy.sh -a ./testdata/single_locus_align_dir -d ./testdata -m SINGLE_LOCUS_FILES -o sixth_Extensiphy_run
```

The `-m` flag allows you to specify a directory with any number of input options. by using -m SINGLE_LOCUS_FILES, we are indicating that the alignment option (-a) will point to a directory containing multiple single locus alignment files that share all the sample taxon names. It is VERY important that all the taxa labels have the same names or this function will not work. This run will take the single loci MSA files, check for loci longer than the cut-off of 700 nucleotides and construct a concatenated alignment of those loci. A file capturing the length and positions of those loci can be found in the file:

```
loci_positions.csv
```

This file is a comma delimited file capturing the loci's position in the concatenated alignment, the loci's file name and the loci's length. This will be useful if you decide to split your concatenated multiple sequence alignment back into single locus (gene) alignment files or just want to know how long each locus in your concatenated alignment is.


### Updating a species tree produced from single locus files

If you have already produced a species-tree from multiple single locus alignments, Extensiphy can take that tree as input along with your separate locus files. Run the following command to read in some multiple single locus alignments and the tree corresponding to the relationships inferred from all of the combined loci. We'll also give our new Extensiphy run output folder a new name so we can distinguish it from our old run.

```bash
$ ./extensiphy.sh -u PHYLO -a ./testdata/single_locus_align_dir -d ./testdata -t ./testdata/combo.tre -m SINGLE_LOCUS_FILES -o split_Extensiphy_run

```

### Outputting updated single locus alignment files

Maybe you also want to output single locus alignment files that have been updated with your new query sequences. Run this command to do that:

```bash
$ ./extensiphy.sh -u PHYLO -a ./testdata/single_locus_align_dir -d ./testdata -t ./testdata/combo.tre -m SINGLE_LOCUS_FILES -g SINGLE_LOCUS_FILES -o locus_out_Extensiphy_run

```

The -g flag allows you to specify the format you wish your output alignments to take (CONCAT_MSA or SINGLE_LOCUS_FILES). When selecting SINGLE_LOCUS_FILES, a concatenated alignment file is produced AS WELL as a folder that contains the separated single locus alignment files in fasta format. Lets take a look at them.

From the Extensiphy folder, run:

```bash
$ ls ./locus_out_Extensiphy_run/outputs/updated_single_loci/

single_locus_1_.fasta
```

You'll notice that there is only a single file here corresponding to a single locus alignment. This is because one of the loci we tried to input into Extensiphy was only 4 nucleotides long, far short of the 700 nucleotide cutoff for using Extensiphy!


This concludes the tutorial. Hopefully you understand a little more about using Extensiphy and how to apply Extensiphy to your use-case and data.
