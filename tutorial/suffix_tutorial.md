---

title: Suffix Tutorial

author: Jasper Toscani Field

output: html_document

---
# File Suffix Tutorial

File suffixes are very important in bioinformatics in general and when using Extensiphy in particular.
This tutorial is here to help you understand how to think about and use suffixes.

## What is a suffix?

A suffix is the end of a file name. If you have a file with the name

```
random_file.txt
```

then the

```
.txt
```
part is the suffix.  

The beginning of the file name is called the prefix:
```
random_file
```

The suffix could be anything like `.doc`, `.xls`, `.csv`, `.fastq`, etc.


## Paired file suffixes

While a suffix is the end of a file, it can be more than just the `.fastq` at the end of the file.  

Lets look at an example  

This FASTA formatted file has a pretty simple suffix

```
example_sequence.fasta
```

Lets take a look at some different files you might encounter during your bioinformatics journey.
These paired-end FASTQ formatted files also have suffixes. What do you think they are?

```
taxon_name_R1.fastq
taxon_name_R2.fastq
```

Ok, its definitely correct to say that both of these files share the same suffix: `.fastq`.  
Heres the tricky part. Because these are paired files, the suffix for each individual file could also be:

```
_R1.fastq
_R2.fastq
```

Yeah, its a little confusing at first. Lets look at these files some more.  
These files share a name:

```
taxon_name
```

This shared name could refer to a taxon name, a gene name or a sequencing run.  
The important thing to remember is that programs that use paired files and request that you specify a suffix
ultimately want you to slice the file names into shared components (in this case `taxon_name`) and separate file suffixes (in this case `_R1.fastq` and `_R2.fastq`).  


## An Example

Lets pose a hypothetical situation based on a common bioinformatics job.  
You want to align some raw reads to a sequence to find new variant nucleotides that belong to the organism that you just sequenced.  
Maybe you have an existing genome belonging to a chicken and you want to align another bird's sequences to that genome. Standard bioinformatics stuff.  
Your hypothetical command is below.

```
our_hypothetical_program --input_alignment file_name.fasta --input_read_1 taxon_name_R1.fastq --input_read_2 taxon_name_R2.fastq --suffix_1 _R2.fastq --suffix_2 _R2.fastq
```
Ok, thats a lot of information, lets break it down.
* `out_hypothetical_program` is the name of the program we're running. We call this with a specific goal in mind.  
* `--input_alignment` tells the program that the next word in the command is the input sequence alignment for our program. This is a program specific option and `file_name.fasta` is the file that this option points to.  
* `sequence_name` is the prefix of some of our input files and presumably, an important identifier like a taxon or sequence run name.  
* `input_reads_1` and `input_reads_2` point to the read files using their complete names. This is important because your telling the computer where to find the input data.
* `suffix_1` and `suffix_2` are exactly what we've been describing. The suffixes that you want to identify in your file names.  

Here's the important part of the interaction between the complete read file names and the suffixes. The program is going to subtract the suffixes from your raw read file names, identifying that file identifier is `taxon_name`.  

You'll output a file with a name like
```
taxon_name.bam
```

This file will contain all the reads that align to the genome. Pretty cool huh?


## It can be a little more complicated

Now, it can be even trickier than that.  

Our files are:
```
raw_reads_R1.fastq
raw_reads_R2.fastq
```

And we described these portions of the files as our suffixes:
```
_R1.fastq
_R2.fastq
```

There is some nuance in what the suffixes are.
We could use the taxon name:
```
raw_reads_
```

Notice that we kept the underscore `_` at the end of the name.  
In this scenario, the suffixes would now be:
```
R1.fastq
R2.fastq
```

You have a lot of control over describing how the file name is split into prefix and suffix.
You'll have to make sure that some sort of identifier between the two separate files is kept in the suffix, in this case `R1` and `R2`.  

This is the end of the prefix-suffix tutorial. Hopefully this has helped you understand how to work with suffixes when a program requires them.
