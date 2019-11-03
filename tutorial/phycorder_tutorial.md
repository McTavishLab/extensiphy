---

title: Phycorder Tutorial

author: Jasper Toscani Field

output: html_document

---
#Phycorder tutorial

Hello!

This tutorial will walk you through installing and using **Phycorder**, a program for rapidly updating phylogenies and multiple sequence alignments.

##Description

Phycorder is a tool for updating an existing multiple sequence alignment and a phylogeny inferred from that alignment. Say you built a phylogeny for a group of bacteria during an outbreak and then received some new sequencing data that you wish to quickly incorporate into the phylogeny. Phycorder makes it convenient to do this while also ensuring you can easily do this again any time you acquire new data.

####Dependencies

Unfortunately, Phycorder requires some dependencies. You know what they say about not reinventing the wheel. We'll walk through the basics of installation and adding the installed programs to your path so that Phycorder can use them.

Using Phycorder is limited to Linux at the moment. Using Ubuntu will ensure the smoothest performance. If you want to use another distro, you'll have to make sure you install analogous one-liners and all that. You have been warned.

Dependencies (Separate programs you'll need to install):

1. [Python 3](https://www.python.org/)
2. [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)
3. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)
4. [Seqtk](https://github.com/lh3/seqtk)
5. [Samtools](http://www.htslib.org/)
6. [Bcftools](http://www.htslib.org/)
7. [Fastx](http://hannonlab.cshl.edu/fastx_toolkit/download.html)
8. [Dendropy](https://dendropy.org/)

Additionally, Phycorder comes with an additional pipeline for generating a phylogenetic tree from scratch: **Gon\_phyling**. These programs are not required for running Phycorder itself but Gon\ling can be useful if you have a lot of data and aren't interested in hand selecting the loci/genes you include in your alignment. Gon\_phyling's dependencies are as follows:

1. [PARSNP](https://harvest.readthedocs.io/en/latest/content/parsnp.html)
2. [Spades](https://github.com/ablab/spades)
3. [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
4. [RAxMLHPC](https://github.com/stamatak/standard-RAxML)


####Pathing

Phycorder will need to automatically look for these programs in your computers PATH. If you're new to the inner workings of computers, think of your PATH as a set of programs or locations on your computer that your computer automatically knows the location of. The following is a basic tutorial on adding programs to your PATH.

