#! /bin/bash
# Comprehensive test of Extensiphy options and outputs


set -e
set -u
set -o pipefail

# establishes the path to find the phycorder directory
TEST_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
EP_DIR=$(cd "$( dirname "${TEST_DIR/..}" )" && pwd)

test_help_menu () {

  echo "${TEST_DIR}"
  echo "${EP_DIR}"

  ep_menu_test=$(${EP_DIR}/extensiphy.sh -h)
  ep_menu_expected_results='Correct version of bfctools found.
Correct version of samtools found.
seqtk found
bwa-mem2 found
fastx toolkit found
vcfutils.pl found


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
 (-n) Set size of loci size cutoff used as input or output (Options: int number)(DEFAULT: 700)'

 echo ${ep_menu_test}
 echo ${ep_menu_expected_results}

 if [[ ${ep_menu_test} != ${ep_menu_expected_results} ]]
 then
   echo "PROBLEM WITH HELP MENU"
   exit
 fi

}

test_basic_run () {



}

test_help_menu
