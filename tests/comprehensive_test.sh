#! /bin/bash
# Comprehensive test of Extensiphy options and outputs


set -e
set -u
set -o pipefail

# establishes the path to find the phycorder directory
TEST_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
EP_DIR=$(cd "$( dirname "${TEST_DIR/..}" )" && pwd)

test_help_menu () {

  # echo "${TEST_DIR}"
  # echo "${EP_DIR}"

  ep_menu_test=$(${EP_DIR}/extensiphy.sh -h)
  ep_menu_expected_results_1='Extensiphy is a program for quickly adding genomic sequence data to multiple sequence alignments and phylogenies.'

 #echo "${ep_menu_test}"
 #echo "${ep_menu_expected_results}"

 if [[ ${ep_menu_test} != *${ep_menu_expected_results_1}* ]]
 then
   echo "PROBLEM WITH HELP MENU"
   exit
 else
   echo "Help Menu Test: PASSED"
 fi

}

test_help_menu

touch logfile.txt
# Test that EP runs and produces and updated alignment
# Tests flags: -a. -d, -1, -2, -u ALIGN, -o
# Examine: extended.aln
../extensiphy.sh -u ALIGN -a ../testdata/combo.fas -d ../testdata -1 _R1.fq -2 _R2.fq -o first_extensiphy_run >> logfile.txt 2>&1

# Test that EP runs and produces and alignment and phylogeny
# Tests flags: -a, -d, -1, -2 , -u PHYLO, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -u PHYLO -a ../testdata/combo.fas -d ../testdata -1 _R1.fq -2 _R2.fq -o second_extensiphy_run >> logfile.txt 2>&1

# Test that EP updates an existing phylogeny
# Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -u PHYLO -a ../testdata/combo.fas -t ../testdata/combo.tre -d ../testdata -1 _R1.fq -2 _R2.fq -o third_extensiphy_run

# Test that EP bootstraps a phylogeny thats being updated
# # Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -b, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL, RAxML_bootstrap.consensusFULL_bootstrap, RAxML_bipartitionsBranchLabels.majority_rule_bootstrap_consensus
../extensiphy.sh -u PHYLO -b ON -a ../testdata/combo.fas -t ../testdata/combo.tre -d ../testdata -1 _R1.fq -2 _R2.fq -o fourth_extensiphy_run

# Tests selecting a specific reference taxon
# Tests flags: -a, -d, -1, -2, -u PHYLO -o, -t, -r
# Examine:  extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -a ../testdata/combo.fas -u PHYLO -d ../testdata -t ../testdata/combo.tre -1 _R1.fq -2 _R2.fq -o fifth_Extensiphy_run -r taxon_11

# Tests using single locus alignment files as inputs instead of concatenated
# Tests flags: -a, -d, -1, -2, -u ALIGN, -o, -m
# Examine:  extended.aln
../extensiphy.sh -u ALIGN -a ../testdata/single_locus_align_dir -d ../testdata -m SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o sixth_Extensiphy_run

# Tests building a tree from single locus alignment files
# Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -o -m SINGLE_LOCUS_FILES
# Examine:  extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -u PHYLO -a ../testdata/single_locus_align_dir -d ../testdata -t ../testdata/combo.tre -m SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o split_Extensiphy_run

# Tests outputting updated single locus alignments
# Tests flags: -a, -d, -t, -m SINGLE_LOCUS_FILES, -u PHYLO, -g SINGLE_LOCUS_FILES, -1, -2, -o
# Examine:  extended.aln, RAxML_bestTree.consensusFULL, /locus_out_Extensiphy_run/RESULTS/updated_single_loci/single_locus_1_.fasta
../extensiphy.sh -u PHYLO -a ../testdata/single_locus_align_dir -d ../testdata -t ../testdata/combo.tre -m SINGLE_LOCUS_FILES -g SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o locus_out_Extensiphy_run
