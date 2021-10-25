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
../extensiphy.sh -a ../testdata/combo.fas -d ../testdata -1 _R1.fq -2 _R2.fq -u ALIGN -o first_ep_run >> logfile.txt 2>&1

# Test that EP runs and produces and alignment and phylogeny
# Tests flags: -a, -d, -1, -2 , -u PHYLO, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -a ../testdata/combo.fas -d ../testdata -1 _R1.fq -2 _R2.fq -u PHYLO -o second_ep_run >> logfile.txt 2>&1
