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

# test_basic_run () {
#
#
#
# }

test_help_menu
