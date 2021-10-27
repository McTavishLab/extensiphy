#! /bin/bash
# Comprehensive test of Extensiphy options and outputs
# establishes the path to find the phycorder directory

set -e
#-e Exit immediately if a pipeline (see Pipelines), which may consist of a
#single simple command (see Simple Commands), a list (see Lists), or a compound
# command (see Compound Commands) returns a non-zero status. The shell does not
# exit if the command that fails is part of the command list immediately
# following a while or until keyword, part of the test in an if statement,
# part of any command executed in a && or || list except the command following
# the final && or ||, any command in a pipeline but the last, or if the
#command’s return status is being inverted with !. If a compound command other
#than a subshell returns a non-zero status because a command failed while -e was
# being ignored, the shell does not exit. A trap on ERR, if set,
#is executed before the shell exits.

set -u
#-u Treat unset variables and parameters other than the special parameters
# ‘@’ or ‘*’ as an error when performing parameter expansion.
#An error message will be written to the standard error, and
#a non-interactive shell will exit.

set -o pipefail
# -o pipefail If set, the return value of a pipeline is the value of the last
#(rightmost) command to exit with a non-zero status, or zero if all commands in
#the pipeline exit successfully. This option is disabled by default.


touch test_results.txt

../extensiphy.sh -h > test_help.txt
## Question do we want the help mneu to crash if programs not installed?

# Check that some expected text is there
VAR=$(grep 'Extensiphy is a program for quickly adding genomic sequence' test_help.txt | wc -l)

if [[ $VAR -eq 1 ]]
then
  echo "text found in help output."
  echo "test help menu: PASSED" >> test_results.txt
else
  echo "text found in help output"
  echo "test help menu: FAILED" >> test_results.txt
fi


touch logfile.txt

###############################################################################
# Test that EP runs and produces and updated alignment
# Tests flags: -a. -d, -1, -2, -u ALIGN, -o
# Examine: extended.aln
../extensiphy.sh -u ALIGN -a ../testdata/combo.fas -d ../testdata -1 _R1.fq -2 _R2.fq -o first_extensiphy_run >> logfile.txt 2>&1

ALIGN=./first_extensiphy_run/RESULTS/extended.aln

num_seqs=$(grep -c ">" ./first_extensiphy_run/RESULTS/extended.aln)
num_lines=$(wc -l ./first_extensiphy_run/RESULTS/extended.aln)
num_chars=$(wc ./first_extensiphy_run/RESULTS/extended.aln | cut -f11 -d ' ')

if [[ ${num_seqs} -eq 23 ]]
then
  echo "text found in help output."
  echo "test simple alignment update: PASSED" >> test_results.txt
else
  echo "text found in help output"
  echo "test simple alignment update: FAILED" >> test_results.txt
fi

###############################################################################
# Test that EP runs and produces and alignment and phylogeny
# Tests flags: -a, -d, -1, -2 , -u PHYLO, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -u PHYLO -a ../testdata/combo.fas -d ../testdata -1 _R1.fq -2 _R2.fq -o second_extensiphy_run >> logfile.txt 2>&1

ALIGN=./second_extensiphy_run/RESULTS/extended.aln
PHYLO=./second_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL
num_lines=$(wc -l ./second_extensiphy_run/RESULTS/extended.aln)
check_tree=$(grep -c ":0.0;" ./second_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL)

if [ ${num_seqs} == 23 ] && [ ${check_tree} -eq 1 ]
then
  echo "text found in help output."
  echo "test alignment update and phylo build: PASSED" >> test_results.txt
else
  echo "text found in help output"
  echo "test alignment update and phylo build: FAILED" >> test_results.txt
fi

###############################################################################
# Test that EP updates an existing phylogeny
# Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -u PHYLO -a ../testdata/combo.fas -t ../testdata/combo.tre -d ../testdata -1 _R1.fq -2 _R2.fq -o third_extensiphy_run >> logfile.txt 2>&1

ALIGN=./third_extensiphy_run/RESULTS/extended.aln
PHYLO=./third_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL
num_lines=$(wc -l ./third_extensiphy_run/RESULTS/extended.aln)
check_tree=$(grep -c ":0.0;" ./third_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL)

if [ ${num_seqs} == 23 ] && [ ${check_tree} -eq 1 ]
then
  echo "text found in help output."
  echo "test alignment update and phylo update: PASSED" >> test_results.txt
else
  echo "text found in help output"
  echo "test alignment update and phylo update: FAILED" >> test_results.txt
fi


###############################################################################
# Test that EP bootstraps a phylogeny thats being updated
# # Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -b, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL, RAxML_bootstrap.consensusFULL_bootstrap, RAxML_bipartitionsBranchLabels.majority_rule_bootstrap_consensus
../extensiphy.sh -u PHYLO -b ON -a ../testdata/combo.fas -t ../testdata/combo.tre -d ../testdata -1 _R1.fq -2 _R2.fq -o fourth_extensiphy_run >> logfile.txt 2>&1

ALIGN=./fourth_extensiphy_run/RESULTS/extended.aln
PHYLO=./fourth_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL
BOOT=./fourth_extensiphy_run/RESULTS/RAxML_bipartitions.majority_rule_bootstrap_consensus
num_lines=$(wc -l ./fourth_extensiphy_run/RESULTS/extended.aln)
check_tree=$(grep -c ":0.0;" ./fourth_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL)
check_boot=$(grep -c "taxon_17" ./fourth_extensiphy_run/RESULTS/RAxML_bipartitions.majority_rule_bootstrap_consensus)

if [ ${num_seqs} == 23 ] && [ ${check_tree} -eq 1 ] && [ ${check_boot} -gt 1 ]
then
  echo "text found in help output."
  echo "test alignment update and bootstrap phylo update: PASSED" >> test_results.txt
else
  echo "text found in help output"
  echo "test alignment update and bootstrap phylo update: FAILED" >> test_results.txt
fi



################################################################################
# Tests selecting a specific reference taxon
# Tests flags: -a, -d, -1, -2, -u PHYLO -o, -t, -r
# Examine:  extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -a ../testdata/combo.fas -u PHYLO -d ../testdata -t ../testdata/combo.tre -1 _R1.fq -2 _R2.fq -o fifth_extensiphy_run -r taxon_11 >> logfile.txt 2>&1

ALIGN=./fifth_extensiphy_run/RESULTS/extended.aln
PHYLO=./fifth_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL
num_lines=$(wc -l ./fifth_extensiphy_run/RESULTS/extended.aln)
check_tree=$(grep -c ":0.0;" ./fifth_extensiphy_run/RESULTS/RAxML_bestTree.consensusFULL)


if [ ${num_seqs} == 23 ] && [ ${check_tree} -eq 1 ]
then
  echo "text found in help output."
  echo "test alignment update and bootstrap phylo update for specific reference: PASSED" >> test_results.txt
else
  echo "text found in help output"
  echo "test alignment update and bootstrap phylo update for specific reference: FAILED" >> test_results.txt
fi

exit
################################################################################
# Tests using single locus alignment files as inputs instead of concatenated
# Tests flags: -a, -d, -1, -2, -u ALIGN, -o, -m
# Examine:  extended.aln
../extensiphy.sh -u ALIGN -a ../testdata/single_locus_align_dir -d ../testdata -m SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o sixth_extensiphy_run >> logfile.txt 2>&1



################################################################################
# Tests building a tree from single locus alignment files
# Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -o -m SINGLE_LOCUS_FILES
# Examine:  extended.aln, RAxML_bestTree.consensusFULL
../extensiphy.sh -u PHYLO -a ../testdata/single_locus_align_dir -d ../testdata -t ../testdata/combo.tre -m SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o split_extensiphy_run >> logfile.txt 2>&1




###############################################################################
# Tests outputting updated single locus alignments
# Tests flags: -a, -d, -t, -m SINGLE_LOCUS_FILES, -u PHYLO, -g SINGLE_LOCUS_FILES, -1, -2, -o
# Examine:  extended.aln, RAxML_bestTree.consensusFULL, /locus_out_Extensiphy_run/RESULTS/updated_single_loci/single_locus_1_.fasta
../extensiphy.sh -u PHYLO -a ../testdata/single_locus_align_dir -d ../testdata -t ../testdata/combo.tre -m SINGLE_LOCUS_FILES -g SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o locus_out_Extensiphy_run >> logfile.txt 2>&1
