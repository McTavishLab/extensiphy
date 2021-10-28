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
touch logfile.txt


./ep_menu_test.sh >> test_results.txt


###############################################################################
# Test that EP runs and produces and updated alignment
# Tests flags: -a. -d, -1, -2, -u ALIGN, -o
# Examine: extended.aln

./ep_update_align_test.sh >> test_results.txt

###############################################################################
# Test that EP runs with single-end read data
# Tests flags: -a, -d, -1, -s SE, -u ALIGN, -o
# Examine: extended.aln
./ep_single_end_reads_test.sh >> test_results.txt


###############################################################################
# Test that EP runs and produces and alignment and phylogeny
# Tests flags: -a, -d, -1, -2 , -u PHYLO, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL

./ep_build_phylo_test.sh >> test_results.txt

###############################################################################
# Test that EP updates an existing phylogeny
# Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL

./ep_update_phylo_test.sh >> test_results.txt

###############################################################################
# Test that EP bootstraps a phylogeny thats being updated
# # Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -b, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL, RAxML_bootstrap.consensusFULL_bootstrap, RAxML_bipartitionsBranchLabels.majority_rule_bootstrap_consensus

./ep_boot_test.sh >> test_results.txt

################################################################################
# Tests selecting a specific reference taxon
# Tests flags: -a, -d, -1, -2, -u PHYLO -o, -t, -r
# Examine:  extended.aln, RAxML_bestTree.consensusFULL

./ep_ref_test.sh >> test_results.txt

################################################################################
# Tests using single locus alignment files as inputs instead of concatenated
# Tests flags: -a, -d, -1, -2, -u PHYLO, -o, -m
# Examine:  extended.aln, RAxML_bestTree.consensusFULL
# ../extensiphy.sh -u ALIGN -a ../testdata/single_locus_align_dir -d ../testdata -m SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o sixth_extensiphy_run >> logfile.txt 2>&1

./ep_single_locus_in_test.sh >> test_results.txt

###############################################################################
# Tests outputting updated single locus alignments
# Tests flags: -a, -d, -t, -m SINGLE_LOCUS_FILES, -u PHYLO, -g SINGLE_LOCUS_FILES, -1, -2, -o
# Examine:  extended.aln, RAxML_bestTree.consensusFULL, /locus_out_Extensiphy_run/RESULTS/updated_single_loci/single_locus_1_.fasta
ep_single_locus_out_test.sh >> test_results.txt

cat test_results.txt
