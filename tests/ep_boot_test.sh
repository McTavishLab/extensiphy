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

###############################################################################
# Test that EP bootstraps a phylogeny thats being updated
# # Tests flags: -a, -d, -1, -2, -u PHYLO, -t, -b, -o
# Examine: extended.aln, RAxML_bestTree.consensusFULL, RAxML_bootstrap.consensusFULL_bootstrap, RAxML_bipartitionsBranchLabels.majority_rule_bootstrap_consensus
../extensiphy.sh -u PHYLO -b ON -a ../testdata/combo.fas -t ../testdata/combo.tre -d ../testdata -1 _R1.fq -2 _R2.fq -o ep_test_five >> logfile.txt 2>&1

ALIGN=./ep_test_five/RESULTS/extended.aln
PHYLO=./ep_test_five/RESULTS/RAxML_bestTree.consensusFULL
BOOT=./ep_test_five/RESULTS/RAxML_bipartitions.majority_rule_bootstrap_consensus
num_lines=$(wc -l ./ep_test_five/RESULTS/extended.aln)
num_seqs=$(grep -c ">" ./ep_test_five/RESULTS/extended.aln)
check_tree=$(grep -c ":0.0;" ./ep_test_five/RESULTS/RAxML_bestTree.consensusFULL)
check_boot=$(grep -c "taxon_17" ./ep_test_five/RESULTS/RAxML_bipartitions.majority_rule_bootstrap_consensus)

if [ ${num_seqs} == 23 ] && [ ${check_tree} -eq 1 ] && [ ${check_boot} -eq 1 ]
then
  echo "test alignment update and bootstrap phylo update: PASSED" >> test_results.txt
else
  echo "test alignment update and bootstrap phylo update: FAILED" >> test_results.txt
fi
