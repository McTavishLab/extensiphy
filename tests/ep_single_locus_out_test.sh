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


../extensiphy.sh -u PHYLO -a ../testdata/single_locus_align_dir -d ../testdata -t ../testdata/combo.tre -m SINGLE_LOCUS_FILES -g SINGLE_LOCUS_FILES -1 _R1.fq -2 _R2.fq -o ep_test_single_locus_out >> logfile.txt 2>&1

ALIGN=./ep_test_single_locus_out/RESULTS/extended.aln
PHYLO=./ep_test_single_locus_out/RESULTS/RAxML_bestTree.consensusFULL
LOCUS_ALIGN=./ep_test_single_locus_out/RESULTS/updated_single_loci/
num_seqs=$(grep -c ">" ./ep_test_single_locus_out/RESULTS/extended.aln)
single_locus_align=$(cat ./ep_test_single_locus_out/RESULTS/updated_single_loci/single_locus_1_.fasta | wc -l)
num_lines=$(wc -l ./ep_test_single_locus_out/RESULTS/extended.aln)
check_tree=$(grep -c ":0.0;" ./ep_test_single_locus_out/RESULTS/RAxML_bestTree.consensusFULL)

if [ ${num_seqs} -eq 23 ] && [ ${check_tree} -eq 1 ] && [ ${single_locus_align} -eq 46]
then
  echo "test output single locus alignment files: PASSED"
else
  echo "test output single locus alignment files: FAILED"
fi
