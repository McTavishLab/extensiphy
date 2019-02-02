#! /bin/bash
# TEST to compare "traditional" methods of loci selection and phylogenetic analyses
# and rapid updating using phycorder

set -e
set -u
set -o pipefail

# get phycorder directory location for pathing
PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

mkdir first_5_for_phycorder

$PHYCORDER/../taxon_splitter.py -m --max_num 5 --taxa_dir $PHYCORDER/fastq_files/ --new_dir $PHYCORDER/first_5_for_phycorder

$PHYCORDER/../gon_phyling.sh $PHYCORDER/gon_phy_first_5.cfg

mkdir remaining_reads_for_phycorder

cp ./fastq_files/*.fastq ./remaining_reads_for_phycorder/

$PHYCORDER/../taxon_splitter.py -d --max_num 5 --taxa_dir $PHYCORDER/remaining_reads_for_phycorder/

$PHYCORDER/../multi_map.sh $PHYCORDER/phycorder_main_run.cfg

$PHYCORDER/../gon_phyling.sh $PHYCORDER/gon_phy_main_run.cfg

mkdir pipeline_results

cp ./phycorder_compare_test/combine_and_infer/RAxML_bestTree.consensusFULL ./pipeline_results/

cp ./fastq_files/trimmed_reads/spades_output/genomes_for_parsnp/alignment_fixing/RAxML_bestTree.core_genome_run.out ./pipeline_results/

tree_compare_results.py --working_dir $PHYCORDER/pipeline_results --gon_phy_tree RAxML_bestTree.core_genome_run.out --phycorder_tree RAxML_bestTree.consensusFULL