#! /bin/bash
# TEST to compare "traditional" methods of loci selection and phylogenetic analyses
# and rapid updating using phycorder

set -e
set -u
set -o pipefail

# get phycorder directory location for pathing
PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

compare_path = function_compare_test

cd $PHYCORDER

cd function_compare_test

mkdir first_5_for_phycorder

$PHYCORDER/taxon_splitter.py -m --max_num 5 --taxa_dir $PHYCORDER/function_compare_test/fastq_files --new_dir $PHYCORDER/function_compare_test/first_5_for_phycorder

$PHYCORDER/gon_phyling.sh $PHYCORDER/function_compare_test/gon_phy_first_5.cfg

$PHYCORDER/multi_map.sh $PHYCORDER/function_compare_test/phycorder_main_run.cfg

$PHYCORDER/gon_phyling.sh $PHYCORDER/function_compare_test/gon_phy_main_run.cfg
