#! /bin/bash
# TEST to compare "traditional" methods of loci selection and phylogenetic analyses
# and rapid updating using phycorder

set -e
set -u
set -o pipefail

# get phycorder directory location for pathing
PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

mkdir first_5_for_phycorder

$PHYCORDER/../taxon_splitter.py -m --max_num 5 --taxa_dir $PHYCORDER/fastq_files --new_dir $PHYCORDER/first_5_for_phycorder

$PHYCORDER/../gon_phyling.sh $PHYCORDER/gon_phy_first_5.cfg

$PHYCORDER/../multi_map.sh $PHYCORDER/phycorder_main_run.cfg

$PHYCORDER/../gon_phyling.sh $PHYCORDER/gon_phy_main_run.cfg
