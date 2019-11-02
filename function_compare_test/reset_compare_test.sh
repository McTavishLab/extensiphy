#! /bin/bash
# script resets phycorder test space for fast rerunning

PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

rm -r $PHYCORDER/remaining_reads_for_phycorder

rm -r $PHYCORDER/pipeline_results

rm -r $PHYCORDER/phycorder_compare_test

rm -r $PHYCORDER/first_5_for_phycorder

cd $PHYCORDER/fastq_files

rm -r trimmed*
