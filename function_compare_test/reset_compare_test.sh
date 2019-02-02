#! /bin/bash
# script resets phycorder test space for fast rerunning

PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

rm -r ./remaining_reads_for_phycorder

rm -r ./pipeline_results

rm -r ./phycorder_compare_test

rm -r ./first_5_for_phycorder

cd $PHYCORDER/fastq_files

rm -r trimmed*
