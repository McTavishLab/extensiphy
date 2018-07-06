#! /bin/bash

PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ALIGN=$1

TREE=$2

READ_DIR=$3

OUTPUT=$4

THREADS=$5

cd $READ_DIR
READ_LOC=$(pwd)
for i in $(ls *R1_.fastq); do
    "$PHYCORDER/map_to_align.sh -a $ALIGN -t $TREE -p "$READ_LOC/$i" -e "$READ_LOC{$i%R1_.fastq}R2_.fastq" -t $THREADS -o $OUTPUT > "$PHYCORDER/multi_map_dev.log"
done

wait

cd $OUTPUT

cat *aligned_cns.fas > $OUTPUT/extended.aln

raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -T $THREADS -s extended.aln -t $TREE -p 12345 -n consensusFULL
