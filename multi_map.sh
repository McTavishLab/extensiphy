#! /bin/bash

PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

ALIGN=$1

TREE=$2

READ_DIR=$3

OUTPUT=$4

cd $READ_DIR
READ_LOC=$(pwd)
for i in $(ls *R1_.fastq); do
    "$PHYCORDER/map_to_align.sh -a $ALIGN -t $TREE -p "$READ_loc/$i" -e "{$i%R1_.fastq}R2_.fastq" -o $OUTPUT > "$PHYCORDER/multi_map_dev.log"
done
