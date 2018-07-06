#! /bin/bash

PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

while getopts ":a:t:p:o:n:r:m:b:w:c:h" opt; do
  case $opt in
    a) align="$OPTARG"
    ;;
    t) tree="$OPTARG"
    ;;
    p) read_dir="$OPTARG"
    ;;
    o) outdir="$OPTARG"
	;;
	n) nam="$OPTARG"
    ;;
    r) read_align="$OPTARG"
    ;;
    m) map="$OPTARG"
    ;;
    b) re_map="$OPTARG"
    ;;
    w) wre_map="$OPTARG"
    ;;
    c) threads="$OPTARG"
    ;;
    h) echo  "alignment in fasta format (-a), tree in Newick format (-t), and reads in fastq (-p -e paired_end_base_filenames or -s single_end_base_filename required)"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z "$align" ] || [ -z "$tree" ]; then
   "alignment (-a), tree (-t), and reads (-p or -s required)"
   exit
fi

#Ttest if files actually exist
#Check to make sure mapping has occured if re-mapping

if [ -f "$align" ]; then
    printf "Alignment is %s\n" "$align"
  else
    printf "Alignment $align not found. Exiting\n" >&2
    exit
fi
if [ -f "$tree" ]; then
    printf "Tree is %s\n" "$tree"
  else
    printf "Tree $tree not found. Exiting\n" >&2
    exit
fi
# if [ $PE -eq 1 ]; then
#   if [ -f ${read_one} ]; then
#      printf "Paired end reads \n"
#      printf "read one is ${read_one}\n"
#   else
#     printf "read one ${read_one} not found. Exiting\n" >&2
#     exit
# fi
#   if [ -f ${read_two} ]; then
#      printf "Paired end reads \n"
#      printf "read two is ${read_two}\n"
#   else
#     printf "read two ${read_two} not found. Exiting\n" >&2
#     exit
# fi
# fi
printf "Argument out is %s\n" "$outdir"
printf "Argument name is %s\n" "$nam"
printf "Argument map is %s\n" "$map"
printf "Argument re_mapis %s\n" "$re_map"

mkdir -p $outdir

#cd $READ_DIR
#READ_LOC=$(pwd)
for i in $(ls *R1_.fastq); do
    $PHYCORDER/map_to_align.sh -a $align -t $tree -p "$read_dir/$i" -e "$read_dir{$i%R1_.fastq}R2_.fastq" -t $threads -o $outdir > "$PHYCORDER/multi_map_dev.log"
done

wait

cd $outdir

cat *aligned_cns.fas > $outdir/extended.aln

raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -T $threads -s extended.aln -t $tree -p 12345 -n consensusFULL
