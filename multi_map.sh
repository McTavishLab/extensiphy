#! /bin/bash
# inputs: a sequence alignment and a tree inferred from that alignment.
# inputs: a directory of paired end reads for new taxa to be added to the alignment and corresponding tree.
# example command: multi_map.sh -example.aln -t example.tre -p example_read_dir -c number_of_threads

set -e
set -u
set -o pipefail


PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $PHYCORDER/phycorder.cfg

# while getopts ":a:t:p:o:n:m:b:w:c:r:h" opt; do
#   case $opt in
#     a) align="$OPTARG"
#     ;;
#     t) tree="$OPTARG"
#     ;;
#     p) read_dir="$OPTARG"
#     ;;
#     o) outdir="$OPTARG"
#   	;;
#   	n) nam="$OPTARG"
#     ;;
#     m) map="$OPTARG"
#     ;;
#     b) re_map="$OPTARG"
#     ;;
#     w) wre_map="$OPTARG"
#     ;;
#     c) threads="$OPTARG"
#     ;;
#     r) phycorder_runs="$OPTARG"
#     ;;
#     h) echo  "alignment in fasta format (-a), tree in Newick format (-t), directory of reads (-p), number of reads your computer can process at one time (-r)"
#     exit
#     ;;
#     \?) echo "Invalid option -$OPTARG" >&2
#     ;;
#   esac
# done
#
# if [ -z "$align" ] || [ -z "$tree" ]; then
#    "alignment (-a), tree (-t), and reads (-p or -s required)"
#    exit
# fi
#
# #Ttest if files actually exist
# #Check to make sure mapping has occured if re-mapping
#
# if [ -f "$align" ]; then
#     printf "Alignment is %s\n" "$align"
#   else
#     printf "Alignment $align not found. Exiting\n" >&2
#     exit
# fi
# if [ -f "$tree" ]; then
#     printf "Tree is %s\n" "$tree"
#   else
#     printf "Tree $tree not found. Exiting\n" >&2
#     exit
# fi
# if [ -d "$read_dir" ]; then
#     printf "Directory of reads is %s\n" "$read_dir"
#   else
#     printf "Directory of reads $read_dir not found. exiting\n" >&2
#     exit
# fi
#
# #this is a hack that is in both scripts!! need to be passed between
# r1_tail="R1_001.fastq"
# r2_tail="R2_001.fastq"

mkdir -p $outdir
cd $outdir

ls ${read_dir}/*$r1_tail | split -a 5 -l $phycorder_runs

#if [ $threads -ge $num_files ]
# then
  printf "Number of cores allocated enough to process all read sets\n"
  printf "Beginning Phycorder runs\n"

for j in $(ls xa*); do
  for i in $(cat $j); do
      base=$(basename $i $r1_tail)
      echo $base
      echo $i
      echo $PHYCORDER
      echo $align
      echo $tree x
      echo $i
      echo ${base}${r2_tail}
      echo $threads
      echo "${base}_output_dir"
      echo "$PHYCORDER/map_to_align.sh -a $align -t $tree -p $i -e ${i%$r1_tail}$r2_tail -c $threads -o ${base}output_dir > parallel-$base-dev.log &"
      time $PHYCORDER/map_to_align.sh -a $align -t $tree -p $i -e ${i%$r1_tail}$r2_tail -c $threads -o ${base}output_dir > parallel-$base-dev.log &
      #wait
      printf "adding new map_to_align run"
  done
  wait
done


 # else
   # section handles times when you are adding more tips to the tree than available processors

#    echo "Number of taxa being added to alignment and tree are greater than number of processors\n"
#    echo "Beginning job-number controlled Phycorder run\n"
#    echo "THIS DOES NOT APPEAR TO BE WORKING. only more cores than new runs, for now\n"
#    for i in $(cat readnames.txt); do
#     base=$(basename $i $r1_tail)
# #     # time $PHYCORDER/map_to_align.sh -a $align -t $tree -p "$read_dir"/"$i" -e "$read_dir"/"${i%R1_.fastq}R2_.fastq" -c $threads -o "$i"_"output_dir" > "$PHYCORDER/multi_map_dev.log" &
#      time $PHYCORDER/map_to_align.sh -a $align -t $tree -p $i  -e ${i%$r1_tail}$r2_tail -c $threads -o "${base}_output_dir" > "multi_map_dev.log" &
#      printf "adding new map_to_align run\n"
#      while [ $(jobs | wc -l) -ge $threads ] ; do sleep 1 ; done
#    done
# fi
   printf "Individual Phycorder runs finished. Combining aligned query sequences and adding them to starting alignment\n"

   mkdir -p combine_and_infer

   mkdir -p phycorder-dev-logs

   for i in $(ls -d *_output_dir); do
      cp $i/*_align.fas combine_and_infer
   done


   cat combine_and_infer/*.fas $align > combine_and_infer/extended.aln

   printf "Extended alignment file creaded (extended.aln), using previous tree as starting tree for phylogenetic inference\n"

   raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s combine_and_infer/extended.aln -t $tree -p 12345 -n consensusFULL

   printf "Multiple taxa update of phylogenetic tree complete\n"
   printf "Moving run logs into phycorder-dev-logs"
   cd ..

   mv *-dev.log $outdir/phycorder-dev-logs
