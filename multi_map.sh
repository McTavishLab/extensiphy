#! /bin/bash
# inputs: a sequence alignment and a tree inferred from that alignment.
# inputs: a directory of paired end reads for new taxa to be added to the alignment and corresponding tree.
# example command: ./multi_map.sh ./sample_phycorder.cfg

set -e
set -u
set -o pipefail


PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# changing location of .cfg file to a variable
# for easy use of multiple config num_files
source $1


mkdir -p $outdir
#cd $outdir
#outdir_path=$(pwd)

# cd $read_dir
#
# for i in $(ls *$r1_tail); do
#   echo ">$i" >> "original_file_names_phycorder.txt"
# done
#
# mv "original_file_names_phycorder.txt" $outdir/

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
      echo "$PHYCORDER/map_to_align.sh -a $align -t $tree -p $i -e ${i%$r1_tail}$r2_tail -1 $r1_tail -2 $r2_tail -c $threads -o ${base}output_dir > parallel-$base-dev.log &"
      echo "Time for $j Phycorder run:"
      time $PHYCORDER/map_to_align.sh -a $align -t $tree -p $i -e ${i%$r1_tail}$r2_tail -1 $r1_tail -2 $r2_tail -c $threads -o ${base}output_dir > parallel-$base-dev.log &
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

   wd=$(pwd)

   for i in $(ls -d *output_dir); do
     cd $i
     count=$(ls *_align.fas | wc -l)
     if [ $count -gt 0 ]; then
       cp *_align.fas $wd/combine_and_infer/
     else
       echo "$i"
       continue
     fi
     cd ..
   done

   # if [ -f "$name_dict" ]; then
   #   printf "entering renaming step based on previous dictionary non containing the new taxa"
   #
   #   $PHYCORDER/name_parser.py -u --dict_file $name_dict --newtaxa_dir combine_and_infer
   #
   #   cat combine_and_infer/OTU*.fas $align > combine_and_infer/extended.aln
   # else
   printf "skipping renaming step"
   cat combine_and_infer/*.fas $align > combine_and_infer/extended.aln
   # fi

   printf "Extended alignment file creaded (extended.aln), using previous tree as starting tree for phylogenetic inference\n"

   cd combine_and_infer

   # strip the unnecessary information from the taxa names in the alignment.
   # this assumes you've used the renaming tool to rename all of the reads for this experiment

   sed -i 's/_$//g' extended.aln

   # for i in $(cat < extended.aln); do
   #   echo "${i%_*}" >> extended2.aln
   # done

   # rm extended.aln
   # mv extended2.aln extended.aln


   # time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s extended.aln -t $tree -p 12345 -n consensusFULL

   time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s extended.aln -t $tree -p 12345 -n consensusFULL

   time raxmlHPC-PTHREADS -s extended.aln -n consensusFULL_bootstrap -m GTRGAMMA  -p 12345 -T $threads -N 100 -b 12345

   time raxmlHPC-PTHREADS -z RAxML_bootstrap.consensusFULL_bootstrap -t RAxML_bestTree.consensusFULL -f b -T $threads -m GTRGAMMA -n majority_rule_bootstrap_consensus

   printf "Multiple taxa update of phylogenetic tree complete\n"
   # printf "Moving run logs into phycorder-dev-logs"
   # cd ..
   #
   # mv *-dev.log $outdir/phycorder-dev-logs
