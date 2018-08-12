#! /bin/bash
# inputs: a sequence alignment and a tree inferred from that alignment.
# inputs: a directory of paired end reads for new taxa to be added to the alignment and corresponding tree.
# example command: multi_map.sh -example.aln -t example.tre -p example_read_dir -c number_of_threads

set -e
set -u
set -o pipefail


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
if [ -d "$read_dir" ]; then
    printf "Directory of reads is %s\n" "$read_dir"
  else
    printf "Directory of reads $read_dir not found. exiting\n" >&2
    exit
fi


cd $read_dir

num_files=$(ls -1 *R1_.fastq | wc -l)

if [ $threads -ge $num_files]
then
  printf "Number of cores allocated enough to process all read sets\n"
  printf "Beginning Phycorder runs\n"

  for i in $(ls *R1_.fastq); do
      time $PHYCORDER/map_to_align.sh -a $align -t $tree -p "$read_dir"/"$i" -e "$read_dir"/"${i%R1_.fastq}R2_.fastq" -c $threads -o "$i"_"output_dir" > "$PHYCORDER/multi_map_dev.log" &
      printf "adding new map_to_align run"
  done

  wait

  printf "Individual Phycorder runs finished. Combining aligned query sequences and adding them to starting alignment\n"

  mkdir combine_and_infer

  for i in $(ls -d *_output_dir); do
    cp $i/*R1*.fas ./
  done

  cp $align $read_dir

  cat *.fas > extended.aln

  printf "Extended alignment file creaded (extended.aln), using previous tree as starting tree for phylogenetic inference\n"

  raxmlHPC-PTHREADS-AVX -m GTRGAMMA -T $threads -s extended.aln -t $tree -p 12345 -n consensusFULL

  printf "Multiple taxa update of phylogenetic tree complete\n"
else
  # section handles times when you are adding more tips to the tree than available processors
  echo "Number of taxa being added to alignment and tree are greater than number of processors\n"
  echo "Beginning job-number controlled Phycorder run\n"

  for i in $(ls *R1_.fastq); do
    time $PHYCORDER/map_to_align.sh -a $align -t $tree -p "$read_dir"/"$i" -e "$read_dir"/"${i%R1_.fastq}R2_.fastq" -c $threads -o "$i"_"output_dir" > "$PHYCORDER/multi_map_dev.log" &
    printf "adding new map_to_align run\n"
    while [ $(jobs | wc -l) -ge $threads ] ; do sleep 1 ; done
  done

  wait

  printf "Individual Phycorder runs finished. Combining aligned query sequences and adding them to starting alignment\n"

  mkdir combine_and_infer

  for i in $(ls -d *_output_dir); do
    cp $i/*R1*.fas ./
  done

  cp $align $read_dir

  cat *.fas > extended.aln

  printf "Extended alignment file creaded (extended.aln), using previous tree as starting tree for phylogenetic inference\n"

  raxmlHPC-PTHREADS-AVX -m GTRGAMMA -T $threads -s extended.aln -t $tree -p 12345 -n consensusFULL

  printf "Multiple taxa update of phylogenetic tree complete\n"

fi
