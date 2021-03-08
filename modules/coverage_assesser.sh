#! /bin/bash
set -e
set -u
set -o pipefail

# establishes the path to find the phycorder directory
PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

while getopts ":e:s:h" opt; do
  case $opt in
    e) ep_output="$OPTARG"
    ;;
    s) samtools_loc="$OPTARG"
    ;;
    h) printf  " Finds coverage of Extensiphy consensus sequences\n
    \n (-e) Extensiphy output directory (full path)
    \n (-s) Full path to samtools file\n"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done



samtools_cov_location=$samtools_loc
multi_folder_link_dir=$ep_output

# ls ${multi_folder_link_dir}

for i in $(ls -d $multi_folder_link_dir/*output_dir );
do
    file=$(basename $i output_dir)
    printf "\n$file\n"
    cd $i
    file_path=$(pwd)
    ln -sf ${file_path}/best_sorted.bam $multi_folder_link_dir/$file.bam
    cd $multi_folder_link_dir
done

touch coverage_outputs.txt

for i in $(ls $multi_folder_link_dir/*.bam );
do
    taxon=$(basename $i .bam)
    echo "$taxon" >> coverage_outputs.txt
    $samtools_cov_location coverage -m $i >> coverage_outputs.txt
done