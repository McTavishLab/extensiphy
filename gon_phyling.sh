#!/usr/bin/env bash
# written by Jasper Toscani Field

# user specifies the directory, currently hard coding for files that end with R1_001.fastq.gz
# this will change in future version but currently, users must alter their files to end in R1_001.fastq.gz and R2_001.fastq.gz as we are assuming paired end reads

GON_PHYLING=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $1

# WD=$(pwd)
# while getopts ":b:d:g:r:c:1:2:h" opt; do
#   case $opt in
#     b) bootstrapping="$OPTARG"
#     ;;
#     d) read_dir="$OPTARG"
# 	  ;;
# 	  g) ref_genome="$OPTARG"
#     ;;
#     r) runs="$OPTARG"
#     ;;
#     c) threads="$OPTARG"
#     ;;
#     1) r1_tail="$OPTARG"
#     ;;
#     2) r2_tail="$OPTARG"
#     ;;
#     h) echo  "alignment in fasta format (-a), tree in Newick format (-t), and reads in fastq (-p -e paired_end_base_filenames or -s single_end_base_filename required)"
#     exit
#     ;;
#     \?) echo "Invalid option -$OPTARG" >&2
#     ;;
#   esac
# done


printf "ref_genome = $ref_genome"

# go to the directory containing the reads
# this is specified in the .cfg file
cd $read_dir

# match files to eachother and enter them into the processing programs
for i in $(ls *$r1_tail); do
   echo fastq "$i" "${i%$r1_tail}$r2_tail"
done

printf "made it through changing into the read directory and fastqc"

mkdir trimmed_reads

trim_pwd=$(pwd)

ls *$r1_tail | split -d -l $gon_phy_runs

for j in $(ls x*); do
  for i in $(cat $j); do

# for i in $(ls *$r1_tail); do
    bbduk.sh in1=$trim_pwd/"$i" in2=$trim_pwd/"${i%$r1_tail}$r2_tail" ref=adapters ktrim=r trimq=10 out=$trim_pwd/trimmed_reads/"$i" out2=$trim_pwd/trimmed_reads/"${i%$r1_tail}$r2_tail" stats=trimmed_stats"$i".out &
  done
  wait
done

wait

printf "made it through bbduk step"

# begin moving trimmed reads into a seperate directory so those files can be worked on
#mkdir trimmed_reads

#mv clean* ./trimmed_reads

cd ./trimmed_reads

# begin assembly section with spades
check=$(ls *.gz | wc -l)


if [[ $check -ge 1 ]]; then

  gzip -d *.gz

else
  ls
fi


mkdir spades_output

ls *$r1_tail | split -d -l $gon_phy_runs

for j in $(ls x*); do
  for i in $(cat $j); do

#cd ./spades_output
    printf "Read processing complete. Beginning spades.py assembly"
    if [[ $threads =~ ^-?[0-9]+$ ]]; then

    # for loop to make a seperate output directory for each set of reads
    	# for i in $(ls *$r1_tail); do


    # spades.py -1 <first/left read file> -2 <second/right read file> -t <threads> -o <output directory>
    # for i in $(ls *R1_001.fastq); do
    		printf "%s threads selected" "$threads"
    		time spades.py -1 "$i" -2 "${i%$r1_tail}$r2_tail" -t $threads -o ./spades_output/$i &
    	# done
    else
    	# for i in $(ls *$r1_tail); do

    		printf "No extra threads requested, default: 1"
    		time spades.py -1 "$i" -2 "${i%$r1_tail}$r2_tail" -t 1 -o ./spades_output/$i &

    	# done
    fi
  done
  wait
done
cd ./spades_output

### TODO ADD REPETITIVE SEQUENCE MASKING FUNCTION HERE


mkdir genomes_for_parsnp

for i in $( ls -d *);
do
        echo $i
        cd $i ;
        pwd ;

        cp contigs.fasta  ../genomes_for_parsnp/$i.fasta
        cd .. ;

done

#location for repetitive sequence masker

cd ./genomes_for_parsnp

mkdir excess_files

mkdir masked_genomes

#
# if [ "$repetitive" -eq 1 ]; then
# 	echo "BEGINNING REPETITIVE SEQUENCE MASKING"
# 	for i in $( ls *.fasta);
# 		do
# 		trf409.legacylinux64 $i 2 7 7 80 10 50 500 -h -m
# 		mv *.mask ./masked_genomes
#
# 	done
#
# 	parsnp -c -p 6 -d ./masked_genomes -r $ref_genome
#
#         mkdir alignment_fixing
#
#         cp ./P*/parsnp.xmfa ./alignment_fixing/
#
#         cd ./alignment_fixing
#
#
# # having to harcode the par selecting parsnp_splitter.py as the file cant seem to be found in path despite other files
# # in same folder being found in path
#         $GON_PHYLING/parsnp_splitter.py parsnp.xmfa
#
# # call RAxML for phylogenetic analysis on loci
#         if [[ "$var" =~ ^-?[0-9]+$ ]]; then
#
#                 printf "%s threads requested" "$threads"
#                 raxmlHPC-PTHREADS -f a -p 23456 -s ./combo.fas -x 23456 -# 100 -m GTRGAMMA -n core_genome_run.out -T $threads
#
#         else
#                 printf "no additional threads requested, using default"
#                 raxmlHPC-PTHREADS -f a -p 23456 -s ./combo.fas -x 23456 -# 100 -m GTRGAMMA -n core_genome_run.out
#         fi
#
# else
echo "SKIPPING REPETITIVE SEQUENCE MASKING AND PROCEEDING WITH PARSNP"

printf "$ref_genome is ref genome selection"

# CHECKING FOR REFERENCE USE
if [ $ref_genome == "NONE" ]; then

  parsnp -c -p $threads -d ./ -r !

elif [ $ref_genome != "NONE" ]; then
  parsnp -c -p $threads -d ./ -r $ref_genome

fi

mkdir alignment_fixing

cp ./P*/parsnp.xmfa ./alignment_fixing/

cd ./alignment_fixing


# old parsnp processing command
# $GON_PHYLING/parsnp_splitter.py parsnp.xmfa
if [ $output_type == "LOCUS" ]; then

	mkdir locus_msa_files

       cat parsnp.xmfa | grep -Po "cluster\d+" | sort | uniq > ./locus_msa_files/locus_IDs.txt
	
	cd ./locus_msa_files

	echo pwd
	cat ./locus_IDs.txt | split -d -l $threads

	

	for j in $(ls x*); do
		for i in $(cat $j); do
			$GON_PHYLING/locus_splitter.py --align_file ../parsnp.xmfa --out_file ./$i-.fasta --locus_id $i
		done
		wait
	done
	cd ..

elif [ $output_type == "LOCI" ]; then
	printf "YOU SELECTED TO NOT SPLIT THE LOCI INTO SEPERATE FILES\n"
fi

# new parallel parsnp processing
# grab the number of taxa in the parsnp.xmfa output file
$GON_PHYLING/parsnp_taxa_count.sh > taxa_count.txt

# split each number of taxa into a set that can be iterrated over and processed in parallel
cat taxa_count.txt | split -a 10 -l $threads

for j in $(ls xa*); do
  for i in $(cat $j); do

    $GON_PHYLING/parallel_parsnp_splitter.py parsnp.xmfa $i &

  done
  wait
done

# combine the split parsnp sequence files into a single file
cat parsnp_chunk-*-.fa > combo.fas

# remove the first newline character and turn combo.txt into a fasta file
sed -i '1d' combo.fas

# trim down the taxa names
sed -i 's/_[^_]*//2g' combo.fas

if [ $bootstrapping == "ON" ]; then

  # call RAxML for phylogenetic analysis on loci
  if [[ $threads =~ ^-?[0-9]+$ ]]; then

  	printf "%s threads requested" "$threads"
  	time raxmlHPC-PTHREADS -f a -p 23456 -s ./combo.fas -x 23456 -# 100 -m GTRGAMMA -n core_genome_run.out -T $threads

    # raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s ./combo.fas -p 12345 -n core_genome_run.out

  else
  	printf "no additional threads requested, using default"
  	time raxmlHPC-PTHREADS -f a -p 23456 -s ./combo.fas -x 23456 -# 100 -m GTRGAMMA -n core_genome_run.out

    # raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s ./combo.fas -p 12345 -n core_genome_run.out
  fi

elif [ $bootstrapping == "OFF" ]; then

  # call RAxML for phylogenetic analysis on loci
  if [[ $threads =~ ^-?[0-9]+$ ]]; then

    printf "%s threads requested" "$threads"
    # raxmlHPC-PTHREADS -f a -p 23456 -s ./combo.fas -x 23456 -# 100 -m GTRGAMMA -n core_genome_run.out -T $threads

    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s ./combo.fas -p 12345 -n core_genome_run.out

  else
    printf "no additional threads requested, using default"
    # raxmlHPC-PTHREADS -f a -p 23456 -s ./combo.fas -x 23456 -# 100 -m GTRGAMMA -n core_genome_run.out

    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s ./combo.fas -p 12345 -n core_genome_run.out
  fi

else
  printf "Switch bootstrapping option to 'ON' or 'OFF' and re-run program."

fi

printf "Assembly and inference of genomes complete."
