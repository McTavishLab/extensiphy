#!/usr/bin/env bash
# written by Jasper Toscani Field

# user specifies the directory, currently hard coding for files that end with R1_001.fastq.gz
# this will change in future version but currently, users must alter their files to end in R1_001.fastq.gz and R2_001.fastq.gz as we are assuming paired end reads

set -e
set -u
set -o pipefail

GON_PHYLING=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#source $1

#Check for dependencies
if [ $(which parsnp | wc -l) -lt 1 ]
    then
        printf "Requires parsnp" >&2
     #   exit 0
    else
        printf "Correct version of parsnp found.\n"
fi
if [  $(which spades.py | wc -l) -lt 1 ] #TODO steup for greater than 1.2? this  is a sloppppy approach
    then
        printf "Requires spades.py" >&2
      #  exit 0
    else
        printf "Correct version of spades.py found.\n"
fi
if [ $(which raxmlHPC-PTHREADS | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "raxmlHPC-PTHREADS not found\n" >&2
    else
        printf "raxmlHPC-PTHREADS found\n"
fi
if [ $(which bbduk.sh | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "bbduk.sh not found. Install and/or add to path\n" >&2
    else
        printf "bbduk.sh found\n"
fi

printf "\n\n"

intermediate="KEEP"
ref_genome="NONE"
gon_phy_runs="2"
threads="2"
bootstrapping="OFF"
output_type="LOCI"
loci_positions="gon_phy_locus_positions.csv"
r1_tail="R1.fastq"
r2_tail="R2.fastq"
loci_len="700"

WD=$(pwd)
while getopts ":b:d:g:r:c:1:2:l:o:n:i:h" opt; do
   case $opt in
     b) bootstrapping="$OPTARG"
     ;;
     d) read_dir="$OPTARG"
     ;;
     g) ref_genome="$OPTARG"
     ;;
     r) gon_phy_runs="$OPTARG"
     ;;
     c) threads="$OPTARG"
     ;;
     1) r1_tail="$OPTARG"
     ;;
     2) r2_tail="$OPTARG"
     ;;
     l) loci_positions="$OPTARG"
     ;;
     o) output_type="$OPTARG"
     ;;
     i) intermediate="$OPTARG"
     ;;
     n) loci_len="$OPTARG"
     ;;
     h) printf  "gon_phyling is a program to assemble sets of paired-end short \
     high-throughput reads into genomes, automatically select homologous loci and \
     construct a phylogenetic tree.\n\n\n \
     EXAMPLE COMMAND: \n\n /path/to/gon_phyling.sh -d /path/to/read_directory -1 [READSET 1 SUFFIX] -2 [READSET 2 SUFFIX] \
     \n\n, INPUT OPTIONS: \
     \n (-d) directory of paired end reads. All output folders and files will be contained here \
     \n (-g) the name of the genome you wish to use as a reference during loci selection (if any)(DEFAULT: NONE) \
     \n (-1, -2) suffixes of paired-end input files in read directory (DEFAULT: -1 R1.fastq -2 R2.fastq) \
     \n\n OUTPUT \
     \n (-b) bootstrapping setting. Do you want to perform 100 boostrap replicates and add the support values to the best tree? (DEFAULT: OFF) \
     \n (-o) output type. Output either a concatenated multiple sequence alignment only or also output separate loci alignment files (DEFAULT: LOCI) (OPTIONS: LOCI, LOCUS) \
     \n (-l) Locus position file. Use if selecting -o LOCUS. Outputs a csv file tracking the loci names and their positions within the concatenated MSA (DEFAULT: gon_phy_locus_positions.csv) \
     \n\n RUNNING PROGRAM \
     \n (-r) gon_phyling runs. This is the number of genomes assembled at a single time (DEFAULT: 2) \
     \n (-c) Threads for each gon_phyling run. Figure out how many cores you have available and input [# of threads x # of parrallel genome assemblies] = cores you can allocate. (DEFAULT: 2), \
     \n (-i) Clean up intermediate output files option to save HD space (Options: CLEAN or KEEP)(DEFAULT: KEEP) \
     \n (-n) Cutoff of loci length when combining or separating loci as input or output (Options: int number)(DEFAULT: 700) \
     \n\n OUTPUT FILES \
     \n After a run, you will recieve the files: \n\ncombo.fas \nRAxML_best_tree.core_genome.out \
     \n\n these are your concatenated MSA and phylogenetic tree, respectively\n"
     exit
     ;;
     \?) echo "Invalid option -$OPTARG" >&2
     ;;
   esac
 done

tmp_read_dir=$(realpath $read_dir)
#tmp_ref_genome=$(realpath $ref_genome)

read_dir=$tmp_read_dir
#ref_genome=$tmp_ref_genome

printf "ref_genome = $ref_genome"

# go to the directory containing the reads
# this is specified in the .cfg file
cd $read_dir

# match files to eachother and enter them into the processing programs
for i in $(ls *$r1_tail); do
   echo fastq "$i" "${i%$r1_tail}$r2_tail"
done

#printf "made it through changing into the read directory and fastqc"

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

#printf "made it through bbduk step"

# begin moving trimmed reads into a seperate directory so those files can be worked on
#mkdir trimmed_reads

#mv clean* ./trimmed_reads

cd ./trimmed_reads

printf "\nCHANGED DIR TO TRIMMED_READS\n"

# begin assembly section with spades
#check=$(ls *.gz | wc -l)


#if [[ $check -ge 1 ]]; then

#  gzip -d *.gz

#else
#  ls
#fi

#check=$(ls -l | grep -o ".gz" | wc -l)

#if [[ $check != 0 ]]; then

#  printf "\n\n FOUND .GZ FILES!!\n\n"
#  printf "\n$check\n"
#  gzip -d *.gz

#else
#  printf "\nNO .GZ FILES FOUND!\n"
#fi


mkdir spades_output

ls *$r1_tail | split -d -l $gon_phy_runs

for j in $(ls x*); do
  for i in $(cat $j); do
    
    printf "\n$i\n"
    name_base=$(basename $i $r1_tail)
    printf "\nname_base\n"
#cd ./spades_output
    printf "Read processing complete. Beginning spades.py assembly"
    if [[ $threads =~ ^-?[0-9]+$ ]]; then

    # for loop to make a seperate output directory for each set of reads
    	# for i in $(ls *$r1_tail); do


    # spades.py -1 <first/left read file> -2 <second/right read file> -t <threads> -o <output directory>
    # for i in $(ls *R1_001.fastq); do
    		printf "%s threads selected" "$threads"
    		#time spades.py -1 "$i" -2 "${i%$r1_tail}$r2_tail" -t $threads -o ./spades_output/$i &
		time spades.py -1 "$i" -2 "${i%$r1_tail}$r2_tail" -t $threads -o ./spades_output/$name_base &
    	# done
    else
    	# for i in $(ls *$r1_tail); do

    		printf "No extra threads requested, default: 1"
    		#time spades.py -1 "$i" -2 "${i%$r1_tail}$r2_tail" -t 1 -o ./spades_output/$i &
		time spades.py -1 "$i" -2 "${i%$r1_tail}$r2_tail" -t 1 -o ./spades_output/$name_base &

    	# done
    fi
  done
  wait
done
cd ./spades_output

### TODO ADD REPETITIVE SEQUENCE MASKING FUNCTION HERE


mkdir genomes_for_parsnp

#for i in $( ls -d *);
#do
#        echo $i
#        cd $i ;
#        pwd ;
#
#        cp contigs.fasta  ../genomes_for_parsnp/$i.fasta
#        cd .. ;

#done

printf "\nOrganizing assembled genomes for Parsnp locus selection\n"

for i in $(ls -d */);
do
	printf "\n$i\n"
        slash_strip=$(basename $i /)
        if [[ $i == "genomes_for_parsnp/" ]]; then
                printf "\n$i SKIPPED\n"
        else
                printf "\nmoving contigs from $i\n"
                cd $i
                #cp contigs.fasta  ../genomes_for_parsnp/$slash_strip.fasta
		#ln -s contigs.fasta  ../genomes_for_parsnp/$slash_strip.fasta
		mv contigs.fasta  ../genomes_for_parsnp/$slash_strip.fasta
                cd ..
        fi
        #cd $i

        #if [[ $dir_check != 0 ]]; then
        #pwd ;

        #cp contigs.fasta  ../genomes_for_parsnp/$i.fasta
        #cd .. ;

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
#echo "SKIPPING REPETITIVE SEQUENCE MASKING AND PROCEEDING WITH PARSNP"

# STRIP ADDED UNNECESSARY FILE FORMAT INFO FROM FILE NAMES
#for i in $(ls -1); do mv $i $(echo $i | sed "s/_$r1_tail.fasta//") ; done

printf "$ref_genome is ref genome selection"

# CHECKING FOR REFERENCE USE
if [ $ref_genome == "NONE" ]; then

  parsnp -c -p $threads -d ./ -r ! 

elif [ $ref_genome != "NONE" ]; then
  
  tmp_ref_genome=$(realpath $ref_genome)
  ref_genome=$tmp_ref_genome
  parsnp -c -p $threads -d ./ -r $ref_genome

fi

mkdir alignment_fixing

sed -i -e 's/.fasta//g' ./P*/parsnp.xmfa

#symlinking this causes problems. Fix later
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

	
  if [ $loci_len != "700" ]; then
    printf "\nSplitting loci using user specificed cutoff of ${loci_len} nucleotides\n"
    for j in $(ls x*); do
		  for i in $(cat $j); do
			  $GON_PHYLING/modules/locus_splitter.py --align_file ../parsnp.xmfa --out_file ./$i-.fasta --locus_id $i --locus_size $loci_len
			  #$GON_PHYLING/limit_len_locus_splitter.py --align_file ../parsnp.xmfa --out_file ./$i-.fasta --locus_id $i --locus_size 1000
		  done
		  wait
	  done
  
  elif [ $loci_len == "700" ]; then
  
    printf "\nSplitting loci using default minimum cutoff of 700 nucleotides.\n"
	  for j in $(ls x*); do
		  for i in $(cat $j); do
			  $GON_PHYLING/modules/locus_splitter.py --align_file ../parsnp.xmfa --out_file ./$i-.fasta --locus_id $i --locus_size 700
			  #$GON_PHYLING/limit_len_locus_splitter.py --align_file ../parsnp.xmfa --out_file ./$i-.fasta --locus_id $i --locus_size 1000
		  done
		  wait
	  done
  fi

	#$GON_PHYLING/locus_combiner.py --msa_folder ./ --suffix .fasta --out_file ../combo.fas --position_dict_file $loci_positions 	
	$GON_PHYLING/modules/new_locus_combiner.py --msa_folder ./ --suffix .fasta --out_file ../combo.fas --position_csv_file $loci_positions --len_filter 50
	cd ..

	# STRIP TAILS FROM NAMES FOR CLEANER TAXON NAMES

elif [ $output_type == "LOCI" ]; then
	printf "YOU SELECTED TO NOT SPLIT THE LOCI INTO SEPERATE FILES\n"


# new parallel parsnp processing
# grab the number of taxa in the parsnp.xmfa output file
$GON_PHYLING/modules/parsnp_taxa_count.sh > taxa_count.txt

# split each number of taxa into a set that can be iterrated over and processed in parallel
#cat taxa_count.txt | split -a 10 -l $threads
cat taxa_count.txt | split -d -l $threads

sed -i -e 's/.fasta//g' ./parsnp.xmfa

for j in $(ls x*); do
  for i in $(cat $j); do

    $GON_PHYLING/modules/parallel_parsnp_splitter.py parsnp.xmfa $i &

  done
  wait
done

# combine the split parsnp sequence files into a single file
cat parsnp_chunk-*-.fa > combo.fas

# remove the first newline character and turn combo.txt into a fasta file
sed -i '1d' combo.fas

printf "\nREMOVING REFERENCE IDENTIFIER\n"
# remove .ref from any seq names
sed -i -e 's/.ref//g' combo.fas


fi


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

printf "\nAssembly and inference of genomes complete.\n"

if [ $intermediate != "KEEP" ]; then

	printf "\nCleaning up intermediate files\n"
	cd $read_dir/trimmed_reads
	rm *$r1_tail*
	rm *$r2_tail*
	cd ./spades_output
	for i in $(ls -d ./*/); do
        if [ $i != "./genomes_for_parsnp/" ]; then
                #printf "\n$i\n"
		rm -r "$i"
        fi
	done
	cd ./genomes_for_parsnp
	rm -r ./P_*
	rm ./*.fasta
	cd ./alignment_fixing
	rm ./*chunk*
	printf "\nFinished cleaning up intermediate files\n"

fi	
