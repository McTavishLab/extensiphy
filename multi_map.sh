#! /bin/bash
# inputs: a sequence alignment and a tree inferred from that alignment.
# inputs: a directory of paired end reads for new taxa to be added to the alignment and corresponding tree.
# example command: ./multi_map.sh ./sample_phycorder.cfg

set -e
set -u
set -o pipefail

# establishes the path to find the phycorder directory
PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# changing location of .cfg file to a variable
# for easy use of multiple config num_files
#source $1


############################################################

#Check for dependencies
if [ $(which bcftools | wc -l) -lt 1 ]
    then
        printf "Requires bcftools" >&2
        exit 0
    else
        printf "Correct version of bfctools found.\n"
fi
if [  $(which samtools | wc -l) -lt 1 ] #TODO steup for greater than 1.2? this  is a sloppppy approach
    then
        printf "Requires samtools" >&2
        exit 0
    else
        printf "Correct version of samtools found.\n"
fi
if [ $(which seqtk | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "seqtk not found\n" >&2
	exit 0
    else
        printf "seqtk found\n"
fi
if [ $(which hisat2 | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "hisat2 not found. Install and/or add to path\n" >&2
	exit 0
    else
        printf "hisat2 found\n"
fi
if [ $(which raxmlHPC-PTHREADS | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "raxmlHPC not found. Install and/or alias or add to path\n" >&2
	exit 0
    else
        printf "raxmlHPC found\n"
fi
if [ $(which fastx_collapser | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "fastx toolkit not found. Install and/or add to path\n" >&2
	exit 0
    else
        printf "fastx toolkit found\n"
fi
if [ $(which vcfutils.pl | wc -l) -lt 1 ] #TODO needs different install than bcftools?
    then
        printf "vcfutils.pl not found. Install and/or add to path\n" >&2
	exit 0
    else
        printf "vcfutils.pl found\n"
fi

printf "\n\n"

outdir="rapup_run"
threads=0
r1_tail="R1.fq"
r2_tail="R2.fq"
end_setting="PE"
phycorder_runs=2
align_type="CONCAT_MSA"
output_type="CONCAT_MSA"
single_locus_suffix=".fasta"
loci_positions="loci_positions.csv"
bootstrapping="OFF"
tree="NONE"
ref_select="RANDOM"

while getopts ":a:t:o:c:p:e:1:2:m:d:g:s:f:b:r:h" opt; do
  case $opt in
    a) align="$OPTARG"
    ;;
    t) tree="$OPTARG"
    ;;
    o) outdir="$OPTARG"
    ;;
    c) threads="$OPTARG"
    ;;
    p) phycorder_runs="$OPTARG"
    ;;
    e) end_setting="$OPTARG"
    ;;
    1) r1_tail="$OPTARG"
    ;;
    2) r2_tail="$OPTARG"
    ;;
    m) align_type="$OPTARG"
    ;;
    d) read_dir="$OPTARG"
    ;;
    g) output_type="$OPTARG"
    ;;
    s) single_locus_suffix="$OPTARG"
    ;;
    f) loci_positions="$OPTARG"
    ;;
    b) bootstrapping="$OPTARG"
    ;;
    r) ref_select="$OPTARG"
    ;;
    h) printf  " RapUp is a program for quickly adding genomic sequence data to multiple sequence alignments and phylogenies. View the README for more specific information. Inputs are generally a multiple sequence file in .fasta format and a directory of .fastq paired-end read sequences.\n\n\n EXAMPLE COMMAND:\n\n /path/to/multi_map.sh -a /path/to/alignment_file -d /path/to/directory_of_reads [any other options]\n\n (-a) alignment in fasta format,\n (-d) directory of paired end fastq read files for all query taxa,\n (-t) tree in Newick format produced from the input alignment that you wish to update with new sequences or specify NONE to perform new inference (DEFAULT: NONE),\n (-m) alignment type (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA) (DEFAULT: CONCAT_MSA),\n (-o) directory name to hold results (DEFAULT: creates rapup_run),\n (-r) Selected a reference sequence from the alignment file for read mapping or leave as default and a random reference will be chosen (DEFAULT: RANDOM),\n (-p) number of taxa to process in parallel,\n (-c) number of threads per taxon being processed,\n (-e) set read-type as single end (SE) or pair-end (PE) (DEFAULT: PE)\n (-1, -2) suffix (ex: R1.fastq or R2.fastq) for both sets of paired end files (DEFAULTS: R1.fq and R2.fq),\n (-g) output format (CONCAT_MSA or SINGLE_LOCUS_FILES) (DEFAULT: CONCAT_MSA),\n (-s) specify the suffix (.fa, .fasta, etc) (DEFAULT: .fasta),\n (-b) bootstrapping tree ON or OFF (DEFAULT: OFF)\n\n\n if using single locus MSA files as input,\n (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv),\n"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z "$align" ]; then
   "alignment file (-a) required"
   exit
fi

#Ttest if files actually exist
#Check to make sure mapping has occured if re-mapping

#if [ -f "$align" ]; then
#    printf "Alignment is %s\n" "$align"
#  else
#    printf "Alignment $align not found. Exiting\n" >&2
#    exit
#fi


if [ $threads -eq 0 ]; then
     threads=2
     #printf "Threads not set, defaulting to 2" >&2
else
     #echo "num threads is"
     #echo $threads
     :
fi

if [ -d $outdir ]; then
        printf "Output folder exists. Choose a different name.\n"
        exit
fi

if [ ! -f $align ]; then
	printf "\nAlignment file doesn't exist or pathing is incorrect.\n"
	exit
fi

if [ $tree != "NONE" ]; then
	if [ ! -f "$tree" ]; then
		printf "\nTree file doesn't exist or pathing is incorrect.\n"
		exit
	fi
fi

if [ ! -d $read_dir ]; then
        printf "Directory of read sequences can't be found. Check name and pathing.\n"
        exit
fi

tmp_align=$(realpath $align)
tmp_tree=$(realpath $tree)
tmp_read_dir=$(realpath $read_dir)

align=$tmp_align
#tree=$tmp_tree
read_dir=$tmp_read_dir

if [ $tree == "NONE" ]; then
	tmp_tree="NONE"
elif [ $tree != "NONE" ]; then
	tmp_tree=$(realpath $tree)
fi

tree=$tmp_tree


# CHECK READ FILES SUFFIX
printf "\nBeginning check of read and suffix accuracy\n"
if [ "$end_setting" == "PE" ]; then
	if ls $read_dir | grep -q "$r1_tail" || ls $read_dir | grep -q "$r2_tail"
	then
		:
		#printf "\nRead suffixes for paired end reads found in specified read directory. Continuing with analysis.\n"
	else
		printf "\nRead suffixes for paired end reads not found in specified read directory. Check your read suffixes and try again.\n"
		exit
	fi

elif [ "$end_setting" == "SE" ]; then
	if ls $read_dir | grep -q "$r1_tail"
		then
		:
                #printf "\nRead suffixes for single end reads found in specified read directory. Continuing with analysis.\n"
        else
                printf "\nRead suffixes for single end reads not found in specified read directory. Check your read suffixes and try again.\n"
                exit
        fi
	
fi


###########################################################
printf "\n###################################################\n"
printf "alignment type = $align_type\n"
printf "alignment file = $align\n"
printf "tree file = $tree\n"
printf "directory of reads = $read_dir\n"
printf "reference selection = $ref_select\n"
printf "number of RapUp runs = $phycorder_runs\n"
printf "number of threads per RapUp run = $threads\n"
printf "suffix for left reads (if paired end or single end) = $r1_tail\n"
printf "suffix for right reads (if paired end only) = $r2_tail\n"
printf "output directory = $outdir\n"
printf "#################################################\n"


#if [ -d $outdir ]; then
#	printf "Output folder exists. Choose a different name.\n"
#	exit       
#fi


mkdir -p $outdir

cd $outdir

workd=$(pwd)

touch $workd/rapup_dev_log.txt

if [ $align_type == "PARSNP_XMFA" ]; then

	mkdir locus_msa_files

	cat $align | grep -Po "cluster\d+" | sort | uniq > ./locus_msa_files/locus_IDs.txt

        cd ./locus_msa_files

        #echo pwd
        cat ./locus_IDs.txt | split -d -l $threads



        for j in $(ls x*); do
                for i in $(cat $j); do
			$PHYCORDER/locus_splitter.py --align_file $align --out_file ./$i-.fasta --locus_id $i --locus_size 1000 >> $workd/rapup_dev_log.txt 2>&1
                done
                wait
        done


	$PHYCORDER/new_locus_combiner.py --msa_folder ./ --suffix .fasta --out_file ../combo.fas --position_csv_file $workd/$loci_positions --suffix $single_locus_suffix --len_filter 1000 >> $workd/rapup_dev_log.txt 2>&1
	align=$( realpath ../combo.fas)

	printf "\nNew alignment file produced\n"
	printf "$align"

	cd ..

elif [ $align_type == "SINGLE_LOCUS_FILES" ]; then
	$PHYCORDER/new_locus_combiner.py --msa_folder $align --suffix $single_locus_suffix --out_file $workd/combo.fas --position_csv_file $loci_positions --len_filter 1000 >> $workd/rapup_dev_log.txt 2>&1
	#printf "$outdir\n"
	#printf "$PHYCORDER\n"	

	align=$( realpath $workd/combo.fas )
	loci_positions=$( realpath $loci_positions)

elif [ $align_type == "CONCAT_MSA" ]; then

	printf "Concatenated Multiple Sequence Alignment selected as input. Assuming that documentation has been read\n"
	printf "Assuming loci are all longer than 1000bp. Continuing with rapid updating\n"

fi


###########################

if [[ ! -z $(grep "-" $align) ]]; then
  printf "\nGAP FOUND BEFORE REMOVAL\n"
else
  printf "\nNO GAPS FOUND BEFORE REMOVAL\n";
fi

#printf "sed 's/-//g' <$align > $new_out/ref_nogap.fas"

#printf "current wd\n"
#pwd
#pull all the gaps from the aligned taxa bc mappers cannot cope.
sed 's/-//g' <$align > $workd/ref_nogap.fas


if [[ ! -z $(grep "-" ./ref_nogap.fas) ]]; then
  printf "\nGAP FOUND AFTER REMOVAL!\n"
else
  printf "\nNO GAPS FOUND AFTER REMOVAL\n";
fi

if [ $ref_select != "RANDOM" ]; then

        $PHYCORDER/ref_producer.py -s --align_file $align --ref_select $ref_select --out_file $workd/best_ref_gaps.fas >> $workd/rapup_dev_log.txt 2>&1

        $PHYCORDER/ref_producer.py -s --align_file $workd/ref_nogap.fas --ref_select $ref_select --out_file $workd/best_ref.fas >> $workd/rapup_dev_log.txt 2>&1


elif [ $ref_select == "RANDOM" ]; then
	# PRODUCE SINGLE REFERENCE SEQUENCES, BOTH WITH AND WITHOUT GAPS FOR ALIGNMENT AND EVENTUAL INCLUSION OF NEW SEQUENCES INTO ALIGNMENT
	$PHYCORDER/ref_producer.py -r --align_file $align --out_file $workd/best_ref_gaps.fas >> $workd/rapup_dev_log.txt 2>&1

	$PHYCORDER/ref_producer.py -r --align_file $workd/ref_nogap.fas --out_file $workd/best_ref.fas >> $workd/rapup_dev_log.txt 2>&1

fi

	# BUILD REFERENCE ALIGNMENT LIBRARY
hisat2-build --threads $threads $workd/best_ref.fas $workd/best_ref >> $workd/rapup_dev_log.txt 2>&1


# ls ${read_dir}/*$r1_tail | split -a 5 -l $phycorder_runs

#ls ${read_dir}/*$r1_tail | split -a 10 -l $phycorder_runs

# SPLIT UP READ FILES NAMES INTO SEPARATE FILES THAT CONSTITUTE RAPUP RUNS
# X01 IS A SINGLE RUN, X02, etc
ls ${read_dir}/*$r1_tail | split -d -l $phycorder_runs

printf "\nNumber of cores allocated enough to process all read sets\n"
printf "\nBeginning RapUp runs\n"


for j in $(ls x*); do
for i in $(cat $j); do
    base=$(basename $i $r1_tail)
    #echo $base
    #echo $i
    #echo $PHYCORDER
    #echo "$workd ####################################################################################################\n"
    #echo $align
    #echo $tree x
    #echo $i
    #echo ${base}${r2_tail}
    #echo $threads
    #echo $align_type
    #echo "${base}_output_dir"
    #echo "$PHYCORDER/map_to_align.sh -a $outdir/best_ref.fas -t $tree -p $i -e ${i%$r1_tail}$r2_tail -1 $r1_tail -2 $r2_tail -c $threads -o ${base}output_dir > parallel-$base-dev.log &"
    #echo "Time for $j Phycorder run:"
    time $PHYCORDER/map_to_align.sh -a $workd/best_ref.fas -t $tree -p $i -e ${i%$r1_tail}$r2_tail -1 $r1_tail -2 $r2_tail -c $threads -d "$workd" -g $workd/best_ref_gaps.fas -o ${base}output_dir >> $workd/rapup_dev_log.txt 2>&1 & 
    #> rapup-dev-logs/parallel-$base-dev.log 2> rapup-dev-logs/parallel-$base-dev-err.log &
    #wait
    printf "\nadding new map_to_align run\n"
done
wait
done

printf "\nIndividual Phycorder runs finished. Combining aligned query sequences and adding them to starting alignment\n"

mkdir -p combine_and_infer


wd=$(pwd)

# loop through phycorder run directories and move finished fasta files to /combine_and_infer/
# for tree inference
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

printf "\nskipping renaming step\n"
cat combine_and_infer/*.fas $align > combine_and_infer/extended.aln
# fi

printf "\nExtended alignment file creaded (extended.aln), using previous tree as starting tree for phylogenetic inference\n"

cd combine_and_infer

# strip the unnecessary information from the taxa names in the alignment.
# this assumes you've used the renaming tool to rename all of the reads for this experiment
# COMMENTED OUT DUE TO CHANGES TO HOW NAME HANDLING IS WORK MOVING FORWARD
#sed -i -e 's/_[^_]*//2g' extended.aln

INFER=$(pwd)

 # handling of bootstrapping

printf "\nAlignment updating complete. Moving to phylogenetic inference.\n"
if [ $bootstrapping == "ON" ]; then

# handles whether user wants to use a starting tree or not
  if [ $tree == "NONE" ]; then
    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -p 12345 -n consensusFULL >> $workd/rapup_dev_log.txt 2>&1

    time raxmlHPC-PTHREADS -s extended.aln -n consensusFULL_bootstrap -m GTRGAMMA  -p 12345 -T $threads -N 100 -b 12345 >> $workd/rapup_dev_log.txt 2>&1

    time raxmlHPC-PTHREADS -z RAxML_bootstrap.consensusFULL_bootstrap -t RAxML_bestTree.consensusFULL -f b -T $threads -m GTRGAMMA -n majority_rule_bootstrap_consensus >> $workd/rapup_dev_log.txt 2>&1

    printf "\nMultiple taxa update of phylogenetic tree complete\n"

  elif [ $tree != "NONE" ]; then

   time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -t $tree -p 12345 -n consensusFULL >> $workd/rapup_dev_log.txt 2>&1

   time raxmlHPC-PTHREADS -s extended.aln -n consensusFULL_bootstrap -m GTRGAMMA  -p 12345 -T $threads -N 100 -b 12345 >> $workd/rapup_dev_log.txt 2>&1

   time raxmlHPC-PTHREADS -z RAxML_bootstrap.consensusFULL_bootstrap -t RAxML_bestTree.consensusFULL -f b -T $threads -m GTRGAMMA -n majority_rule_bootstrap_consensus >> $workd/rapup_dev_log.txt 2>&1

   printf "\nMultiple taxa update of phylogenetic tree complete\n"

 fi

elif [ $bootstrapping == "OFF" ]; then

  # handles whether user wants to use a starting tree or not
  if [ $tree == "NONE" ]; then

    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -p 12345 -n consensusFULL >> $workd/rapup_dev_log.txt 2>&1

  elif [ $tree != "NONE" ]; then

    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -t $tree -p 12345 -n consensusFULL >> $workd/rapup_dev_log.txt 2>&1

    printf "\nMultiple taxa update of phylogenetic tree complete\n"

  fi

else
  printf "\nSwitch bootstrapping option to 'ON' or 'OFF' and re-run program.\n"

fi

output_dir=$(pwd)

# handling of multiple single locus MSA files as input
# this will be expanded as HGT detection is added
# for now, it serves as a SNP check
if [ $output_type == "SINGLE_LOCUS_FILES" ]; then

	$PHYCORDER/locus_position_identifier.py --out_file_dir $INFER/updated_single_loci --position_csv_file $loci_positions --concatenated_fasta $INFER/extended.aln >> $workd/rapup_dev_log.txt 2>&1

	printf "\nMultiple single locus MSA file handling selected\n"
	printf "\nAlignment file is: "$output_dir/"extended.aln\n"
	printf "\nTree file is: "$output_dir/"RAxML_bestTree.consensusFULL\n"
	printf "\nSingle locus alignment files are in: "$output_dir/"updated_single_loci\n"

elif [ $output_type == "CONCAT_MSA" ]; then
	printf "\nSingle concatenated loci MSA file handling selected\n"
	printf "\nAlignment file is: "$output_dir/"extended.aln\n"
        printf "\nTree file is: "$output_dir/"RAxML_bestTree.consensusFULL\n"



fi

   # printf "Moving run logs into phycorder-dev-logs"
   # cd ..
   #
   # mv *-dev.log $outdir/phycorder-dev-logs
