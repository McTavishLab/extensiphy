#! /bin/bash
# inputs: a sequence alignment and a tree inferred from that alignment.
# inputs: a directory of paired end reads for new taxa to be added to the alignment and corresponding tree.


set -e
set -u
set -o pipefail

# establishes the path to find the phycorder directory
PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# # # keep track of the last executed command
# current_command=$BASH_COMMAND
# trap last_command=$current_command DEBUG
# # # echo an error message before exiting
# trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

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
if [ $(which bwa-mem2 | wc -l) -lt 1 ]
    then
        printf "bwa-mem2 not found. Install and/or add to path\n" >&2
	exit 0
    else
        printf "bwa-mem2 found\n"
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

outdir="EP_output"
threads=0
r1_tail="R1.fq"
r2_tail="R2.fq"
end_setting="PE"
phycorder_runs=2
align_type="CONCAT_MSA"
output_type="CONCAT_MSA"
single_locus_suffix=".fasta"
loci_positions="loci_positions.csv"
loci_len="700"
bootstrapping="OFF"
tree="NONE"
ref_select="RANDOM"
intermediate="KEEP"
use="ALIGN"
outputs_dir="RESULTS"
storage="intermediate_files"

while getopts ":a:t:o:c:p:e:1:2:m:d:g:s:f:n:b:r:i:u:h" opt; do
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
    n) loci_len="$OPTARG"
    ;;
    b) bootstrapping="$OPTARG"
    ;;
    r) ref_select="$OPTARG"
    ;;
    i) intermediate="$OPTARG"
    ;;
    u) use="$OPTARG"
    ;;
    h) printf  " Extensiphy is a program for quickly adding genomic sequence data to multiple sequence alignments and phylogenies. \n \
    View the README for more specific information. \n \
    Inputs are generally a multiple sequence file in fasta format and a directory of \n \
    Fastq paired-end read sequences. \
    \n\n\n EXAMPLE COMMAND: \
    \n\n /path/to/extensiphy.sh -u ALIGN -a /path/to/alignment_file -d /path/to/directory_of_reads [any other options] \
    \n\n REQUIRED FLAGS \
    \n (-a) alignment in fasta format, \
    \n (-d) directory of paired end fastq read files for all query taxa, \
    \n (-u) produce only an updated alignment or perform full phylogenetic estimation (ALIGN or PHYLO) (DEFAULT: ALIGN)
    \n\n OPTIONAL FLAGS \
    \n (-t) tree in Newick format produced from the input alignment that you wish to update with new sequences or specify NONE to perform new inference (DEFAULT: NONE), \
    \n (-m) alignment type (SINGLE_LOCUS_FILES, PARSNP_XMFA or CONCAT_MSA) (DEFAULT: CONCAT_MSA), \
    \n (-o) directory name to hold results (DEFAULT: creates EP_output), \
    \n (-i) clean up intermediate output files to save HD space (Options: CLEAN, KEEP)(DEFAULT: KEEP), \
    \n (-r) Select a reference sequence from the alignment file for read mapping or leave as default and the first sequence in the alignment used (DEFAULT: RANDOM), \
    \n (-p) number of taxa to process in parallel, \
    \n (-c) number of threads per taxon being processed, \
    \n (-e) set read-type as single end (SE) or pair-end (PE) (DEFAULT: PE) \
    \n (-1, -2) suffix (ex: R1.fastq or R2.fastq) for both sets of paired end files (DEFAULTS: R1.fq and R2.fq), \
    \n (-g) output format (CONCAT_MSA or SINGLE_LOCUS_FILES) (DEFAULT: CONCAT_MSA), \
    \n (-s) specify the suffix (.fa, .fasta, etc) (DEFAULT: .fasta), \
    \n (-b) bootstrapping tree ON or OFF (DEFAULT: OFF) \
    \n\n\n if using single locus MSA files as input, \
    \n (-f) csv file name to keep track of individual loci when concatenated (DEFAULT: loci_positions.csv), \
    \n (-n) Set size of locus minimum size cutoff used as input or output (Options: int number)(DEFAULT: 700) \
    \n"
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

if [ ${use} == "PHYLO" ]; then
# If user is trying to also build a phylogeny, make sure the RAxML is installed.
        if [ $(which raxmlHPC-PTHREADS | wc -l) -lt 1 ]; then
                printf "raxmlHPC not found. Install and/or alias or add to path\n" >&2
	        exit 0
        else
                printf "raxmlHPC found\n"
        fi

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

if [ $align_type == "CONCAT_MSA" ]; then
        if [ ! -f $align ]; then
	        printf "\nAlignment file doesn't exist or pathing is incorrect.\n"
	        exit
        fi

elif [ $align_type == "SINGLE_LOCUS_FILES" ]; then
        if [ ! -d $align ]; then
                printf "\nAlignment file doesn't exist or pathing is incorrect.\n"
	        exit
        fi
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

#test vcffixer module
# $PHYCORDER/tests/misc/module_tests/vcffixer_test.py --ep_path ${PHYCORDER}

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

#CHECK ALIGNMENT IN FASTA FORMAT
#THIS SHOULD BE EXPANDED
if [ $align_type == "CONCAT_MSA" ]; then
        if head -1 $align | grep -q "^>"; then
	        :
        else
	        printf "\nAlignment file $align doesn't appear to be a fasta format file. Please input a fasta file.\n"
	        exit
        fi
fi

#CHECK TREE FILE FOR CORRECT NEWICK FORMAT
if [ $tree != "NONE" ]; then
	if grep -q "^(.*:.*):" $tree ; then
		:
	else
		printf "\nTree file $tree appears to not be in newick format.\n"
		exit
	fi
fi

#echo "$r1_tail"
#echo "$r2_tail"
#echo "$read_dir"

if [ "${end_setting}" == "PE" ]; then
  # echo $( ls $read_dir/*$r1_tail )

  myarray_r1=(`find  $read_dir -maxdepth 1 -name "*$r1_tail"`)
  myarray_r2=(`find  $read_dir -maxdepth 1 -name "*$r2_tail"`)

  count_r1=${#myarray_r1[@]}
  count_r2=${#myarray_r2[@]}

  echo "$count_r1 files found in $read_dir with suffix $r1_tail"
  echo "$count_r2 files found in $read_dir with suffix $r2_tail"

  if  [ $count_r1 == 0 ] ; then echo "No files found in $read_dir with suffix $r1_tail" && exit 1;  fi
  if  [ $count_r2 == 0 ]; then echo "No files found in $read_dir with suffix $r2_tail" && exit 1;  fi

elif [ "${end_setting}" == "SE" ]; then
  echo $( ls $read_dir/*$r1_tail )

  myarray_r1=(`find  $read_dir -maxdepth 1 -name "*$r1_tail"`)

  count_r1=${#myarray_r1[@]}

  echo "$count_r1 files found in $read_dir with suffix $r1_tail"

  if  [ $count_r1 == 0 ] ; then echo "No files found in $read_dir with suffix $r1_tail" && exit 1;  fi

fi

# CHECK READ FILES SUFFIX
#printf "\nBeginning check of read and suffix accuracy\n"
#if [ "${end_setting}" == "PE" ]; then
#	if ls ${read_dir} | grep -q "${r1_tail}" || ls ${read_dir} | grep -q "${r2_tail}"
#	then
#		for read_1 in $( ls -1 $read_dir/*$r1_tail ); do
#			if head -1 "$read_1" | grep -q "^@" && head -3 "$read_1" | tail -1 | grep -q "^+"; then
#				:
#			else
#				printf "\nRead file $read_1 isn't formatted as a fastq file. Check your read file format before proceeding\n"
#				exit
#			fi
#		done
#		for read_2 in $( ls -1 $read_dir/*$r2_tail ); do
#                        if head -1 "$read_2" | grep -q "^@" && head -3 "$read_2" | tail -1 | grep -q "^+"; then
#                                :
#                        else
#                                printf "\nRead file $read_2 isn't formatted as a fastq file. Check your read file format before proceeding\n"
#                                exit
#                        fi
#                done
#
#		:
#		#printf "\nRead suffixes for paired end reads found in specified read directory. Continuing with analysis.\n"
#	else
#		printf "\nRead suffixes for paired end reads not found in specified read directory. Check your read suffixes and try again.\n"
#		exit
#	fi
#
#elif [ "$end_setting" == "SE" ]; then
#	if ls $read_dir | grep -q "$r1_tail"
#		then
#		:
#                #printf "\nRead suffixes for single end reads found in specified read directory. Continuing with analysis.\n"
#        else
#                printf "\nRead suffixes for single end reads not found in specified read directory. Check your read suffixes and try again.\n"
#                exit
#        fi
#
#fi


###########################################################
printf "\n###################################################\n"
printf "alignment type = $align_type\n"
printf "alignment file = $align\n"
printf "tree file = $tree\n"
printf "directory of reads = $read_dir\n"
printf "reference selection = $ref_select\n"
printf "number of Extensiphy runs = $phycorder_runs\n"
printf "number of threads per Extensiphy run = $threads\n"
printf "suffix for left reads (if paired end or single end) = $r1_tail\n"
if [ ${end_setting} == "PE" ]; then
  printf "suffix for right reads (if paired end only) = $r2_tail\n"
fi
printf "output directory = $outdir\n"
printf "#################################################\n"


#if [ -d $outdir ]; then
#	printf "Output folder exists. Choose a different name.\n"
#	exit
#fi

echo "${outdir}"
mkdir -p ${outdir}

tmp_outdir=$(realpath ${outdir})
outdir=${tmp_outdir}
echo "${outdir}"
cd ${outdir}

mkdir ${outdir}/${outputs_dir}
mkdir ${outdir}/${storage}

workd=$(pwd)

touch ${workd}/ep_dev_log.txt

if [ ${align_type} == "PARSNP_XMFA" ]; then

	mkdir locus_msa_files

	cat $align | grep -Po "cluster\d+" | sort | uniq > ./locus_msa_files/locus_IDs.txt

        cd ./locus_msa_files

        #echo pwd
        cat ./locus_IDs.txt | split -d -l $threads

        if [ $loci_len == "700" ]; then

                for j in $(ls x*); do
                        for i in $(cat $j); do
                    $PHYCORDER/modules/locus_splitter.py --align_file $align --out_file ./$i-.fasta --locus_id $i --locus_size 700 >> $workd/ep_dev_log.txt 2>&1
                        done
                        wait
                done


            $PHYCORDER/modules/new_locus_combiner.py --msa_folder ./ --suffix .fasta --out_file ../combo.fas --position_csv_file $workd/$loci_positions --suffix $single_locus_suffix --len_filter 700 >> $workd/ep_dev_log.txt 2>&1

        elif [ $loci_len != "700" ]; then
                for j in $(ls x*); do
                        for i in $(cat $j); do
                    $PHYCORDER/modules/locus_splitter.py --align_file $align --out_file ./$i-.fasta --locus_id $i --locus_size $loci_len >> $workd/ep_dev_log.txt 2>&1
                        done
                        wait
                done


            $PHYCORDER/modules/new_locus_combiner.py --msa_folder ./ --suffix .fasta --out_file ../combo.fas --position_csv_file $workd/$loci_positions --suffix $single_locus_suffix --len_filter $loci_len >> $workd/ep_dev_log.txt 2>&1

        fi

        align=$( realpath ../combo.fas)

    printf "\nNew alignment file produced\n"
    printf "$align"

    cd ..

elif [ $align_type == "SINGLE_LOCUS_FILES" ]; then
        if [ $loci_len == "700" ]; then
                printf "\nFiltering and combining input single locus alignments by default length of $loci_len\n"
	        $PHYCORDER/modules/new_locus_combiner.py --msa_folder $align --suffix $single_locus_suffix --out_file $workd/combo.fas --position_csv_file $loci_positions --len_filter 700 >> $workd/ep_dev_log.txt 2>&1
	        #printf "$outdir\n"
	        #printf "$PHYCORDER\n"

	        align=$( realpath $workd/combo.fas )
	        loci_positions=$( realpath $loci_positions)

        elif [ $loci_len != "700" ]; then
                printf "\nFiltering and combining input single locus alignments by user specified length of $loci_len\n"
                $PHYCORDER/modules/new_locus_combiner.py --msa_folder $align --suffix $single_locus_suffix --out_file $workd/combo.fas --position_csv_file $loci_positions --len_filter $loci_len >> $workd/ep_dev_log.txt 2>&1
	        #printf "$outdir\n"
	        #printf "$PHYCORDER\n"

	        align=$( realpath $workd/combo.fas )
	        loci_positions=$( realpath $loci_positions)
        fi

elif [ $align_type == "CONCAT_MSA" ]; then

	printf "Concatenated Multiple Sequence Alignment selected as input.\n"
	printf "Assuming loci are all longer than 700bp. Continuing with rapid updating\n"

fi


###########################

# Strip newline characters from sequences so each sequence is contiguous
# Also strips gap characters ('-') from each sequence
$PHYCORDER/modules/degen_fixer.py --align_file $align --output $workd/ref_nogap.fas


if [ $ref_select != "RANDOM" ]; then

        printf "\nGenerating a reference from taxon $ref_select\n"
        $PHYCORDER/modules/ref_producer.py -s --align_file $align --ref_select $ref_select --out_file $workd/best_ref_gaps.fas >> $workd/ep_dev_log.txt 2>&1

        $PHYCORDER/modules/ref_producer.py -s --align_file $workd/ref_nogap.fas --ref_select $ref_select --out_file $workd/best_ref.fas >> $workd/ep_dev_log.txt 2>&1


elif [ $ref_select == "RANDOM" ]; then
	# PRODUCE SINGLE REFERENCE SEQUENCES, BOTH WITH AND WITHOUT GAPS FOR ALIGNMENT AND EVENTUAL INCLUSION OF NEW SEQUENCES INTO ALIGNMENT

        printf "\nGenerating a reference from first sequence in alignment\n"
	$PHYCORDER/modules/ref_producer.py -r --align_file $align --out_file $workd/best_ref_gaps.fas >> $workd/ep_dev_log.txt 2>&1

	$PHYCORDER/modules/ref_producer.py -r --align_file $workd/ref_nogap.fas --out_file $workd/best_ref.fas >> $workd/ep_dev_log.txt 2>&1

fi

	# BUILD REFERENCE ALIGNMENT LIBRARY
#hisat2-build --threads $threads $workd/best_ref.fas $workd/best_ref >> $workd/rapup_dev_log.txt 2>&1
bwa-mem2 index $workd/best_ref.fas >> $workd/ep_dev_log.txt 2>&1

# ls ${read_dir}/*$r1_tail | split -a 5 -l $phycorder_runs

#ls ${read_dir}/*$r1_tail | split -a 10 -l $phycorder_runs

# SPLIT UP READ FILES NAMES INTO SEPARATE FILES THAT CONSTITUTE RAPUP RUNS
# X01 IS A SINGLE RUN, X02, etc
ls ${read_dir}/*$r1_tail | split -d -l $phycorder_runs

#printf "\nNumber of cores allocated enough to process all read sets\n"
printf "\nBeginning Extensiphy runs\n"

#TODO: ADD MORE FUNCTIONS TO THIS PARALLEL RUNNING SECTION IF POSSIBLE
# SUCH AS REMOVING FILES IF CLEAN OPTION IS SPECIFIED

wd=$(pwd)
mkdir -p combine_and_infer

echo "end setting == ${end_setting}"

if [ "$end_setting" == "PE" ]; then
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
    			$PHYCORDER/modules/map_to_align.sh -a $workd/best_ref.fas -t $tree -p $i -e ${i%$r1_tail}$r2_tail -1 $r1_tail -2 $r2_tail -c $threads -d "$workd" -g $workd/best_ref_gaps.fas -o ${base}output_dir >> $workd/ep_dev_log.txt 2>&1 &
    			#> rapup-dev-logs/parallel-$base-dev.log 2> rapup-dev-logs/parallel-$base-dev-err.log &
    			#wait
    			printf "\nassembling new taxon\n"
		done
		wait

                for i in $(ls -d *output_dir); do
                        cd $i
                        if [ $(ls -l | grep "_align.fas" | wc -l) -gt 0 ]; then
                                cp *_align.fas $wd/combine_and_infer/
                        elif [ $(ls -l | grep "_align.fas" | wc -l) -gt 0 ]; then
                                echo "$i NO FASTA PRODUCED"
                                continue
                        fi
                        cd ..
                done

                if [ ${intermediate} == "CLEAN" ]; then
                        printf "\nCleaning up intermediate output files.\n"
	                for i in $(ls -d *output_dir); do
		                cd $i
		                for j in $(ls -1); do
			                rm ./$j
		                done
		                cd ..
		                rmdir $i
	                done
                fi

	done

elif [ "$end_setting" == "SE" ]; then
	for j in $(ls x*); do
                for i in $(cat $j); do
                        base=$(basename $i $r1_tail)
                        #echo $base
                        #echo $imv ${outdir}/combine_and_infer/extended.aln ${outdir}/outputs/
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
                        $PHYCORDER/modules/map_to_align.sh -a $workd/best_ref.fas -t $tree -p $i -1 $r1_tail -c $threads -d ${workd} -g $workd/best_ref_gaps.fas -o ${base}output_dir >> $workd/ep_dev_log.txt 2>&1 &
                        #> rapup-dev-logs/parallel-$base-dev.log 2> rapup-dev-logs/parallel-$base-dev-err.log &
                        #wait
                        printf "\nassembling new taxon\n"
                done
                wait

                for i in $(ls -d *output_dir); do
                        cd $i
                        if [ $(ls -l | grep "_align.fas" | wc -l) -gt 0 ]; then
                                cp *_align.fas $wd/combine_and_infer/
                        elif [ $(ls -l | grep "_align.fas" | wc -l) -gt 0 ]; then
                                echo "$i NO FASTA PRODUCED"
                                continue
                        fi
                        cd ..
                        done

                if [ ${intermediate} == "CLEAN" ]; then
                        printf "\nCleaning up intermediate output files.\n"
	                for i in $(ls -d *output_dir); do
		                cd $i
		                for j in $(ls -1); do
			                rm ./$j
		                done
		                cd ..
		                rmdir $i
	                done
                fi
        done
fi


printf "\nIndividual Extensiphy runs finished. Combining aligned query sequences and adding them to starting alignment\n"


if [ $intermediate == "KEEP" ]; then
	printf "\nKeeping intermediate output files\n"
fi

printf "\nskipping renaming step\n"
cat combine_and_infer/*.fas $align > combine_and_infer/extended.aln
# fi

printf "\nExtended alignment file created (extended.aln).\n"

cd ${outdir}/${outputs_dir}

# strip the unnecessary information from the taxa names in the alignment.
# this assumes you've used the renaming tool to rename all of the reads for this experiment
# COMMENTED OUT DUE TO CHANGES TO HOW NAME HANDLING IS WORK MOVING FORWARD
#sed -i -e 's/_[^_]*//2g' extended.aln

INFER=$(pwd)

# move files into appropriate places for later storage or user use
mv ${outdir}/combine_and_infer/extended.aln ${outdir}/${outputs_dir}/
mv ${outdir}/best_ref* ${outdir}/${storage}/
mv ${outdir}/ref_nogap.fas ${outdir}/${storage}/
mv ${outdir}/snps.txt ${outdir}/${storage}/
mv ${outdir}/x* ${outdir}/${storage}/


if [ $use == "ALIGN" ]; then
    printf "\nAlignment updating complete.\n"
    printf "\nAlignment file will be in ${outdir}/RESULTS \n"

elif [ $use == "PHYLO" ]; then
    printf "\nAlignment updating complete.\n"
    printf "\nPhylogenetic estimation selected.\n"
    printf "\n"
     # handling of bootstrapping

    #printf "\nAlignment updating complete. Moving to phylogenetic inference.\n"
    if [ $bootstrapping == "ON" ]; then

    # handles whether user wants to use a starting tree or not
      if [ $tree == "NONE" ]; then
        raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -p 12345 -n consensusFULL >> $workd/ep_dev_log.txt 2>&1

        printf "\nInference of updated phylogenetic tree complete\n"
        printf "\nRunning 100 bootstrap replicates\n"
        raxmlHPC-PTHREADS -s extended.aln -n consensusFULL_bootstrap -m GTRGAMMA  -p 12345 -T $threads -N 100 -x 12345 >> $workd/ep_dev_log.txt 2>&1


        raxmlHPC-PTHREADS -z RAxML_bootstrap.consensusFULL_bootstrap -t RAxML_bestTree.consensusFULL -f b -T $threads -m GTRGAMMA -n majority_rule_bootstrap_consensus >> $workd/ep_dev_log.txt 2>&1

        printf "\nML tree, bootstrap replicates in $outdir/RESULTS.\n"

      elif [ $tree != "NONE" ]; then

       raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -t $tree -p 12345 -n consensusFULL >> $workd/ep_dev_log.txt 2>&1


       printf "\nMultiple taxa update of input tree $tree complete\n"

       printf "\nRunning 100 bootstrap replicates\n"
       raxmlHPC-PTHREADS -s extended.aln -n consensusFULL_bootstrap -m GTRGAMMA  -p 12345 -T $threads -N 100 -x 12345 >> $workd/ep_dev_log.txt 2>&1

       raxmlHPC-PTHREADS -z RAxML_bootstrap.consensusFULL_bootstrap -t RAxML_bestTree.consensusFULL -f b -T $threads -m GTRGAMMA -n majority_rule_bootstrap_consensus >> $workd/ep_dev_log.txt 2>&1

       printf "\nML tree, bootstrap replicates in $outdir/RESULTS.\n"


     fi

    elif [ $bootstrapping == "OFF" ]; then

      # handles whether user wants to use a starting tree or not
      if [ $tree == "NONE" ]; then

        raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -p 12345 -n consensusFULL >> $workd/ep_dev_log.txt 2>&1

        printf "\nInference of updated phylogenetic tree complete\n"
      elif [ $tree != "NONE" ]; then

        raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -t $tree -p 12345 -n consensusFULL >> $workd/ep_dev_log.txt 2>&1

        printf "\nMultiple taxa update of phylogenetic tree $tree complete\n"

      fi

    else
      printf "\nSwitch bootstrapping option to 'ON' or 'OFF' and re-run program.\n"

    fi
#end of phylo inference "if" statement
fi

output_dir=$(pwd)

# handling of multiple single locus MSA files as input
# this will be expanded as HGT detection is added
# for now, it serves as a SNP check
if [ $output_type == "SINGLE_LOCUS_FILES" ]; then

	$PHYCORDER/modules/locus_position_identifier.py --out_file_dir ${INFER}/updated_single_loci --position_csv_file $loci_positions --concatenated_fasta $INFER/extended.aln >> $workd/ep_dev_log.txt 2>&1

	printf "\nMultiple single locus MSA file handling selected\n"
	printf "\nAlignment file is: "$output_dir/"extended.aln\n"
        if [ $use == "PHYLO" ]; then
	    printf "\nTree file is: "$output_dir/"RAxML_bestTree.consensusFULL\n"
        fi
	printf "\nSingle locus alignment files are in: "$output_dir/"updated_single_loci\n"

elif [ $output_type == "CONCAT_MSA" ]; then
	printf "\nSingle concatenated loci MSA file handling selected\n"
	printf "\nAlignment file is: "$output_dir/"extended.aln\n"
        if [ $use == "PHYLO" ]; then
            printf "\nTree file is: "$output_dir/"RAxML_bestTree.consensusFULL\n"
        fi


fi

   # printf "Moving run logs into phycorder-dev-logs"
   # cd ..
   #
   # mv *-dev.log $outdir/phycorder-dev-logs
