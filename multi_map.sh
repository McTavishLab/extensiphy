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
source $1

printf "###################################################\n"
printf "$align_type\n"
printf "$align\n"
printf "$tree\n"
printf "$read_dir\n"
printf "$phycorder_runs\n"
printf "$threads\n"
printf "$r1_tail\n"
printf "$r2_tail\n"
printf "$outdir\n"
printf "#################################################\n"


mkdir -p $outdir

cd $outdir

if [ $align_type == "PARSNP_XMFA" ]; then

	mkdir locus_msa_files

	cat $align | grep -Po "cluster\d+" | sort | uniq > ./locus_msa_files/locus_IDs.txt
       	#cat parsnp.xmfa | grep -Po "cluster\d+" | sort | uniq > ./locus_msa_files/locus_IDs.txt

        cd ./locus_msa_files

        echo pwd
        cat ./locus_IDs.txt | split -d -l $threads



        for j in $(ls x*); do
                for i in $(cat $j); do
			$PHYCORDER/locus_splitter.py --align_file $align --out_file ./$i-.fasta --locus_id $i --locus_size 1000
                        #$GON_PHYLING/locus_splitter.py --align_file ../parsnp.xmfa --out_file ./$i-.fasta --locus_id $i --locus_size 1000
                done
                wait
        done

        # TODO: ADD COLLECTION SCRIPT THAT ASSEMBLES SINGLE LOCUS FILES INTO CONCATENATED FILE

        $PHYCORDER/locus_combiner.py --msa_folder ./ --suffix .fasta --out_file ../combo.fas --position_dict_file $loci_positions

	align=$( realpath ../combo.fas)

	printf "New alignment file produced\n"
	printf "$align"

	cd ..

elif [ $align_type == "SINGLE_LOCUS_FILES" ]; then
	
	$PHYCORDER/msa_producer.py --align_dir $align --len_filter 1000 --out_file $PHYCORDER/$outdir/combo.fas
	printf "$outdir\n"
	printf "$PHYCORDER\n"	

	align=$( realpath $PHYCORDER/$outdir/combo.fas )


elif [ $align_type == "CONCAT_MSA" ]; then

	printf "Concatenated Multiple Sequence Alignment selected as input. Assuming that documentation has been read\n"
	printf "Assuming loci are all longer than 1000bp. Continuing with rapid updating\n"

fi

# cd $outdir

# ls ${read_dir}/*$r1_tail | split -a 5 -l $phycorder_runs

ls ${read_dir}/*$r1_tail | split -a 10 -l $phycorder_runs


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
    echo $align_type
    echo "${base}_output_dir"
    echo "$PHYCORDER/map_to_align.sh -a $align -t $tree -p $i -e ${i%$r1_tail}$r2_tail -1 $r1_tail -2 $r2_tail -c $threads -o ${base}output_dir > parallel-$base-dev.log &"
    echo "Time for $j Phycorder run:"
    time $PHYCORDER/map_to_align.sh -a $align -t $tree -p $i -e ${i%$r1_tail}$r2_tail -1 $r1_tail -2 $r2_tail -c $threads -o ${base}output_dir > parallel-$base-dev.log &
    #wait
    printf "adding new map_to_align run"
done
wait
done

printf "Individual Phycorder runs finished. Combining aligned query sequences and adding them to starting alignment\n"

mkdir -p combine_and_infer

mkdir -p phycorder-dev-logs

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

printf "skipping renaming step"
cat combine_and_infer/*.fas $align > combine_and_infer/extended.aln
# fi

printf "Extended alignment file creaded (extended.aln), using previous tree as starting tree for phylogenetic inference\n"

cd combine_and_infer

# strip the unnecessary information from the taxa names in the alignment.
# this assumes you've used the renaming tool to rename all of the reads for this experiment

sed -i 's/_$//g' extended.aln

INFER=$(pwd)

 # handling of bootstrapping

printf "Alignment updating complete. Moving to phylogenetic inference."
if [ $bootstrapping == "ON" ]; then

# handles whether user wants to use a starting tree or not
  if [ $tree == "NONE" ]; then
    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -p 12345 -n consensusFULL

    time raxmlHPC-PTHREADS -s extended.aln -n consensusFULL_bootstrap -m GTRGAMMA  -p 12345 -T $threads -N 100 -b 12345

    time raxmlHPC-PTHREADS -z RAxML_bootstrap.consensusFULL_bootstrap -t RAxML_bestTree.consensusFULL -f b -T $threads -m GTRGAMMA -n majority_rule_bootstrap_consensus

    printf "Multiple taxa update of phylogenetic tree complete\n"

  elif [ $tree != "NONE" ]; then

   time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -t $tree -p 12345 -n consensusFULL

   time raxmlHPC-PTHREADS -s extended.aln -n consensusFULL_bootstrap -m GTRGAMMA  -p 12345 -T $threads -N 100 -b 12345

   time raxmlHPC-PTHREADS -z RAxML_bootstrap.consensusFULL_bootstrap -t RAxML_bestTree.consensusFULL -f b -T $threads -m GTRGAMMA -n majority_rule_bootstrap_consensus

   printf "Multiple taxa update of phylogenetic tree complete\n"

 fi

elif [ $bootstrapping == "OFF" ]; then

  # handles whether user wants to use a starting tree or not
  if [ $tree == "NONE" ]; then

    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -p 12345 -n consensusFULL

  elif [ $tree != "NONE" ]; then

    time raxmlHPC-PTHREADS -m GTRGAMMA -T $threads -s $INFER/extended.aln -t $tree -p 12345 -n consensusFULL

    printf "Multiple taxa update of phylogenetic tree complete\n"

  fi

else
  printf "Switch bootstrapping option to 'ON' or 'OFF' and re-run program."

fi


# handling of multiple single locus MSA files as input
# this will be expanded as HGT detection is added
# for now, it serves as a SNP check
if [ $output_type == "SINGLE_LOCUS_FILES" ]; then

	ls $INFER/*.fas | split -d -l $phycorder_runs

	for j in $(ls x*); do
		for i in $(cat $j); do
			$PHYCORDER/locus_position_identifier.py --out_file_dir $INFER --position_dict_file $loci_positions --concatenated_fasta $i 
		done
		wait
	done

	echo "Multiple single locus MSA file handling selected"
elif [ $output_type == "CONCAT_MSA" ]; then
	echo "Single concatenated loci MSA file handling selected"



fi

   # printf "Moving run logs into phycorder-dev-logs"
   # cd ..
   #
   # mv *-dev.log $outdir/phycorder-dev-logs
