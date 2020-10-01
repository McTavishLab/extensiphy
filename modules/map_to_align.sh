#!/bin/bash
#inputs: an alignment, a tree, and whole genome reads
#determine how to try differnet alignments
#deal with possibility of multiple samples?
#HOw to deal with single mapping to multiple alignements??!?
#./map_to_align.sh -a example.aln -t tree.tre -p /home/ejmctavish/projects/Exelixis/SISRS/full_aln/datafiles/SRR610374 -o fulltest -n fulltest #
#DEFAULT ARGS
set -e
set -u
set -o pipefail

PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


printf "phycorder directory is %s\n" "$PHYCORDER"


PE=0
outdir=phycorder_run
nam=QUERY
read_align=0
re_map=1
map=1
read_name_prefix=SRR
wre_map=0
threads=0


WD=$(pwd)
while getopts ":a:t:p:e:s:o:n:r:c:1:2:m:d:g:h" opt; do
  case $opt in
    a) align="$OPTARG"
    ;;
    t) tree="$OPTARG"
    ;;
    p) read_one="$OPTARG"
    ;;
    e) PE=2;
       read_two="$OPTARG"
    ;;
    s) read_stub="$OPTARG"
    ;;
    o) outdir="$OPTARG"
	  ;;
	  n) nam="$OPTARG"
    ;;
    r) read_align="$OPTARG"
    ;;
    c) threads="$OPTARG"
    ;;
    1) r1_tail="$OPTARG"
    ;;
    2) r2_tail="$OPTARG"
    ;;
    m) align_type="$OPTARG"
    ;;
    d) master_dir="$OPTARG"
    ;;
    g) gap_align="$OPTARG"
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

if [ $threads -eq 0 ]; then
     threads=2
     printf "Threads not set, defaulting to 2" >&2
else
     echo "num threads is"
     echo $threads
fi

#######################################################################
assert ()                 #  If condition false,
{                         #+ exit from script
                          #+ with appropriate error message.
  E_PARAM_ERR=98
  E_ASSERT_FAILED=99


  if [ -z "$2" ]          #  Not enough parameters passed
  then                    #+ to assert() function.
    return $E_PARAM_ERR   #  No damage done.
  fi

  lineno=$2

  if [ ! $1 ]
  then
    echo "Assertion failed:  \"$1\""
    echo "File \"$0\", line $lineno"    # Give name of file and line number.
    exit $E_ASSERT_FAILED
  # else
  #   return
  #   and continue executing the script.
  fi
}
#######################################################################

printf "master dir == %s\n" "$master_dir"
printf "tail_1 is %s\n" "$r1_tail"
printf "tail_2 is %s\n" "$r2_tail"
printf "Argument out is %s\n" "$outdir"
printf "Argument name is %s\n" "$nam"
printf "Argument threads is %s\n" "$threads"
mkdir -p $outdir
#Check that tipnames in alignemnet are the same as tipnames in tree

echo 'Performing full mapping of reads to all sequences in alignment'
workd=$(pwd)
cd $workd

# TEST OF GAP REMOVAL FINISHED


### TODO PLAY WITH BOWTIE2 --very-fast command to chekc speed up time

#pretend the alignemnt is a set of chromosomes

#this is a hack that is in both scripts!! need to be passed between

#echo "PAIRED ENDS"
base=$(basename $read_one $r1_tail)
echo "basename is $base"
refnam=$(head -n 1 $align)

echo "refname is $refnam"

hisat_idx=$(echo "$align" | cut -f 1 -d '.')

echo "Time for bowtie2-build:"

#TOTDO THINK HARD ABOUT IMPLAICTIONS OF LOCAL VS GLOBAL AIGN!!!
echo "time for bowtie2 mapping:"
printf "\nPE = $PE\n" 
if [ "$PE" -eq 2 ]; then
	printf "\nRunning in Paired-End Mode\n"
	#bowtie2 -p $threads --very-fast -x $outdir/best_ref -1 $read_one -2 $read_two -S $outdir/best_map.sam --no-unal --local
	hisat2 -p $threads --very-fast -x $hisat_idx -1 $read_one -2 $read_two -S $outdir/best_map.sam --no-unal
	printf "\nRun with command: hisat2 -p $threads --very-fast -x $hisat_idx -1 $read_one -2 $read_two -S $outdir/best_map.sam --no-unal\n"

elif [ "$PE" -eq 0 ]; then
	printf "\nRunning in Single-End Mode\n"
	hisat2 -p $threads --very-fast -x $hisat_idx -U $read_one -S $outdir/best_map.sam --no-unal
	printf "\nRun with command: hisat2 -p $threads --very-fast -x $hisat_idx -U $read_one -S $outdir/best_map.sam --no-unal\n"
fi

# samtools faidx $align
# echo '>samtools faidx passed'

samtools view -bS $outdir/best_map.sam > $outdir/best_map.bam
echo '>samtools view passed'
samtools sort $outdir/best_map.bam -o $outdir/best_sorted.bam
echo '>samtools sort passed'
samtools index $outdir/best_sorted.bam
echo '>samtools index passed'

echo '>Time for mpileup step:'
#bcftools mpileup -f $outdir/best_ref.fas $outdir/best_sorted.bam -o $outdir/best_sorted.vcf
bcftools mpileup -f $align $outdir/best_sorted.bam -o $outdir/best_sorted.vcf

echo '>Dropping duplicate positions'
$PHYCORDER/vcf_duplicate_dropper.py --vcf_file $outdir/best_sorted.vcf --out_file $outdir/dupes_removed_best_cns.vcf

bcftools call -c $outdir/dupes_removed_best_cns.vcf -o $outdir/best_cns.vcf

vcfutils.pl vcf2fq $outdir/best_cns.vcf >  $outdir/cns.fq

echo '>samtools mpileup passed'
seqtk seq -a $outdir/cns.fq > $outdir/cns.fa
echo '>seqtk passed'

# make copy of cns.fa for later reference of best reference for mapping
cp $outdir/cns.fa $outdir/cns_ref.fa

# get name of reference in a variable to use for removal
refname_cns=$(head -1 $outdir/cns.fa)

sed -i -e "s/$refname_cns/>${base}/g" $outdir/cns.fa
echo '>sed producing cns.fa passed'

sed -i -e "s/>/>QUERY_${base}_ref_/g" $outdir/cns_ref.fa
echo '>sed producing cns_ref.fa passed'

# use vcffixer.py to find the N's in the cns.fa and replace with the real nucleotides
# from the reads. This seems to be a strange bug from all versions of mpileup
$PHYCORDER/vcffixer.py --vcf_file $outdir/best_sorted.vcf --align_file $outdir/cns.fa --out_file $outdir/cns_fixed.fa

refnam=$(head -n 1 $master_dir/ref_nogap.fas)

#$PHYCORDER/vcf_to_fasta.py --vcf_file $outdir/best_sorted.vcf --seq_name $base --out_file $outdir/cns_fixed.fa

printf ">beginning aligned consensus processing"
#python

$PHYCORDER/align_consensus.py --gapped-ref $gap_align --consensus $outdir/cns_fixed.fa --outfile $outdir/${base}_align.fas



# if [ $wre_map -eq 1 ]
#      then
#          #generate consensus mapped to second best locus to assuage some ref dependence
#          wrefnam=$(sort -rnk3 $outdir/mapping_info | head -2 | tail -1| cut -f1)
#          grep -Pzo '(?s)>'$wrefnam'.*?>' $outdir/ref_nogap.fas |head -n-1 > $outdir/worse_ref_uneven.fas
#          python fastafixer.py $outdir/worse_ref_uneven.fas $outdir/worse_ref.fas
#          echo 'The second best reference found in your alignment was '$wrefnam
#          echo 'mapping reads to '$wrefnam

#          bowtie2-build $outdir/worse_ref.fas $outdir/worse_ref >> $outdir/bowtiebuild.log

#          if [ $PE -eq 1 ];
#              then
#                  bowtie2 -x $outdir/worse_ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/worse_map.sam --no-unal --local
#              else
#                  bowtie2 -x $outdir/worse_ref  -U ${read_stub}.fastq -S $outdir/worse_map.sam --no-unal --local
#          fi
#          samtools faidx $outdir/worse_ref.fas
#          samtools view -bS $outdir/worse_map.sam > $outdir/worse_map.bam
#          samtools sort $outdir/worse_map.bam -o $outdir/worse_sorted.bam
#          samtools index $outdir/worse_sorted.bam
#          samtools mpileup -uf $outdir/worse_ref.fas $outdir/worse_sorted.bam| bcftools call -c | vcfutils.pl vcf2fq >  $outdir/worse_cns.fq
#          seqtk seq -a $outdir/worse_cns.fq > $outdir/worse_cns.fa
#          sed -i -e "s/>/>${nam}${read_one}_/g" $outdir/cns.fa
#          if [ $(diff $outdir/cns.fa $outdir/worse_cns.fa | wc -l | cut -f3) -gt 4 ]
#               then
#                   echo 'Alternate references result in different sequences. Placing both, but investigating differences recommended!'
#                   cat $outdir/cns.fa $outdir/worse_cns.fa > $outdir/mappings.fa
#                   cd $outdir
#                   grep -Pzo '(?s)>'$refnam'.*?>' $align |head -n-1 > $outdir/worse_ref_gaps.fas
#                   align_consensus.py --gapped-ref $outdir/worse_ref_gaps.fas --consensus $outdir/worse_cns.fa --outfile $outdir/worse_aln.fa
#                   cat $outdir/tmp.aln $outdir/worse_align.fa > $outdir/tmp2.aln
#                   raxmlHPC -m GTRGAMMA -s $outdir/tmp.aln -t $tree -p 12345 -n consensusFULL

#                  cd $WD
# #         else
# #             echo 'Using worse reference resulted in identical sequences - only aligning and placing one.'
# #             cd $outdir
# #                 papara -t ${WD}/${tree} -s ${aln_stub}.phy -q cns.fa -n re_consensus
# #                 raxmlHPC -m GTRCAT -f v -s papara_alignment.re_consensus -t ${WD}/$tree -n ${nam}_consensusPC
# #             cd $WD
# # fi
