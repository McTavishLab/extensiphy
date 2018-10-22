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

# BOWTIE2PATH="/home/ejmctavish/Applications/bowtie2-2.3.4.2-linux-x86_64"
printf "phycorder directory is %s\n" "$PHYCORDER"



#Check for dependencies
if [ $(which bcftools | wc -l) -lt 1 ]
    then
        printf "Requires bcftools" >&2
     #   exit 0
    else
        printf "Correct version of bfctools found.\n"
fi
if [  $(which samtools | wc -l) -lt 1 ] #TODO steup for greater than 1.2? this  is a sloppppy approach
    then
        printf "Requires samtools" >&2
      #  exit 0
    else
        printf "Correct version of samtools found.\n"
fi
if [ $(which seqtk | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "seqtk not found\n" >&2
    else
        printf "seqtk found\n"
fi
if [ $(which bowtie2 | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "bowtie2 not found. Install and/or add to path\n" >&2
    else
        printf "bowtie2 found\n"
fi
if [ $(which raxmlHPC-PTHREADS-SSE3 | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "raxmlHPC not found. Install and/or alias or add to path\n" >&2
    else
        printf "raxmlHPC found\n"
fi
if [ $(which fastx_collapser | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "fastx toolkit not found. Install and/or add to path\n" >&2
    else
        printf "fastx toolkit found\n"
fi
if [ $(which vcfutils.pl | wc -l) -lt 1 ] #TODO needs different install than bcftools?
    then
        printf "vcfutils.pl not found. Install and/or add to path\n" >&2
    else
        printf "vcfutils.pl found\n"
fi

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
while getopts ":a:t:p:e:s:o:n:r:c:1:2:h" opt; do
  case $opt in
    a) align="$OPTARG"
    ;;
    t) tree="$OPTARG"
    ;;
    p) PE=1;
       read_one="$OPTARG"
    ;;
    e)
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
# # if [ $PE -eq 1 ]; then
# #   if [ -f ${read_one} ]; then
# #      printf "Paired end reads \n"
# #      printf "read one is ${read_one}\n"
# #   else
# #     printf "read one ${read_one} not found. Exiting\n" >&2
# #     exit
# # fi
# #   if [ -f ${read_two} ]; then
# #      printf "Paired end reads \n"
# #      printf "read two is ${read_two}\n"
# #   else
# #     printf "read two ${read_two} not found. Exiting\n" >&2
# #     exit
# # fi
# fi
printf "tail_1 is %s\n" "$r1_tail"
printf "tail_2 is %s\n" "$r2_tail"
printf "Argument out is %s\n" "$outdir"
printf "Argument name is %s\n" "$nam"
printf "Argument threads is %s\n" "$threads"
mkdir -p $outdir

#Check that tipnames in alignemnet are the same as tipnames in tree

echo 'Performing full mapping of reads to all sequences in alignment'

#pull all the gaps from the aligned taxa bc mappers cannot cope.
sed 's/-//g' <$align >$outdir/ref_nogap.fas

#printf ">THIS IS THE FIRST TEST SECTION!!!!!!!!!!!!\n"
#cat $outdir/ref_nogap.fas


### TODO PLAY WITH BOWTIE2 --very-fast command to chekc speed up time

#pretend the alignemnt is a set of chromosomes

#this is a hack that is in both scripts!! need to be passed between

echo "PAIRED ENDS"
base=$(basename $i $r1_tail)
echo "basename is $base"
mkdir -p ${base}_outdir
# bowtie2 -p $threads --very-fast -x $outdir/ref -1 $i -2 ${i%$r1_tail}$r2_tail -S $outdir/${base}_outdir/full_alignment.sam --no-unal
# printf ">map PE 1 passed"
#
# samtools view -bS $outdir/${base}_outdir/full_alignment.sam > $outdir/${base}_outdir/full_alignment.bam
#
# samtools sort $outdir/${base}_outdir/full_alignment.bam -o $outdir/${base}_outdir/full_sorted.bam
#
# samtools index $outdir/${base}_outdir/full_sorted.bam
#
# samtools idxstats $outdir/${base}_outdir/full_sorted.bam > $outdir/${base}_outdir/mapping_info
#
# if [ $(sort -rnk3 $outdir/${base}_outdir/mapping_info | head -1 | cut -f3) -lt 10 ]; then
#     echo 'LESS THAN TEN READS MAPPED TO ANY TAXON. Try a different input alignment?' >&2
#     exit
# fi
#
# echo '>Refining mapping and calling consensus sequence'
# sort -rnk3 $outdir/${base}_outdir/mapping_info >  $outdir/${base}_outdir/mapping_info_sort
#
refnam=$(head -n 1 ref_nogap.fas)

echo "refname is $refnam"
#
#
grep -Pzo '(?s)'$refnam'.*?(>|\Z)' ref_nogap.fas |head -n-1 > ${base}_outdir/best_ref_uneven.fas
# grep -Pzo '(?s)>'$refnam'.*?(>|\Z)' $align |head -n-1 > $outdir/best_ref_gaps.fas

fold -w 80 $outdir?${base}_outdir/best_ref_uneven.fas > best_ref.fas
# #printf ">going to fastafixer"
# $PHYCORDER/fastafixer.py $outdir/${base}_outdir/best_ref_uneven.fas $outdir/${base}_outdir/best_ref.fas #starightens out line lengths
# echo '>The best reference found in your alignment was '$refnam
# echo '>mapping reads to '$refnam

# bowtie2-build --threads $threads $outdir/${base}_outdir/best_ref.fas $outdir/${base}_outdir/best_ref >> $outdir/${base}_outdir/bowtiebuild.log
bowtie2-build --threads $threads ${base}_outdir/best_ref.fas ${base}_outdir/best_ref >> ${base}_outdir/bowtiebuild.log

#TOTDO THINK HARD ABOUT IMPLAICTIONS OF LOCAL VS GLOBAL AIGN!!!
time bowtie2 -p $threads --very-fast -x ${base}_outdir/best_ref -1 $read_one -2 $read_two -S ${base}_outdir/best_map.sam --no-unal --local


samtools faidx ${base}_outdir/best_ref.fas
echo '>samtools faidx passed'

samtools view -bS ${base}_outdir/best_map.sam > ${base}_outdir/best_map.bam
echo '>samtools view passed'
samtools sort ${base}_outdir/best_map.bam -o ${base}_outdir/best_sorted.bam
echo '>samtools sort passed'
samtools index ${base}_outdir/best_sorted.bam
echo '>samtools index passed'
time samtools mpileup -uf ${base}_outdir/best_ref.fas ${base}_outdir/best_sorted.bam| bcftools call -c | vcfutils.pl vcf2fq >  ${base}_outdir/cns.fq
echo '>samtools mpileup passed'
seqtk seq -a ${base}_outdir/cns.fq > ${base}_outdir/cns.fa
echo '>seqtk passed'

#automatic naming names it to the ref wich is confusing
#sed -i -e "s/>/>${nam}${read_one}_/g" $outdir/cns.fa
sed -i -e "s/>/>QUERY_${base}_ref_/g" ${base}_outdir/cns.fa
echo '>sed producing cns.fa passed'

# TODO this might be a DANGEROUS way to handle this issue. think of alternate
# cat $align | head -1 | cut -f1 | cut -c 2- > $outdir/${base}_outdir/best_ref_gaps_name.fas

#refnam=$(cat $outdir/${base}_outdir/best_ref_gaps_name.fas | head -1)

 refnam=$(head -n 1 ref_nogap.fas)

#pull the aligned reference from the alignement
#grep -Pzo '(?s)>'$refnam'.*?>' $align |head -n-1 > $outdir/${base}_outdir/best_ref_gaps.fas
grep -Pzo '(?s)'$refnam'.*?(>|\Z)' $align |head -n-1 > ${base}_outdir/best_ref_gaps.fas

echo '>grep for refnam passed'

printf ">beginning aligned consensus processing"
#python

$PHYCORDER/align_consensus.py --gapped-ref ${base}_outdir/best_ref_gaps.fas --consensus ${base}_outdir/cns.fa --outfile "${base}_outdir"/"${base}_align.fas"


# cat ${align} $outdir/aligned_cns.fas >  $outdir/extended.aln

# cd $outdir
# run full raxml? tooo sloooo
# raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -s extended.aln -t $tree -p 12345 -n consensusFULL

# cd $WD

#todo strip all fq to fa



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
