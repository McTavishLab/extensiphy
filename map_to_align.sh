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

WD=$(pwd)
while getopts ":a:t:p:e:s:o:n:r:m:b:w:c:h" opt; do
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

# if [ -z "$align" ] || [ -z "$tree" ]; then
#    "alignment (-a), tree (-t), and reads (-p or -s required)"
#    exit
# fi
#
# #Ttest if files actually exist
# #Check to make sure mapping has occured if re-mapping
#
# if [ -f "$align" ]; then
#     printf "Alignment is %s\n" "$align"
#   else
#     printf "Alignment $align not found. Exiting\n" >&2
#     exit
# fi
# if [ -f "$tree" ]; then
#     printf "Tree is %s\n" "$tree"
#   else
#     printf "Tree $tree not found. Exiting\n" >&2
#     exit
# fi
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
printf "Argument out is %s\n" "$outdir"
printf "Argument name is %s\n" "$nam"
printf "Argument map is %s\n" "$map"
printf "Argument re_mapis %s\n" "$re_map"

mkdir -p $outdir

#Check that tipnames in alignemnet are the same as tipnames in tree

echo 'Performing full mapping of reads to all sequences in alignment'

#pull all the gaps from the aligned taxa bc mappers cannot cope.
sed 's/-//g' <$align >$outdir/ref_nogap.fas

printf ">THIS IS THE FIRST TEST SECTION!!!!!!!!!!!!\n"
cat $outdir/ref_nogap.fas


### TODO PLAY WITH BOWTIE2 --very-fast command to chekc speed up time

#pretend the alignemnt is a set of chromosomes
bowtie2-build --threads $threads $outdir/ref_nogap.fas $outdir/ref > $outdir/bowtiebuild.log
printf ">build 1 passed"

if [ $PE -eq 1 ];
	then
	    echo "PAIRED ENDS"
	    bowtie2 -p $threads --very-fast -x $outdir/ref -1 ${read_one} -2 ${read_two} -S $outdir/full_alignment.sam --no-unal
      printf ">map PE 1 passed"
    else
      bowtie2 -p $threads --very-fast -x $outdir/ref -U ${read_one}-S $outdir/full_alignment.sam --no-unal
      printf ">map single 1 passed"
fi

samtools view -bS $outdir/full_alignment.sam > $outdir/full_alignment.bam

samtools sort $outdir/full_alignment.bam -o $outdir/full_sorted.bam

samtools index $outdir/full_sorted.bam

samtools idxstats $outdir/full_sorted.bam > $outdir/mapping_info

printf ">samtools passed"
if [ $(sort -rnk3 $outdir/mapping_info | head -1 | cut -f3) -lt 10 ]; then
    echo 'LESS THAN TEN READS MAPPED TO ANY TAXON. Try a different input alignment?' >&2
    exit
fi
    #TODO this is VERY DANGEROUS


echo '>Refining mapping and calling consensus sequence'
refnam=$(sort -rnk3 $outdir/mapping_info | head -1 | cut -f1)

grep -Pzo '(?s)>'$refnam'.*?>' $outdir/ref_nogap.fas |head -n-1 > $outdir/best_ref_uneven.fas
printf ">THIS IS A TEST SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
cat $outdir/best_ref_uneven.fas
printf ">TEST OVER\n"

printf ">going to fastafixer"
$PHYCORDER/fastafixer.py $outdir/best_ref_uneven.fas $outdir/best_ref.fas #starightens out line lengths
echo '>The best reference found in your alignment was '$refnam
echo '>mapping reads to '$refnam

printf ">TEST SEGMENT BEGINNING\n"
printf ">best_ref_uneven.fas"
cat $outdir/best_ref_uneven.fas
printf ">best_ref.fas"
cat $outdir/best_ref.fas
printf ">TEST OVER\n"

bowtie2-build --threads $threads $outdir/best_ref.fas $outdir/best_ref >> $outdir/bowtiebuild.log

#TOTDO THINK HARD ABOUT IMPLAICTIONS OF LOCAL VS GLOBAL AIGN!!!
if [ $PE -eq 1 ]
	then
	    bowtie2 -p $threads --very-fast -x $outdir/best_ref -1 ${read_one} -2 ${read_two} -S $outdir/best_map.sam --no-unal --local
    else
    	bowtie2 -p $threads --very-fast -x $outdir/best_ref  -U ${read_one} -S $outdir/best_map.sam --no-unal --local
fi

samtools faidx $outdir/best_ref.fas

printf ">TEST SEGMENT BEGINNING\n"
printf ">best_ref.fas"
cat $outdir/best_ref.fas
printf ">TEST FINISHED\n"

samtools view -bS $outdir/best_map.sam > $outdir/best_map.bam
samtools sort $outdir/best_map.bam -o $outdir/best_sorted.bam
samtools index $outdir/best_sorted.bam
samtools mpileup -uf $outdir/best_ref.fas $outdir/best_sorted.bam| bcftools call -c | vcfutils.pl vcf2fq >  $outdir/cns.fq
seqtk seq -a $outdir/cns.fq > $outdir/cns.fa

printf ">TEST SEGMENT BEGINNING\n"
printf ">cns.fa"
cat $outdir/cns.fa
printf ">TEST FINISHED"

#automatic naming names it to the ref wich is confusing
#sed -i -e "s/>/>${nam}${read_one}_/g" $outdir/cns.fa
sed -i -e "s/>/>QUERY_/g" $outdir/cns.fa

printf ">TEST SEGMENT cns.fa"
cat $outdir/cns.fa
printf ">TEST FINISHED\n"

#pull the aligned reference from the alignement
grep -Pzo '(?s)>'$refnam'.*?>' $align |head -n-1 > $outdir/best_ref_gaps.fas
printf ">TEST SEGMENT best_ref_gaps.fas"
cat $outdir/best_ref_gaps.fas
printf ">TEST FINISHED\n"

printf ">beginning aligned consensus processing"
#python

y=${read_one%.fastq}
cns_align=${y##*/}
$PHYCORDER/align_consensus.py --gapped-ref $outdir/best_ref_gaps.fas --consensus $outdir/cns.fa --outfile "$outdir"/"$cns_align.fas"

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
