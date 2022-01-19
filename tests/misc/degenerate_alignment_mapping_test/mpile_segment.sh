#!/bin/bash


while getopts ":a:t:p:e:s:o:n:r:c:1:2:h" opt; do
  case $opt in
    a) read_one="$OPTARG"
    ;;
    e) read_two="$OPTARG"
    ;;
    o) outdir="$OPTARG"
    ;;
    c) threads="$OPTARG"
    ;;
    r) ref="$OPTARG"
    ;;
    h) echo  "paired end reads in read set one (-a), read set two (-e), alignment that has had a bowtie-build reference created for it (-r), and outdirectory (-o) and the number of threads you'd like to use (-c)"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

time bowtie2-build $outdir/$ref $outdir/$ref

echo "time for bowtie2 mapping:"
time bowtie2 -p $threads --very-fast -x $outdir/$ref -1 $read_one -2 $read_two -S $outdir/best_map.sam --no-unal --local


samtools faidx $outdir/$ref
echo '>samtools faidx passed'

samtools view -bS $outdir/best_map.sam > $outdir/best_map.bam
echo '>samtools view passed'
samtools sort $outdir/best_map.bam -o $outdir/best_sorted.bam
echo '>samtools sort passed'
samtools index $outdir/best_sorted.bam
echo '>samtools index passed'

echo '>Time for mpileup step:'
time bcftools mpileup -f $outdir/$ref $outdir/best_sorted.bam -o $outdir/best_sorted.vcf

time bcftools call -c $outdir/best_sorted.vcf -o $outdir/best_cns.vcf

time vcfutils.pl vcf2fq $outdir/best_cns.vcf >  $outdir/cns.fq

echo '>samtools mpileup passed'
seqtk seq -a $outdir/cns.fq > $outdir/cns.fa
echo '>seqtk passed'
