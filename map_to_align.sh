
#!/bin/bash
#inputs: an alignment, a tree, and whole genome reads 

#do do set up paired end vs not flags
#determine how to try differnet alignments
#deal with possibility of multiple samples?
#HOw to deal with single mapping to multiple alignements??!?

#DEFAULT ARGS
PE=0
outdir=EPAome_run
nam=QUERY

while getopts ":a:t:p:s:o:n:" opt; do
  case $opt in
    a) align="$OPTARG"
    ;;
    t) tree="$OPTARG"
    ;;
    p) PE=1;
       read_stub="$OPTARG"
    ;;
    s) read_stub="$OPTARG"
    ;;
    o) outdir="$OPTARG"
	;;
	n) nam="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "Argument alightis %s\n" "$align"
printf "Argument tree is %s\n" "$tree"
printf "Argument PE is %s\n" "$PE"
printf "Argument stub is %s\n" "$read_stub"
printf "Argument out is %s\n" "$outdir"
printf "Argument name is %s\n" "$nam"

#read_stub=~/projects/Exelixis/SISRS/full_aln/datafiles/SRR610374
#fasta = ~/projects/Exelixis/SISRS/fasta/SRR610375.fasta


mkdir -p $outdir
sed 's/-//g' <$align >$outdir/ref_nogap.fas

bowtie2-build $outdir/ref_nogap.fas $outdir/ref
if [ $PE -eq 1 ];
	then 
	    echo "PAIRED ENDS"
	    bowtie2 -x $outdir/ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/full_alignment.sam --no-unal;
    else 
    	bowtie2 -x $outdir/ref -U ${read_stub}.fastq -S $outdir/full_alignment.sam --no-unal;
fi

samtools view -bS $outdir/full_alignment.sam > $outdir/full_alignment.bam
samtools sort $outdir/full_alignment.bam $outdir/full_sorted
samtools index $outdir/full_sorted.bam 
samtools idxstats $outdir/full_sorted.bam > $outdir/mapping_info
nam=`sort -rnk3 $outdir/mapping_info | head -1 | cut -f1`
if ((`sort -rnk3 $outdir/mapping_info | head -1 | cut -f3`<10)); then
    echo 'LESS THAN TEN READS MAPPED TO ANY LOCUS. Try a different input alignment?'
    exit
fi
#assert at least some reads mapped!! 

grep -Pzo '(?s)>'$nam'.*?>' $outdir/ref_nogap.fas |head -n-1 > $outdir/best_ref.fas

echo 'The best reference found in your alignment was '$nam
echo 'mapping reads to '$nam

bowtie2-build $outdir/best_ref.fas $outdir/best_ref

#TOTDO THINK HARD ABOUT IMPLAICTIONS OF LOCAL VS GLOBAL AIGN!!!

if [ $PE -eq 1 ];
	then 
	    echo "PAIRED ENDS"
	    bowtie2 -x $outdir/best_ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/best_map.sam --no-unal
    else 
    	bowtie2 -x $outdir/best_ref  -U ${read_stub}.fastq -S $outdir/best_map.sam --no-unal;
fi

#

samtools view -bS $outdir/best_map.sam > $outdir/best_map.bam
samtools sort $outdir/best_map.bam $outdir/best_sorted
samtools index $outdir/best_sorted.bam 
samtools mpileup -uf $outdir/best_ref.fas $outdir/best_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq > $outdir/cns.fq  

#uhh works better twice?
samtools mpileup -uf $outdir/best_ref.fas $outdir/best_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq > $outdir/cns.fq  
#need better fq to FA solution!!!!

python ~/projects/Exelixis/EPAome/samtoolsfq_to_fa.py $outdir/cns.fq $outdir/cns.fa $nam

~/projects/Exelixis/pagan-msa/src/pagan --ref-seqfile $align -t $tree --queryfile $outdir/cns.fa  --outfile $outdir/contig_alignment

#run RAXML EPA on the alignments

raxmlHPC -m GTRCAT -f v -s $outdir/contig_alignment.fas -t $tree -n $nam_consensusEPA

#run full raxml? tooo sloooo
raxmlHPC -m GTRGAMMA -s $outdir/contig_alignment.fas -t $tree -p 12345 -n consensusFULL



#----------playing with denovo locus assmbly, doesn't work RN-------------------------------------------------------
#grep SRR $outdir/full_alignment.sam | cut -f1 | uniq > $outdir/matches

#python matchgrabber.py analyses/matches $fasta
#seqtk subseq ${read_stub}_R1.fastq   $outdir/matches > $outdir/mapped.fq
#seqtk subseq ${read_stub}_R2.fastq   $outdir/matches >> $outdir/mapped.fq

#assemble all the reads that mapped anywhere in the alignment into contigs
#minia $outdir/mapped.fq 7 3 1000 $outdir/assemb


#also align all mapped reads to the alignment
#~/projects/Exelixis/pagan-msa/src/pagan --ref-seqfile $align -t $tree --queryfile $outdir/mapped.fq --outfile read_alignment

#try aligning contigs too
#~/projects/Exelixis/pagan-msa/src/pagan --ref-seqfile $align -t $tree --queryfile $outdir/assemb.contigs.fa --outfile contig_alignment


#run RAXML EPA on the alignments
#raxmlHPC ­-f v ­-s read_alignment ­-t $tree ­-m GTRCAT ­-n contigs



#run RAXML EPA on the alignments
#raxmlHPC ­-f v ­-s contig_alignment ­-t $tree ­-m GTRCAT ­-n reads