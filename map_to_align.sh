
#!/bin/bash
#inputs: an alignment, a tree, and whole genome reads 



align=$1
tree=tree.tre
#read_stub=~/projects/Exelixis/SISRS/full_aln/datafiles/SRR610374
#fasta = ~/projects/Exelixis/SISRS/fasta/SRR610375.fasta
echo "alignment is "$align
read_stub=SRR1019272
outdir=analyses_new
mkdir -p $outdir
sed 's/-//g' <$align >$outdir/ref_nogap.fas
bowtie2-build $outdir/ref_nogap.fas $outdir/ref
#bowtie2 -x $outdir/ref -1 ${read_stub}_R1.fastq -2 ${read_stub}_R2.fastq -S $outdir/full_alignment.sam --no-unal
bowtie2 -x $outdir/ref -U ${read_stub}.fastq -S $outdir/full_alignment.sam --no-unal

#fastq_to_fasta -in ${read_stub}_R1.fastq -out analyses ${read_stub}_R1.fasta
#fastq_to_fasta -in ${read_stub}_R2.fastq -out analyses ${read_stub}_R2.fasta

#cat analyses ${read_stub}_R1.fasta analyses ${read_stub}_R2.fasta > analyses ${read_stub}_R1.fasta
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
#bowtie2 -x $outdir/best_ref -1 ${read_stub}_R1.fastq -2 ${read_stub}_R2.fastq -S $outdir/best_map.sam --no-unal
bowtie2 -x $outdir/best_ref  -U ${read_stub}.fastq -S $outdir/full_alignment.sam --no-unal

samtools view -bS $outdir/best_map.sam > $outdir/best_map.bam
samtools sort $outdir/best_map.bam $outdir/best_sorted
samtools index $outdir/best_sorted.bam 

##RISKKKK!!! need to make sure goes in N's in case of no alignment/poor quality
samtools mpileup -uf $outdir/best_ref.fas $outdir/best_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq > $outdir/cns.fq  

#need better fq to FA solution!!

~/projects/Exelixis/pagan-msa/src/pagan --ref-seqfile $align -t $tree --queryfile $outdir/cns.fa  --outfile $outdir/contig_alignment



#run RAXML EPA on the alignments
raxmlHPC -m GTRCAT -f v -s $outdir/contig_alignment.fas -t tree.tre -n $runname_consensusEPA



#run full raxml
raxmlHPC ­-s $outdir/contig_alignment -m GTRGAMMA -t $tree -p 12345 -n consensusFULL
#raxmlHPC -m GTRCAT ­­-f v ­-s $outdir/contig_alignment ­-t $tree ­--n contigs


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