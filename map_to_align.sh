#!/bin/bash
#inputs: an alignment, a tree, and whole genome reads
#determine how to try differnet alignments
#deal with possibility of multiple samples?
#HOw to deal with single mapping to multiple alignements??!?
#./map_to_align.sh -a example.aln -t tree.tre -p /home/ejmctavish/projects/Exelixis/SISRS/full_aln/datafiles/SRR610374 -o fulltest -n fulltest #
#DEFAULT ARGS
papara=/home/ejmctavish/projects/Exelixis/papara_nt-2.4/papara
EPAOME=/home/ejmctavish/projects/Exelixis/EPAome

PE=0
outdir=EPAome_run
nam=QUERY
read_align=0
re_map=1
map=1
read_name_prefix=SRR

wre_map=0

WD=$(pwd)
while getopts ":a:t:p:s:o:n:r:m:b:w:" opt; do
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
    r) read_align="$OPTARG"
    ;;
    m) map="$OPTARG"
    ;;
    b) re_map="$OPTARG"
    ;;
    w) wre_map="$OPTARG"
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
printf "Argument read_align is %s\n" "$read_align"
printf "Argument map is %s\n" "$map"
printf "Argument re_mapis %s\n" "$re_map"

mkdir -p $outdir
if [ $map -eq 1 ];
    then
    echo 'Performing full mapping of reads to all sequences in alignment'
    sed 's/-//g' <$align >$outdir/ref_nogap.fas
    bowtie2-build $outdir/ref_nogap.fas $outdir/ref > $outdir/bowtiebuild.log
    if [ $PE -eq 1 ];
    	then 
    	    echo "PAIRED ENDS"
    	    bowtie2 -x $outdir/ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/full_alignment.sam --no-unal
        else 
        	bowtie2 -x $outdir/ref -U ${read_stub}.fastq -S $outdir/full_alignment.sam --no-unal
    fi

    samtools view -bS $outdir/full_alignment.sam > $outdir/full_alignment.bam
    samtools sort $outdir/full_alignment.bam $outdir/full_sorted
    samtools index $outdir/full_sorted.bam 
    samtools idxstats $outdir/full_sorted.bam > $outdir/mapping_info
    if [ $(sort -rnk3 $outdir/mapping_info | head -1 | cut -f3) -lt 10 ]; then
        echo 'LESS THAN TEN READS MAPPED TO ANY LOCUS. Try a different input alignment?'
        exit
    fi
    #assert at least some reads mapped!! 
fi

if [ $read_align -eq 1 ]
    then 
        echo 'Attempting to align and place all mapped reads'
        grep $read_name_prefix $outdir/full_alignment.sam | cut -f1 | uniq > $outdir/matches #ToDo this relies on read names starting with SRR. Need better approach
        if [ $(wc -l $outdir/matches | cut -f1 -d' ') -lt  10 ]; 
            then
               echo 'error in matched read grepping'
               exit
        fi
        if [ $PE -eq 1 ];
            then 
                seqtk subseq ${read_stub}_1.fastq $outdir/matches > $outdir/matches.fq
                seqtk subseq ${read_stub}_2.fastq $outdir/matches > $outdir/matches2.fq
                rnam=$(head -n 1 $outdir/matches | awk '{print substr($0,0,4)}')
                sed -i s/@$rnam/@R2$rnam/ $outdir/matches2.fq #these two lines are so that paired end reads have different names. Adds R2 to beginning of read string
                cat $outdir/matches2.fq >> $outdir/matches.fq
            else 
                seqtk subseq ${read_stub}.fastq $outdir/matches > $outdir/matches.fq
        fi
        fastq_to_fasta -i $outdir/matches.fq -o $outdir/matches.fa
        fastx_collapser < $outdir/matches.fa > $outdir/matches_unique.fa
        aln_stub=$(echo $align | cut -f1 -d.)
        python $EPAOME/fasta_to_phylip.py $align $outdir/$aln_stub.phy
        cd $outdir
            $papara -t ${WD}/${tree} -s ${aln_stub}.phy -q matches_unique.fa -n reads #NOTE, Quotes in trees cause issues. From dendropy or elsewhere?!
            raxmlHPC -m GTRCAT -f v -s papara_alignment.reads -t ${WD}/${tree} -n ${nam}_reads_EPA
        cd $WD
fi    

if [ $re_map -eq 1 ]
    then 
        echo 'Refining mapping and calling consensus sequence'
        refnam=$(sort -rnk3 $outdir/mapping_info | head -1 | cut -f1)
        grep -Pzo '(?s)>'$refnam'.*?>' $outdir/ref_nogap.fas |head -n-1 > $outdir/best_ref_uneven.fas
        python fastafixer.py $outdir/best_ref_uneven.fas $outdir/best_ref.fas
        echo 'The best reference found in your alignment was '$refnam
        echo 'mapping reads to '$refnam

        bowtie2-build $outdir/best_ref.fas $outdir/best_ref >> $outdir/bowtiebuild.log

        #TOTDO THINK HARD ABOUT IMPLAICTIONS OF LOCAL VS GLOBAL AIGN!!!
        if [ $PE -eq 1 ]
        	then 
        	    echo "PAIRED ENDS"
        	    bowtie2 -x $outdir/best_ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/best_map.sam --no-unal --local
            else 
            	bowtie2 -x $outdir/best_ref  -U ${read_stub}.fastq -S $outdir/best_map.sam --no-unal --local
        fi

        samtools faidx $outdir/best_ref.fas
        samtools view -bS $outdir/best_map.sam > $outdir/best_map.bam
        samtools sort $outdir/best_map.bam $outdir/best_sorted
        samtools index $outdir/best_sorted.bam 
        samtools mpileup -uf $outdir/best_ref.fas $outdir/best_sorted.bam| bcftools call -c | vcfutils.pl vcf2fq >  $outdir/cns.fq 
        python ~/projects/Exelixis/EPAome/samtoolsfq_to_fa.py $outdir/cns.fq $outdir/cns.fa $nam
        if [ $wre_map -eq 1 ]
            then
                #generate consensus mapped to second best locus to assuage some ref dependence
                wrefnam=$(sort -rnk3 $outdir/mapping_info | head -2 | tail -1| cut -f1)
                grep -Pzo '(?s)>'$refnam'.*?>' $outdir/ref_nogap.fas |head -n-1 > $outdir/worse_ref_uneven.fas
                python fastafixer.py $outdir/worse_ref_uneven.fas $outdir/worse_ref.fas
                echo 'The second best reference found in your alignment was '$refnam
                echo 'mapping reads to '$refnam

                bowtie2-build $outdir/worse_ref.fas $outdir/worse_ref >> $outdir/bowtiebuild.log

                if [ $PE -eq 1 ];
                    then 
                        echo "PAIRED ENDS"
                        bowtie2 -x $outdir/worse_ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/worse_map.sam --no-unal --local
                    else 
                        bowtie2 -x $outdir/worse_ref  -U ${read_stub}.fastq -S $outdir/worse_map.sam --no-unal --local
                fi
                samtools faidx $outdir/worse_ref.fas 
                samtools view -bS $outdir/worse_map.sam > $outdir/worse_map.bam
                samtools sort $outdir/worse_map.bam $outdir/worse_sorted
                samtools index $outdir/worse_sorted.bam 
                samtools mpileup -uf $outdir/worse_ref.fas $outdir/worse_sorted.bam| bcftools call -c | vcfutils.pl vcf2fq >  $outdir/worse_cns.fq 

                python ~/projects/Exelixis/EPAome/samtoolsfq_to_fa.py $outdir/worse_cns.fq $outdir/worse_cns.fa worse_query
                cat $outdir/cns.fa $outdir/worse_cns.fa > $outdir/mappings.fa
                fastx_collapser < $outdir/mappings.fa > $outdir/mappings_unique.fa
                aln_stub=$(echo $align | cut -f1 -d.)
                python $EPAOME/fasta_to_phylip.py $align $outdir/$aln_stub.phy
                cd $outdir
                  $papara -t ${WD}/${tree} -s ${aln_stub}.phy -q mappings_unique.fa -n consensus 
                  #run RAXML EPA on the alignments
                  raxmlHPC -m GTRCAT -f v -s papara_alignment.consensus -t ${WD}/$tree -n ${nam}_consensusEPA
                cd $WD
            else
                $papara -t ${WD}/${tree} -s ${aln_stub}.phy -q cns.fa -n consensus 
                raxmlHPC -m GTRCAT -f v -s papara_alignment.consensus -t ${WD}/$tree -n ${nam}_consensusEPA
            fi
            #run full raxml? tooo sloooo
       # raxmlHPC -m GTRGAMMA -s $outdir/contig_alignment.fas -t $tree -p 12345 -n consensusFULL

fi




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