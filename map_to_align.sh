#!/bin/bash
#inputs: an alignment, a tree, and whole genome reads
#determine how to try differnet alignments
#deal with possibility of multiple samples?
#HOw to deal with single mapping to multiple alignements??!?
#./map_to_align.sh -a example.aln -t tree.tre -p /home/ejmctavish/projects/Exelixis/SISRS/full_aln/datafiles/SRR610374 -o fulltest -n fulltest #
#DEFAULT ARGS
PHYCORDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

printf "phycorder directory is %s\n" "$PHYCORDER"



#Check for dependencies
if [ $(bcftools -v  | grep 1.2 | wc -l) -lt 1 ]
    then
        printf "Requires bcftools v. 1.2. Exiting\n" >&2 
     #   exit 0
    else
        printf "Correct version of bfctools found.\n"
fi
if [ $(samtools 2>&1 >/dev/null | grep 1.2 | wc -l ) -lt 1 ] #TODO steup for greater than 1.2? this  is a sloppppy approach
    then
        printf "Requires samtools v. 1.2. Exiting\n" >&2 
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
if [ $(which raxmlHPC | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "raxmlHPC not found. Install and/or add to path\n" >&2 
    else
        printf "raxmlHPC found\n"
fi
if [ $(which fastx_collapser | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "fastx toolkit not found. Install and/or add to path\n" >&2 
    else
        printf "fastx toolkit found\n"
fi
if [ $(which papara | wc -l) -lt 1 ] #TODO steup for greater than 1.2?
    then
        printf "papara not found. Install and/or add to path\n" >&2 
    else
        printf "papara  found\n"
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
while getopts ":a:t:p:s:o:n:r:m:b:w:h" opt; do
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
    h) echo  "alignment in fasta format (-a), tree in Newick format (-t), and reads in fastq (-p paired_end_base_filename or -s single_end_base_filename required)"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z "$align" ] || [ -z "$tree" ] || [ -z "$read_stub" ]; then
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
if [ -f "$tree" ]; then  
    printf "Tree is %s\n" "$tree"
  else
    printf "Tree $tree not found. Exiting\n" >&2
    exit
fi
if [ $PE -eq 1 ]; then
  if [ -f "${read_stub}_1.fastq" ]; then      
     printf "Paired end reads \n"
     printf "read one is ${read_stub}_1.fastq\n"
  else
    printf "read one ${read_stub}_1.fastq not found. Exiting\n" >&2 
    exit
  fi
else
  if [ -f "${read_stub}.fastq" ]; then      
    printf "reads are ${read_stub}.fastq\n"
  else
    printf "read one ${read_stub}.fastq not found. Exiting\n" >&2 
    exit
  fi
fi
printf "Argument out is %s\n" "$outdir"
printf "Argument name is %s\n" "$nam"
printf "Argument read_align is %s\n" "$read_align"
printf "Argument map is %s\n" "$map"
printf "Argument re_mapis %s\n" "$re_map"


mkdir -p $outdir
aln_stub=$(echo $align | cut -f1 -d.)
#python $PHYCORDER/fasta_to_phylip.py --input-fasta $align --output-phy $outdir/$aln_stub.phy
python fasta_to_phylip.py --input-fasta  $align --output-phy $outdir/$aln_stub.phy
#Check that tipnames in alignemnet are the same as tipnames in tree

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
    samtools sort $outdir/full_alignment.bam -o $outdir/full_sorted.bam
    samtools index $outdir/full_sorted.bam 
    samtools idxstats $outdir/full_sorted.bam > $outdir/mapping_info
    if [ $(sort -rnk3 $outdir/mapping_info | head -1 | cut -f3) -lt 10 ]; then
        echo 'LESS THAN TEN READS MAPPED TO ANY TAXON. Try a different input alignment?' >&2
        exit
    fi 
    #TODO this is VERY DANGEROUS
fi

#EJM: I think it takes the reads that mapped anywhere and maps them to the best ref
if [ $read_align -eq 1 ]
    then 
        echo 'Attempting to align and place all mapped reads'
        grep $read_name_prefix $outdir/full_alignment.sam | cut -f1 | uniq > $outdir/matches #ToDo this relies on read names starting with SRR. Need better approach
        if [ $(wc -l $outdir/matches | cut -f1 -d' ') -lt  10 ]; 
            then
               echo 'error in matched read grepping' >&2
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
        cd $outdir
            papara -t ${WD}/${tree} -s ${aln_stub}.phy -q matches_unique.fa -n reads #NOTE, Quotes in trees cause issues. From dendropy or elsewhere?!
            raxmlHPC -m GTRCAT -f v -s papara_alignment.reads -t ${WD}/${tree} -n ${nam}_reads_EPA
        cd $WD
fi    

if [ $re_map -eq 1 ]
    then 
        echo 'Refining mapping and calling consensus sequence'
        refnam=$(sort -rnk3 $outdir/mapping_info | head -1 | cut -f1)
        grep -Pzo '(?s)>'$refnam'.*?>' $outdir/ref_nogap.fas |head -n-1 > $outdir/best_ref_uneven.fas
        grep -Pzo '(?s)>'$align'.*?>' $outdir/ref_gaps.fas 
        for i in $(grep -aob '-' hand_best_ref_gaps.fas|sort -rn | cut -d":" -f1);
            do sed -i "s/./&-/$i" ./cns_aln.fa;
        done

        python fastafixer.py $outdir/best_ref_uneven.fas $outdir/best_ref.fas
        echo 'The best reference found in your alignment was '$refnam
        echo 'mapping reads to '$refnam

        bowtie2-build $outdir/best_ref.fas $outdir/best_ref >> $outdir/bowtiebuild.log

        #TOTDO THINK HARD ABOUT IMPLAICTIONS OF LOCAL VS GLOBAL AIGN!!!
        if [ $PE -eq 1 ]
        	then 
        	    bowtie2 -x $outdir/best_ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/best_map.sam --no-unal --local
            else 
            	bowtie2 -x $outdir/best_ref  -U ${read_stub}.fastq -S $outdir/best_map.sam --no-unal --local
        fi

        samtools faidx $outdir/best_ref.fas
        samtools view -bS $outdir/best_map.sam > $outdir/best_map.bam
        samtools sort $outdir/best_map.bam -o $outdir/best_sorted.bam
        samtools index $outdir/best_sorted.bam 
        samtools mpileup -uf $outdir/best_ref.fas $outdir/best_sorted.bam| bcftools call -c | vcfutils.pl vcf2fq >  $outdir/cns.fq
        seqtk seq -a $outdir/cns.fq > $outdir/cns.fa
        sed -i -e "s/>/>${nam}_/g" $outdir/cns.fa
        if [ $wre_map -eq 1 ]
            then
                #generate consensus mapped to second best locus to assuage some ref dependence
                wrefnam=$(sort -rnk3 $outdir/mapping_info | head -2 | tail -1| cut -f1)
                grep -Pzo '(?s)>'$wrefnam'.*?>' $outdir/ref_nogap.fas |head -n-1 > $outdir/worse_ref_uneven.fas
                python fastafixer.py $outdir/worse_ref_uneven.fas $outdir/worse_ref.fas
                echo 'The second best reference found in your alignment was '$wrefnam
                echo 'mapping reads to '$wrefnam

                bowtie2-build $outdir/worse_ref.fas $outdir/worse_ref >> $outdir/bowtiebuild.log

                if [ $PE -eq 1 ];
                    then 
                        bowtie2 -x $outdir/worse_ref -1 ${read_stub}_1.fastq -2 ${read_stub}_2.fastq -S $outdir/worse_map.sam --no-unal --local
                    else 
                        bowtie2 -x $outdir/worse_ref  -U ${read_stub}.fastq -S $outdir/worse_map.sam --no-unal --local
                fi
                samtools faidx $outdir/worse_ref.fas 
                samtools view -bS $outdir/worse_map.sam > $outdir/worse_map.bam
                samtools sort $outdir/worse_map.bam -o $outdir/worse_sorted.bam
                samtools index $outdir/worse_sorted.bam 
                samtools mpileup -uf $outdir/worse_ref.fas $outdir/worse_sorted.bam| bcftools call -c | vcfutils.pl vcf2fq >  $outdir/worse_cns.fq 
                seqtk seq -a $outdir/worse_cns.fq > $outdir/worse_cns.fa
                sed -i -e "s/>/>${nam}_worse_/g" $outdir/worse_cns.fa
                if [ $(diff $outdir/cns.fa $outdir/worse_cns.fa | wc -l | cut -f3) -gt 4 ]
                     then 
                        echo 'Alternate references result in different sequences. Placing both, but investigating differences recommended!'
                        cat $outdir/cns.fa $outdir/worse_cns.fa > $outdir/mappings.fa
                        cd $outdir
                          papara -t ${WD}/${tree} -s ${aln_stub}.phy -q mappings.fa -n multi_consensus 
                          #run RAXML EPA on the alignments
                          raxmlHPC -m GTRCAT -f v -s papara_alignment.multi_consensus -t ${WD}/$tree -n ${nam}_consensusPC
                        cd $WD
                else
                    echo 'Using worse reference resulted in identical sequences - only aligning and placing one.'
                    cd $outdir
                        papara -t ${WD}/${tree} -s ${aln_stub}.phy -q cns.fa -n re_consensus 
                        raxmlHPC -m GTRCAT -f v -s papara_alignment.re_consensus -t ${WD}/$tree -n ${nam}_consensusPC
                    cd $WD
                fi
        else #when is this condition met?
            cd $outdir
              #  papara -t ${WD}/${tree} -s ${aln_stub}.phy -q cns.fa -n fi_consensus 
                mafft --add cns.fa --reorder ${WD}/$align > extended.aln
                raxmlHPC -m GTRCAT -f v -s extended.aln -t ${WD}/$tree -n ${nam}_consensusPC
            cd $WD
        fi
        #run full raxml? tooo sloooo
        raxmlHPC -m GTRGAMMA -s $outdir/contig_alignment.fas -t $tree -p 12345 -n consensusFULL

fi

#todo strip all fq to fa 

 